library(Rmpi)

slave_function <- function()
{
    fate <- integer(1)
    mpi.recv(x = fate, type = 1, source = 0, tag = 0)
    print(paste("received fate ", fate, mpi.comm.rank() ))
    if (as.logical(fate))
        active_slave_function()
    else
        passive_slave_function()
}

passive_slave_function <- function()
{
    run <- TRUE
    print("In passive")
    while (run)
    {
        Sys.sleep(1)
        #mpi.irecv(x = junk, type = 1, source = 0, tag = 0, request = 0)
        
        #if (mpi.test(request = 0) && mpi.get.sourcetag()[2] == 0)
        #    run <<- FALSE;
    }
    
    #mpi.send(x = 1, type = 1, dest = 0, tag = 0)
}

active_slave_function <- function()
{
    print("In active")
	
	mpi.recv(x = task, type = 1, source = 0, tag = 0)
	if(task == 0)
	{
		nb_threads <- 1
		mpi.recv(x = nb_threads, type = 1, source = 0, tag = 0)
		
		
		data_matrix <- mpi.recv.Robj(source = 0, tag = 0)
		sapply(1:ncol(data_matrix), function(i) {
			   sapply((i + 1):ncol(data_matrix), function(j) {
					  cor(data_matrix[ , i], data_matrix[ , j])
					  })
			   })
		cor_matrix <- cor(data_matrix)
		mpi.send.Robj(x = cor_matrix, dest = 0, tag = 0)
	}
    Sys.sleep(5)
    
    mpi.send(x = 0, type = 1, dest = 0, tag = 0)
    
    mpi.send.Robj(obj = paste("active: ", mpi.comm.rank()), dest = 0, tag = 0)
}

master_function <- function()
{
    mpi.spawn.Rslaves(nslaves=8)
    mpi.remote.exec(library(Rmpi))
    all_slaves <- mpi.remote.exec(list(rank = mpi.comm.rank(), host = mpi.get.processor.name()))
    mpi.remote.exec(source("~/mpi.R"))
    active_slaves <- list()
    passive_slaves <- list()
    
    mpi.bcast.cmd(slave_function())
    
    lapply(all_slaves, function(slave)
    {
        if (length(active_slaves) > 0 &&
                sum(sapply(active_slaves, function(i) slave$host == i$host)) > 1)
        {
	    print(paste("changing ", slave$rank,"to passive"))
            passive_slaves[[length(passive_slaves) + 1]] <<- slave
            mpi.send(x = 0, type = 1, dest = slave$rank, tag = 0)
        }
        else
        {
	    print(paste("changing ", slave$rank,"to active"))
            active_slaves[[length(active_slaves) + 1]] <<- slave
            mpi.send(x = 1, type = 1, dest = slave$rank, tag = 0)
        }
    })
    
    active_slaves_on_wait <- rep(TRUE, length(active_slaves))
    first_run <- TRUE 
    while (sum(active_slaves_on_wait))
    {
        Sys.sleep(1)
        
        lapply(seq(length(active_slaves_on_wait)), function(active_slave_index)
        {
            slave_status <- integer(1)
            
            if (active_slaves_on_wait[[active_slave_index]])
            {
		#browser()
		print(paste("waiting for ", active_slaves[[active_slave_index]]$rank))
                if(first_run)
			mpi.irecv(x = slave_status, type = 1, source = active_slaves[[active_slave_index]]$rank, tag = 0, request =active_slaves[[active_slave_index]]$rank)
                if (mpi.test(request = active_slaves[[active_slave_index]]$rank))
                {
                    active_slaves_on_wait[[active_slave_index]] <<- FALSE
                    print(paste("received aknowledgement from ", active_slave_index))
                }
                else
                  print( "Pas recu")
	    }
        })
	first_run <- FALSE
    }
	test_matrix <- matrix(runif(1000**2), 1000, 1000)
    lapply(active_slaves, function(active_slave)
    {
		   mpi.send(x = 0, type = 1, dest = active_slave$rank, tag = 0)
		   mpi.send(x = length(active_slaves), type = 1, dest = active_slave$rank, tag = 0)
		   mpi.send.Robj(test_matrix, active_slave$rank, tag = 0)
    })
	
	lapply(active_slaves, function(active_slave)
	{
	     mpi.recv.Robj(source = active_slave$rank, tag = 0)
		
	})

    #mpi.close.Rslaves()
    #mpi.finalize()
}

