`set.thread.count` <- function(thread_count)
{
    .Call(C_set_thread_count, as.integer(thread_count))
}