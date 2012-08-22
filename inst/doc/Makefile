whoami=$(shell whoami)

# Makefile for package vignettes

# put all PDF targets here, separated by spaces
PDFS= mRMRe.pdf

all: $(PDFS) 

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi *.spl *tikzDictionary

%.pdf: %.Rnw
	"$(R_HOME)/bin/R" CMD Sweave $*
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*
