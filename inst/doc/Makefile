whoami=$(shell whoami)

# Makefile for package vignettes

# put all PDF targets here, separated by spaces
PDFS= mRMRe.pdf

all: $(PDFS) 

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.tikz *.toc *.map *.sh *.dvi *.bib *.spl *tikzDictionary

%.pdf: %.Rnw
	$(R_HOME)/bin/Rscript -e "Sweave('$*.Rnw')"
	pdflatex mRMRe
	bibtex mRMRe
	pdflatex mRMRe
	pdflatex mRMRe

