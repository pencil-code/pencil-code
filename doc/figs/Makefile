# Makefile for transforming pictures to


.SUFFIXES: .eps .pdf

.eps.pdf:
	epstopdf $?
#.dvi.pdf:
#	dvipdf $*.dvi $*.pdf

epsfiles = $(wildcard *.eps)

# ------------------

default: pdf


pdf: $(epsfiles:.eps=.pdf)

# Clean up
clean:
	rm -f *.pdf

# End of file

