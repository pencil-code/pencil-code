# Makefile for generating different output formats of the manual.


.SUFFIXES: .tex .dvi .ps .pdf

.dvi.ps:
	dvips $* -o $*.ps
.ps.pdf:
	ps2pdf $*.ps $*.pdf
#.dvi.pdf:
#	dvipdf $*.dvi $*.pdf

# ------------------

default: manual.dvi

manual.dvi: manual.tex manual.ind
	if [ -z $TEXINPUTS ]; then TEXINPUTS=":./texinputs"; else TEXINPUTS="$(TEXINPUTS):./texinputs"; fi; export TEXINPUTS; latex manual; makeindex manual; latex manual


# End of file
