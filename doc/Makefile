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


pdffigs:
	cd figs; make pdf

manual.dvi: manual.tex
	if [ -z $TEXINPUTS ]; then TEXINPUTS=":./texinputs"; else TEXINPUTS="$(TEXINPUTS):./texinputs"; fi; export TEXINPUTS; latex manual; makeindex manual; latex manual

manual.pdf: manual.tex pdffigs
	if [ -z $TEXINPUTS ]; then TEXINPUTS=":./texinputs"; else TEXINPUTS="$(TEXINPUTS):./texinputs"; fi; export TEXINPUTS; pdflatex manual; makeindex manual; pdflatex manual

# Clean up
clean:
	rm -f manual.ps manual.dvi *.aux *.toc *.log *.idx *.ind *.ilg
# same, but including image directory
cleann: clean
	cd figs; make clean

# End of file
