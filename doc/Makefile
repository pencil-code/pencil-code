# Makefile for generating different output formats of the manual.
#
# Usage:
#  make
#  make manual.ps
#  make body.ps      (excludes appendix and indexes)
#  make appendix.ps  (appendix and indexes only)
#  make faq          (html version of Frequently Asked Questions)

.SUFFIXES: .tex .dvi .ps .pdf

.dvi.ps:
	dvips $* -o $*.ps
.ps.pdf:
	ps2pdf $*.ps $*.pdf
#.dvi.pdf:
#	dvipdf $*.dvi $*.pdf

# ------------------

default: manual.dvi

all: manual.dvi manual.ps.gz manual.pdf

pdffigs:
	cd figs; make pdf

manual.dvi: manual.tex
	if [ -z $TEXINPUTS ]; then TEXINPUTS="./texinputs:"; else TEXINPUTS="./texinputs:$(TEXINPUTS)"; fi; export TEXINPUTS; latex manual; makeindex -s dotted.idxsty manual; makeindex -s dotted.idxsty manual.vidx -o manual.vind; makeindex -s dotted.idxsty manual.fidx -o manual.find; latex manual

# manual.pdf: manual.tex pdffigs
# 	if [ -z $TEXINPUTS ]; then TEXINPUTS=":./texinputs"; else TEXINPUTS="$(TEXINPUTS):./texinputs"; fi; export TEXINPUTS; pdflatex manual; makeindex manual; pdflatex manual
manual.pdf: manual.dvi
	dvipdfm -p a4 manual.dvi

# To generate manual.ps, we need to temporarily replace
# `\usepackage[dvipdfm]{hyperref' by `\usepacakge{hyperref}':
manual.ps: manual.tex
	perl -i.bak1 -0777 -pe \
	  's/^(\s*\\usepackage)\[dvipdfm\](\{hyperref\})/$$1$$2/m' manual.tex
	$(MAKE) manual.dvi
	perl -i.bak2 -0777 -pe \
	  's/^(\s*\\usepackage)(\{hyperref\})/$$1\[dvipdfm\]$$2/m' manual.tex
	dvips -Ppdf -G0 manual -o manual.ps
	rm -f manual.dvi	# or it might get used by next `make manual.pdf'

manual.ps.gz: manual.ps
	gzip -f manual.ps

body.ps: manual.dvi
	last=`perl -000 -ne 'print $$1 if (/=======LAST_BODY_PAGE=======\s*\[([0-9]+)\]/s)' manual.log`; \
	dvips -l "$$last" manual.dvi -o body.ps

appendix.ps: manual.dvi
	first=`perl -000 -ne 'print $$1+1 if (/=======LAST_BODY_PAGE=======\s*\[([0-9]+)\]/s)' manual.log`; \
	echo "first = <$$first>"; \
	dvips -p "$$first" manual.dvi -o appendix.ps

faq: www/faq.html

www/faq.html: www/faq.tex
	(cd www; htlatex faq)

www/faq.tex: manual.tex
	if [ ! -d www ]; then mkdir www; fi
	perl -ne 'print if (1 .. /^\s*\\begin{document}/); \
	          print if (/^\s*\\section{Frequently Asked Questions}/ ... \
	                    /^\s*\\section/)' manual.tex \
	| tac | perl -0777 -pe's/^\s*\\section.*/\\end{document}\n/' | tac > www/faq.tex

# Clean up
clean:
	rm -f manual.ps* body.ps* appendix.ps* manual.dvi *.aux *.toc *.log \
	*.idx *.ind *.ilg *.vidx *.vind *.fidx *.find
# same, but including image directory
cleann: clean
	cd figs; make clean

# End of file
