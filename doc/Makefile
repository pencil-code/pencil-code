# Makefile for generating different output formats of the manual.
#
# Usage:
#  make
#  make manual.ps
#  make manual.pdf
#  make body.ps      (excludes appendix and indexes)
#  make appendix.ps  (appendix and indexes only)
#  make faq          (html version of Frequently Asked Questions)

LATEX = latex -src
MAY_HAVE_CHANGED='may have changed. Rerun to get'

.SUFFIXES:            # clear all generic rules and suffixes
.SUFFIXES: .tex .dvi .ps .pdf

.tex.dvi:
	dvips $* -o $*.ps
.dvi.ps:
	dvips $* -o $*.ps
.ps.pdf:
	ps2pdf $*.ps $*.pdf
#.dvi.pdf:
#	dvipdfm $*.dvi $*.pdf

# ------------------

default: manual.dvi

all: manual.dvi manual.ps.gz manual.pdf

pdffigs:
	cd figs; make pdf

manual.dvi: manual.tex
	@if [ -z $${TEXINPUTS} ]; then \
	    TEXINPUTS="./texinputs:"; \
	else \
	    TEXINPUTS="./texinputs:$(TEXINPUTS)"; \
	fi; \
	export TEXINPUTS; \
	echo "TEXINPUTS=$${TEXINPUTS}"; \
	$(LATEX) manual \
	    && makeindex -s dotted.idxsty manual \
	    && makeindex -s dotted.idxsty manual.vidx -o manual.vind \
	    && makeindex -s dotted.idxsty manual.fidx -o manual.find \
	    && $(LATEX) manual \
	    && if ( fgrep -q $(MAY_HAVE_CHANGED) manual.log ); then \
	           latex $*; \
	       fi

# manual.pdf: manual.tex pdffigs
# 	if [ -z $TEXINPUTS ]; then TEXINPUTS=":./texinputs"; else TEXINPUTS="$(TEXINPUTS):./texinputs"; fi; export TEXINPUTS; pdflatex manual; makeindex manual; pdflatex manual
manual.pdf: manual.dvi
	dvipdfm -p a4 manual.dvi
	rm -f manual.dvi

# To generate manual.ps, we need to temporarily replace
# `\usepackage[dvipdfm]{hyperref' by `\usepacakge{hyperref}':
# [A really ugly hack -- this should be done using am dynamic include file]
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

# To generate body/appendix.ps, we need to temporarily replace
# `\usepackage[dvipdfm]{hyperref' by `\usepacakge{hyperref}':
body.ps: manual.tex
	perl -i.bak1 -0777 -pe \
	  's/^(\s*\\usepackage)\[dvipdfm\](\{hyperref\})/$$1$$2/m' manual.tex
	$(MAKE) manual.dvi
	perl -i.bak2 -0777 -pe \
	  's/^(\s*\\usepackage)(\{hyperref\})/$$1\[dvipdfm\]$$2/m' manual.tex
	last=`perl -000 -ne 'print $$1 if (/=======LAST_BODY_PAGE=======\s*\[([0-9]+)\]/s)' manual.log`; \
	dvips -l "$$last" manual.dvi -o body.ps
	rm -f manual.dvi	# or it might get used by next `make manual.pdf'

appendix.ps: manual.tex
	perl -i.bak1 -0777 -pe \
	  's/^(\s*\\usepackage)\[dvipdfm\](\{hyperref\})/$$1$$2/m' manual.tex
	$(MAKE) manual.dvi
	perl -i.bak2 -0777 -pe \
	  's/^(\s*\\usepackage)(\{hyperref\})/$$1\[dvipdfm\]$$2/m' manual.tex
	first=`perl -000 -ne 'print $$1+1 if (/=======LAST_BODY_PAGE=======\s*\[([0-9]+)\]/s)' manual.log`; \
	echo "first = <$$first>"; \
	dvips -p "$$first" manual.dvi -o appendix.ps
	rm -f manual.dvi	# or it might get used by next `make manual.pdf'

faq: www/faq.html

www/faq.html: www/faq.tex
	cd www; if [ -z "`type htlatex 2> /dev/null`" ]; \
	        then htlatex="htlatex.sh"; \
	        else htlatex="htlatex"; fi; $$htlatex faq

www/faq.tex: manual.tex
# 1. Extract header
# 2. Extract \section{... FAQ} till next section/part
# 3. Comment out \usepackage{titlesec} (which has problems with htlatex)
# 4. replace last \section/\part with \end{document}
	if [ ! -d www ]; then mkdir www; fi
	perl -n -e 'print if (1 .. /^\s*\\begin{document}/);' \
	        -e 'print if (/^\s*\\section.*Frequently Asked Questions/ \
	                      ... \
	                      /^\s*\\(section|part)/)' manual.tex \
	| perl -p -e 's/^(\\usepackage\[[^\]]*]\s*\{titlesec\b)/%% $$1/;' \
	          -e 's/^(\s*\\titleformat)/%% $$1/' \
	| tac \
	| perl -0777 -p -e 's/^\s*\\(section|part).*/\\end{document}\n/' \
	| tac \
	> www/faq.tex

# Clean up
clean:
	rm -f manual.ps* body.ps* appendix.ps* manual.dvi *.aux *.toc *.log \
	*.idx *.ind *.ilg *.out *.vidx *.vind *.fidx *.find manual.image.*
# same, but including image directory
cleann: clean
	rm -f *.pdf
	cd figs; make clean

# End of file
