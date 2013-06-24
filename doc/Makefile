# Makefile for generating different output formats of the manual.
#
# Usage:
#  make
#  make manual.ps
#  make manual.pdf
#  make body.ps      (excludes appendix and indexes)
#  make appendix.ps  (appendix and indexes only)
#  make faq          (html version of Frequently Asked Questions)
#
# Notes:
#  - The poor buggers under Solaris should use `gmake' instead of
#    `make', as the latter behaves completely crazy.
#  - If dvipdfm does a bad job for you (e.g. omits all graphics),
#    run `make manual_pdf' instead of `make manual.pdf'

LATEX = latex -src
DVIPS = dvips -Ppdf -G0
PS2PDF = ps2pdf14 -sPAPERSIZE=a4
DVIPDFM = dvipdfm -p a4

MAY_HAVE_CHANGED='may have changed. Rerun to get'

.SUFFIXES:            # clear all generic rules and suffixes
.SUFFIXES: .tex .dvi .ps .pdf

.tex.dvi:
	$(LATEX) $*
	if [ `fgrep -c $(MAY_HAVE_CHANGED) $*.log` != '' ]; then \
	    $(LATEX) $*; \
	fi
.dvi.ps:
	$(DVIPS) $* -o $*.ps
.ps.pdf:
	$(PS2PDF) $*.ps $*.pdf
#.dvi.pdf:
#	$(DVIPDFM) $*.dvi $*.pdf

# ------------------

default: inlinedoc manual.pdf manual.ps.gz

all:  manual.pdf manual.dvi manual.ps.gz

pdffigs:
	cd figs; make pdf

manual_dvi: manual.tex inlinedoc-diag.tex
	@TEXINPUTS="./texinputs:${TEXINPUTS}"; \
	export TEXINPUTS; \
	echo "TEXINPUTS=$${TEXINPUTS}"; \
	$(LATEX) manual \
	    && makeindex -s dotted.idxsty manual \
	    && makeindex -s dotted.idxsty manual.vidx -o manual.vind \
	    && makeindex -s dotted.idxsty manual.fidx -o manual.find \
	    && $(LATEX) manual \
	    && if [ `fgrep -c $(MAY_HAVE_CHANGED) manual.log` != '' ]; then \
	           $(LATEX) manual; \
	       fi

manual.dvi: driver_dvi manual_dvi
	@true  # dummy cmd, so make won't apply generic rule at end

# By now (Aug 2007), on Ubuntu and RedHat the route through .ps seems to
# be the only way to get hyperlinks.
manual.pdf: driver_ps2pdf manual.ps
	$(PS2PDF) manual.ps

# This is what we used previously: dvipdfm (works fine on Gentoo, but
# probably not much else)
manual_pdf: driver_ manual_dvi		# don't use driver_dvipdfm...
	dvipdfm -p a4 manual.dvi

manual.ps: driver_ps2pdf manual_dvi
	$(DVIPS) manual -o manual.ps

manual.ps.gz: manual.ps
	gzip -f manual.ps

inlinedoc-diag.tex: inlinedoc

inlinedoc:
	./extract-diag-doc

driver_:
	@printf "%s\n" \
	    '% Auto-generated file' \
	    '\def\mydriver{}' \
	    > driver_switch.tex
driver_dvi:
	@printf "%s\n" \
	    '% Auto-generated file' \
	    '\def\mydriver{dvips}' \
	    > driver_switch.tex
driver_ps2pdf:
	@printf "%s\n" \
	    '% Auto-generated file' \
	    '\def\mydriver{ps2pdf}' \
	    > driver_switch.tex
driver_dvips:
	@printf "%s\n" \
	    '% Auto-generated file' \
	    '\def\mydriver{dvips}' \
	    > driver_switch.tex
driver_dvipdfm:
	@printf "%s\n" \
	    '% Auto-generated file' \
	    '\def\mydriver{dvipdfm}' \
	    > driver_switch.tex


body.ps: manual.tex
	$(MAKE) manual.dvi
	last=`perl -000 -ne 'print $$1 if (/=======LAST_BODY_PAGE=======\s*\[([0-9]+)\]/s)' manual.log`; \
	dvips -l "$$last" manual.dvi -o body.ps
	rm -f manual.dvi	# or it might get used by next `make manual.pdf'

appendix.ps: manual.tex
	$(MAKE) manual.dvi
	first=`perl -000 -ne 'print $$1+1 if (/=======LAST_BODY_PAGE=======\s*\[([0-9]+)\]/s)' manual.log`; \
	echo "first = <$$first>"; \
	dvips -p "$$first" manual.dvi -o appendix.ps
	rm -f manual.dvi	# or it might get used by next `make manual.pdf'

faq: www/faq.html

# FAQ as HTML
# Note: If you get
#   --- warning --- Can't find/open file tex4ht.env | .tex4ht
#   --- error --- Illegal storage address
# try
#   sudo texhash
# [didn't help for me, texlive-2012.20120611-3~ubuntu12.04.1], or
#   ln -s /etc/tex4ht/tex4ht.env ~/.tex4ht

www/faq.html: www/faq.tex
	cd www; \
	if [ -z "`type htlatex 2> /dev/null`" ]; then \
		htlatex="htlatex.sh"; \
	else \
		htlatex="htlatex"; \
	fi; \
	env TEXINPUTS="$${TEXINPUTS}:.." $$htlatex faq


www/faq.tex: manual.tex
# 1. Extract header
# 2. Extract \section{... FAQ} till next section/part
# 3. Comment out \usepackage{titlesec} (which has problems with htlatex)
# 4. replace last \section/\part with \end{document}
	if [ ! -d www ]; then mkdir www; fi
	perl -n -e 'print if (1 .. /^\s*\\begin{document}/);' \
	        -e 'print if (/^\s*\\section.*Frequently Asked Questions/ ... /^\s*\\(section|part)/)' manual.tex \
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
# same, but including image directory and automatically generated inline
# doc files
cleanall: cleann
cleann: clean
	rm -f *.pdf inlinedoc-*
	cd figs; make clean

# End of file
