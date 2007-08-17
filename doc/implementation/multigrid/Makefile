# Makefile

## $Revision: 1.1 $
## $Date: 2007-08-17 04:40:43 $

TARGET=multigrid
MAY_HAVE_CHANGED='Rerun to get .* \(right\|correct\)'
LATEX=latex -src
PDFLATEX=pdflatex

## Generic rules:
.SUFFIXES: .tex .pdf .dvi .ps

.tex.dvi:
	$(LATEX) $* || (rm -f $*.dvi; false)
	makeindex $*
	if ( grep -q $(MAY_HAVE_CHANGED) $*.log ); then $(LATEX) $*; fi
	if ( grep -q $(MAY_HAVE_CHANGED) $*.log ); then $(LATEX) $*; fi

.tex.pdf:
	$(PDFLATEX) $* || (rm -f $*.pdf; false)
	makeindex $*
	if ( grep -q $(MAY_HAVE_CHANGED) $*.log ); then $(PDFLATEX) $*; fi
	if ( grep -q $(MAY_HAVE_CHANGED) $*.log ); then $(PDFLATEX) $*; fi

.dvi.ps:
	dvips $*

## Specific rules
all: $(TARGET).pdf

$(TARGET).pdf: *.tex
$(TARGET).dvi: *.tex
$(TARGET).ps: $(TARGET).dvi

test.dvi: *.tex
test.pdf: *.tex
test.ps: test.dvi

clean:
	rm -f *.{aux,log,out,toc,ilg} comment.cut

prv_clean:
	rm -rf *.prv _region_.*

cleann: clean prv_clean
	bash -c "rm -f  $(TARGET).{idx,blg,bbl,dvi,pdf,ps,ind,ilg,aux,toc,out,idx,log}"
cleanall: cleann


# End of file
