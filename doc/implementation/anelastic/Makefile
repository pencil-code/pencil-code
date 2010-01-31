## $Id$

TARGET=anelastic_hurdle
MAY_HAVE_CHANGED='Rerun to get .* \(right\|correct\)'
LATEX=latex -src
PDFLATEX=pdflatex

## Generic rules:
.SUFFIXES: .tex .pdf .dvi .ps

.tex.dvi:
	$(LATEX) $* || (rm -f $*.dvi; false)
	if ( grep -q $(MAY_HAVE_CHANGED) $*.log ); then $(LATEX) $*; fi
	if ( grep -q $(MAY_HAVE_CHANGED) $*.log ); then $(LATEX) $*; fi

.tex.pdf:
	$(PDFLATEX) $* || (rm -f $*.pdf; false)
	if ( grep -q $(MAY_HAVE_CHANGED) $*.log ); then $(PDFLATEX) $*; fi
	if ( grep -q $(MAY_HAVE_CHANGED) $*.log ); then $(PDFLATEX) $*; fi

.dvi.ps:
	dvips $*

## Specific rules
all: $(TARGET).pdf

$(TARGET).pdf: *.tex
$(TARGET).dvi: *.tex
$(TARGET).ps: $(TARGET).dvi

clean:
	bash -c "rm -f *.{aux,log,out,toc,ilg,pdf} comment.cut"

cleann: clean
	bash -c "rm -f $(TARGET).{idx,blg,bbl,dvi,pdf,ps,ind,ilg,aux,toc,out,idx,log}"

cleanall: cleann


# End of file
