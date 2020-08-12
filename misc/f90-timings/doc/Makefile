.SUFFIXES: .tex .dvi .ps .pdf .ps.gz

.ps.ps.gz:
	gzip -c $*.ps > $*.ps.gz

.dvi.ps:
	dvips $* -o $*.ps

default: timings.ps.gz timings.pdf

timings.ps.gz: timings.ps
timings.ps: timings.dvi

PDFFIGS=
PSFFIGS=

TEXSOURCES=timings.tex

timings.dvi: $(TEXSOURCES) $(PSFIGS)
	make -C fig
	( latex --jobname timings_dvi --interaction nonstopmode timings >& /dev/null ; true) 
	( latex --jobname timings_dvi --interaction nonstopmode timings >& /dev/null ; true) 
	latex --jobname timings_dvi timings
	mv timings_dvi.dvi timings.dvi

timings.pdf: $(TEXSOURCES) $(PDFFIGS)
	make -C fig pdffigs
	( pdflatex --jobname timings_pdf --interaction nonstopmode timings >& /dev/null ; true)
	( pdflatex --jobname timings_pdf --interaction nonstopmode timings >& /dev/null ; true)
	pdflatex --jobname timings_pdf timings
	mv timings_pdf.pdf timings.pdf

submit.tar.gz: $(TEXSOURCES) $(PSFIGS)
	tar zcvf submit.tar.gz timings.tex \
             fig/img_0009.ps fig/pcomp_spec1.ps fig/pcomp_spec2.ps \
             fig/img_0024.ps fig/pcomp_spec1dt.ps  \
             fig/img_0049.ps fig/pcomp_spec1dt_nu.ps fig/pspec512b3.ps \
             fig/img_0498.ps fig/pcomp_spec1dt_res.ps 

clean:
	make -C fig clean
	rm -f timings.pdf timings.ps timings.ps.gz *.aux *.log *.dvi submit.tar.gz
