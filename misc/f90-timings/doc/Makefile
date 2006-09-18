.SUFFIXES: .tex .dvi .ps .pdf .ps.gz

.ps.ps.gz:
	gzip -c $*.ps > $*.ps.gz

.dvi.ps:
	dvips $* -o $*.ps

default: paper.ps.gz paper.pdf

paper.ps.gz: paper.ps
paper.ps: paper.dvi

PDFFIGS=
PSFFIGS=

TEXSOURCES=paper.tex

paper.dvi: $(TEXSOURCES) $(PSFIGS)
	make -C fig
	( latex --jobname paper_dvi --interaction nonstopmode paper >& /dev/null ; true) 
	( latex --jobname paper_dvi --interaction nonstopmode paper >& /dev/null ; true) 
	latex --jobname paper_dvi paper
	mv paper_dvi.dvi paper.dvi

paper.pdf: $(TEXSOURCES) $(PDFFIGS)
	make -C fig pdffigs
	( pdflatex --jobname paper_pdf --interaction nonstopmode paper >& /dev/null ; true)
	( pdflatex --jobname paper_pdf --interaction nonstopmode paper >& /dev/null ; true)
	pdflatex --jobname paper_pdf paper
	mv paper_pdf.pdf paper.pdf

submit.tar.gz: $(TEXSOURCES) $(PSFIGS)
	tar zcvf submit.tar.gz paper.tex \
             fig/img_0009.ps fig/pcomp_spec1.ps fig/pcomp_spec2.ps \
             fig/img_0024.ps fig/pcomp_spec1dt.ps  \
             fig/img_0049.ps fig/pcomp_spec1dt_nu.ps fig/pspec512b3.ps \
             fig/img_0498.ps fig/pcomp_spec1dt_res.ps 

clean:
	make -C fig clean
	rm -f paper.pdf paper.ps paper.ps.gz *.aux *.log *.dvi submit.tar.gz
