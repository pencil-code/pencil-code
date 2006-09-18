.SUFFIXES: .tex .eps .pdf .png .ps

PDF_FIGS=

.ps.pdf:
	epstopdf $*.ps

.eps.pdf:
	epstopdf $*.eps

.png.ps:
	convert $*.ps $*.png

default: link_pngs pdffigs

link_pngs:
	for file in ../png/*.png; do ln -sf $$file .; done

pdffigs: $(PDF_FIGS) 

clean:
	rm -f $(PDF_FIGS)

img_0009.ps: img_0019.png
img_0024.ps: img_0024.png
img_0049.ps: img_0049.png
img_0498.ps: img_0498.png
pcomp_spec1.pdf: pcomp_spec1.ps
pcomp_spec2.pdf: pcomp_spec2.ps
pcomp_spec1dt.pdf: pcomp_spec1dt.ps
pcomp_spec1dt_nu.pdf: pcomp_spec1dt_nu.ps
pcomp_spec1dt_res.pdf: pcomp_spec1dt_res.ps
pspec512b3.pdf: pspec512b3.ps

