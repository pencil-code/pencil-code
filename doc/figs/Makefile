# Makefile for transforming EPS pictures to PDF


.SUFFIXES: .eps .pdf

.eps.pdf:
	epstopdf $?

# ------------------

default: pdf


pdf: *.eps
	epstopdf_all $?


# Clean up
clean:
	rm -f *.pdf

# End of file

