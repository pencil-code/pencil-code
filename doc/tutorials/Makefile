default: hello

hello:  hello_disk.tex
	latex hello_disk
	latex hello_disk
	dvips -o hello_disk.ps hello_disk.dvi
	ps2pdf hello_disk.ps hello_disk.pdf
	xpdf hello_disk.pdf &
clean:
	rm -f *.log *.aux *.dvi *.ps *.pdf
