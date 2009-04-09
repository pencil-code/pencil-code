default: hello

hello:  hello_disk.tex
	latex hello_disk
	latex hello_disk
	dvips -R0 hello_disk.dvi -o hello_disk.ps
	ps2pdf hello_disk.ps hello_disk.pdf
	#xpdf hello_disk.pdf &
clean:
	rm -f *.log *.aux *.dvi *.ps *~
