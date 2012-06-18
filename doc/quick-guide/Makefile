default: compile

compile:  quick_start.tex
	latex quick_start
	latex quick_start
	dvips -R0 quick_start.dvi -o quick_start.ps
	ps2pdf quick_start.ps quick_start.pdf
	#xpdf hello_disk.pdf &
clean:
	rm -f *.log *.aux *.dvi *.ps *~
