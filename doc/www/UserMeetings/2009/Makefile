default: prop

funding_proposal.pdf: funding_proposal.tex
	pdflatex funding_proposal.tex
	pdflatex funding_proposal.tex

prop: funding_proposal_dfg.tex
	pdflatex funding_proposal_dfg.tex
	pdflatex funding_proposal_dfg.tex
	xpdf funding_proposal_dfg.pdf &

clean: 
	rm -f *.log *.aux *.dvi *.ps *.pdf *out
