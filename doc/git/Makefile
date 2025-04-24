default: git-best-practises.pdf git-getting-started.pdf

git-getting-started.pdf: git-getting-started.tex
	pdflatex git-getting-started
	pdflatex git-getting-started

git-best-practises.pdf: git-best-practises.tex git-better-history.pdf git-wrong-history.pdf
	pdflatex git-best-practises
	pdflatex git-best-practises

# Clean up
clean:
	*.aux *.toc *.log \
	*.idx *.ind *.ilg *.out *.vidx *.vind *.fidx *.find
# same, but including image directory and automatically generated inline
# doc files
cleanall: cleann
cleann: clean
	rm -f *.pdf inlinedoc-*
	cd figs; make clean
