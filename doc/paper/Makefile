all:
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/latex.template
	pandoc \
    -V repository="https://github.com/pencil-code/pencil-code" \
    -V archive_doi="http://dx.doi.org/10.5281/zenodo.3961647" \
    -V graphics="true" \
	-V journal_name="Journal of Open Source Software" \
    -V issue="34" \
    -V volume="4" \
    -V page="1265" \
    -V logo_path="joss-logo.png" \
    -V year="2020" \
    -V submitted="17 September 2020" \
    -V citation_author="The Pencil Code Collaboration" \
    -V paper_title="The Pencil Code, a modular MPI code for partial differential equations and particles: multipurpose and multiuser-maintained" \
    paper.md --template latex.template --filter pandoc-citeproc --bibliography paper.bib -o paper.pdf -V geometry:margin=1in

tex:
	wget -N https://raw.githubusercontent.com/openjournals/whedon/master/resources/latex.template
	pandoc \
    -V repository="https://github.com/pencil-code/pencil-code" \
    -V archive_doi="http://dx.doi.org/10.5281/zenodo.3961647" \
    -V graphics="true" \
	-V journal_name="Journal of Open Source Software" \
    -V issue="submitted" \
    -V page="1265" \
    -V logo_path="joss-logo.png" \
    -V year="2020" \
    -V submitted="17 September 2020" \
    -V citation_author="The Pencil Code Collaboration" \
    -V paper_title="The Pencil Code, a modular MPI code for partial differential equations and particles: multipurpose and multiuser-maintained" \
    paper.md --template latex.template --filter pandoc-citeproc --bibliography paper.bib -o paper.tex -V geometry:margin=1in
