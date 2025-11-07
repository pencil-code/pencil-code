FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y gcc gfortran pip openmpi-bin libopenmpi-dev git libxrender1 x11-apps python3-tk
RUN pip3 install numpy sphinx-autoapi sphinxcontrib.images sphinxcontrib.bibtex sphinx-rtd-theme sphinx rstcloth six bibtexparser latexcodec matplotlib
WORKDIR /app/pencil-code/doc/readthedocs
RUN  git config --global --add safe.directory /app/pencil-code
