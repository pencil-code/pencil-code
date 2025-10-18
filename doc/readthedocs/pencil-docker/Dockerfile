FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y gcc gfortran pip openmpi-bin libopenmpi-dev git libxrender1 x11-apps python3-tk gnudatalanguage python3-gdl vim
RUN pip3 install numpy scipy matplotlib setuptools==68.2.2 eqtools==1.4.0 h5py plotly retrying astropy Cython dill pexpect tqdm vtk sphinx-fortran sphinxcontrib.images GitPython sphinx-rtd-theme pyvista f90nml scikit-learn pandas flatdict
ARG UID
ARG GID
ARG USER
ARG GROUP
RUN groupadd -g $GID -o $GROUP
RUN useradd -m -u $UID -g $GID -o -s /bin/bash $USER
USER $USER
WORKDIR /app
COPY execute-test /app/
ENV PENCIL_HOME=/app/pencil-code
ENV PYTHONPATH=/app/pencil-code/python:$PYTHONPATH
ENV PATH=/app:$PATH
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
RUN echo 'alias idl=gdl' >> ~/.bashrc
RUN echo 'source $PENCIL_HOME/sourceme.sh' >> ~/.bashrc
RUN echo 'export GDL_PATH=+/usr/share/gnudatalanguage/lib:$IDL_PATH' >> ~/.bashrc
WORKDIR /app/pencil-code
