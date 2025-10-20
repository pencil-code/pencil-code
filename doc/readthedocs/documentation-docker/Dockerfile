FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y gcc gfortran pip openmpi-bin libopenmpi-dev git libxrender1 x11-apps python3-tk
RUN pip3 install numpy scipy matplotlib setuptools==68.2.2 eqtools==1.4.0 h5py plotly retrying astropy Cython dill pexpect tqdm vtk sphinxcontrib.images GitPython sphinx-rtd-theme pyvista f90nml scikit-learn pandas flatdict sphinx rstcloth
WORKDIR /app/pencil-code/doc/readthedocs
RUN  git config --global --add safe.directory /app/pencil-code
