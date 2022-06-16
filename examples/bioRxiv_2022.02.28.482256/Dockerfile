FROM buildpack-deps:focal

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    cmake=3.16.3-1ubuntu1 \
    libmpich-dev=3.3.2-2build1 \
    mpich=3.3.2-2build1 \
    doxygen=1.8.17-0ubuntu2 \
    libboost-dev=1.71.0.0ubuntu2 \
    libgsl-dev=2.5+dfsg-6build1 \
    fonts-humor-sans=1.0-3 \
    cython3=0.29.14-0.1ubuntu3 \
    python3-dev=3.8.2-0ubuntu2 \
    python3-pip=20.0.2-5ubuntu1.6 \
    python3-scipy=1.3.3-3build1 \
    python3-pandas=0.25.3+dfsg-7 \
    python3-seaborn=0.10.0-1 \
    python3-astropy=4.0-4ubuntu1 \
    python3-yaml=5.3.1-1ubuntu0.1 \
    python3-sympy=1.5.1-2.1 \
    ipython3=7.13.0-1 \
    jupyter=4.6.3-3 \
    jupyter-notebook=6.0.3-2 \
    antlr4=4.7.2-2 \
    libhdf5-dev=1.10.4+repack-11ubuntu1 \
    bison=2:3.5.1+dfsg-1 \
    flex=2.6.4-6.2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

    

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10 && \
    update-alternatives --install /usr/bin/ipython ipython /usr/bin/ipython3 10


# install mpi4py and h5py (the deb packages may depend on libopenmpi which we do not want)
RUN pip install mpi4py==3.1.3 && \
    pip install h5py==3.7.0  # 2.10.0

# install matplotlib to a newer versions
RUN pip install matplotlib==3.5.2


# Install NEST 3.3
RUN git clone --depth 1 -b v3.3 https://github.com/nest/nest-simulator /usr/src/nest-simulator && \
  mkdir nest-build && \
  cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/nest/ \
        -Dwith-boost=ON \
        -Dwith-ltdl=ON \
        -Dwith-gsl=ON \
        -Dwith-readline=ON \
        -Dwith-python=ON \
        -Dwith-mpi=ON \
        -Dwith-openmp=ON \
        -S /usr/src/nest-simulator \
        -B nest-build && \
  make -j4 -C nest-build && \
  make -C nest-build install

# clean up install/build files
RUN rm -r /usr/src/nest-simulator
RUN rm -r nest-build

# Add NEST binary folder to PATH
RUN echo "source /opt/nest/bin/nest_vars.sh" >> root/.bashrc


# ---- install NESTML -----
RUN pip install git+https://github.com/nest/nestml.git@a3a1b0d9fdb26e53aae45d8520c1f5b7dc97eaa3#egg=pynestml


# ---- install NEURON ----
RUN git clone --depth 1 -b 8.0.2 https://github.com/neuronsimulator/nrn.git /usr/src/nrn
RUN mkdir nrn-bld

RUN cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/nrn/ \
  -DCURSES_NEED_NCURSES=ON \
  -DNRN_ENABLE_INTERVIEWS=OFF \
  -DNRN_ENABLE_MPI=ON \
  -DNRN_ENABLE_RX3D=OFF \
  -DNRN_ENABLE_PYTHON=ON \
  -S /usr/src/nrn \
  -B nrn-bld

RUN cmake --build nrn-bld --parallel 4 --target install

# add nrnpython to PYTHONPATH
ENV PYTHONPATH /opt/nrn/lib/python:${PYTHONPATH}

# add NEURON binaries to PATH
ENV PATH /opt/nrn/bin:${PATH}

# clean up
RUN rm -r /usr/src/nrn
RUN rm -r nrn-bld

# ---- pip install some additional things
RUN pip install pymoo==0.4.2.2
RUN pip install git+https://github.com/NeuralEnsemble/parameters@b95bac2bd17f03ce600541e435e270a1e1c5a478#egg=parameters


# ---- install LFPykernels (main branch) -----
# RUN pip install git+https://github.com/LFPy/LFPykernels@29d3fcb40101c0bedca775a2985893de312488ad#egg=lfpykernels
RUN pip install LFPy==2.2.6
RUN pip install git+https://github.com/LFPy/LFPykernels@revision#egg=lfpykernels


# If running with Singularity, run the below line in the host.
# PYTHONPATH set here doesn't carry over:
# export SINGULARITYENV_PYTHONPATH=/opt/nest/lib/python3.8/site-packages
# Alternatively, run "source /opt/local/bin/nest_vars.sh" while running the container
