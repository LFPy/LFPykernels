FROM buildpack-deps:focal

RUN apt-get update && \
    apt-get install -y \
    cmake \
    libmpich-dev \
    mpich \
    doxygen \
    libboost-dev \
    libgsl-dev \
    cython3 \
    python3-dev \
    python3-pip \
    python3-numpy \
    python3-scipy \
    python3-matplotlib \
    python3-pandas \
    ipython3 \
    jupyter \
    libhdf5-dev \
    bison flex \
    python3-mpi4py python3-h5py


RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10 && \
    update-alternatives --install /usr/bin/ipython ipython /usr/bin/ipython3 10


# Install NEST 3.1
RUN git clone --depth 1 -b v3.1 https://github.com/nest/nest-simulator && \
  mkdir nest-build && \
  cd  nest-build && \
  cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/nest/ \
        -Dwith-boost=ON \
        -Dwith-ltdl=ON \
        -Dwith-gsl=ON \
        -Dwith-readline=ON \
        -Dwith-python=ON \
        -Dwith-mpi=ON \
        -Dwith-openmp=ON \
        ../nest-simulator && \
  make -j4 && \
  make install && \
  cd ..

# clean up install/build files
RUN rm -r nest-simulator
RUN rm -r nest-build

# Add NEST binary folder to PATH
RUN echo "source /opt/nest/bin/nest_vars.sh" >> root/.bashrc


# ---- install NESTML -----
RUN pip install git+https://github.com/nest/nestml.git@a3a1b0d9fdb26e53aae45d8520c1f5b7dc97eaa3#egg=pynestml


# ---- install NEURON ----
RUN git clone --depth 1 -b 8.0.0 https://github.com/neuronsimulator/nrn.git
RUN mkdir nrn-bld && cd nrn-bld

RUN cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/nrn/ \
  -DCURSES_NEED_NCURSES=ON \
  -DNRN_ENABLE_INTERVIEWS=OFF \
  -DNRN_ENABLE_MPI=ON \
  -DNRN_ENABLE_RX3D=OFF \
  -DNRN_ENABLE_PYTHON=ON \
  ../nrn

RUN cmake --build . --parallel 4 --target install && \
  cd ..

# add nrnpython to PYTHONPATH
ENV PYTHONPATH /opt/nrn/lib/python:${PYTHONPATH}

# clean up
RUN rm -r nrn
RUN rm -r nrn-bld


# ---- install LFPykernels -----
RUN pip install git+https://github.com/LFPy/LFPykernels@main#egg=lfpykernels

# ---- pip install some additional things
RUN pip install pymoo==0.4.2.2
RUN pip install git+https://github.com/NeuralEnsemble/parameters@master#egg=parameters



# If running with Singularity, run the below line in the host.
# PYTHONPATH set here doesn't carry over:
# export SINGULARITYENV_PYTHONPATH=/opt/nest/lib/python3.8/site-packages
# Alternatively, run "source /opt/local/bin/nest_vars.sh" while running the container
