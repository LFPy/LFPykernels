# Codes for Hagen et al., (2022)

These codes correspond to results shown in the peer-reviewed manuscript:

Hagen E, Magnusson SH, Ness TV, Halnes G, Babu PN, et al. (2022)
Brain signal predictions from multi-scale networks using a linearized framework.
PLOS Computational Biology 18(8): e1010353. <https://doi.org/10.1371/journal.pcbi.1010353>

Bibtex format:

    @article{10.1371/journal.pcbi.1010353,
        doi = {10.1371/journal.pcbi.1010353},
        author = {Hagen, Espen AND Magnusson, Steinn H. AND Ness, Torbjørn V. AND Halnes, Geir AND Babu, Pooja N. AND Linssen, Charl AND Morrison, Abigail AND Einevoll, Gaute T.},
        journal = {PLOS Computational Biology},
        publisher = {Public Library of Science},
        title = {Brain signal predictions from multi-scale networks using a linearized framework},
        year = {2022},
        month = {08},
        volume = {18},
        url = {https://doi.org/10.1371/journal.pcbi.1010353},
        pages = {1-51},
        number = {8},
    }

If you use or refer to this work, please cite it as above.
Adaptations or modifications of this work should comply with the provided `LICENSE` file at the top level of this repository.  

## Requirements

For the preprint above, we mainly used the following software and versions thereof:

    GCC 11.2.0,
    mpich 3.4.2,
    Python 3.9.6,
    ipython 7.27.0,
    jupyter-notebook 6.0.3,
    numpy 1.21.3,
    scipy 1.7.1,
    matplotlib 3.4.3,
    pandas 1.3.4,
    seaborn 0.11.2,
    pymoo 0.4.2.2,
    mpi4py 3.1.3,
    h5py 3.5.0,
    NEURON 8.0.2,
    MEAutility 1.5.0,
    LFPykit 0.4,
    LFPy 2.2.6,
    LFPykernels 0.1.rc8 (https://github.com/LFPy/LFPykernels, git SHA: 4fd79ab),
    NEST 3.1 (https://github.com/nest/nest-simulator, git SHA: 512022e54),
    NESTML 4.0-post-dev (https://github.com/nest/nestml, git SHA: 0b251ec),
    parameters 0.2.1 (https://github.com/NeuralEnsemble/parameters, git SHA:b95bac2).

On any particular system or Python virtual environment, similar software versions may likely give similar results but exact replication of results is unlikely.
Note that simulations of networks and corresponding hybrid scheme simulations **do** require access to an HPC resource (see details in paper).
Postprocessing and plotting can be done on a normal laptop.

## Running

Make sure that above requirements are met (for a self-contained solution see contents of `/.Dockerfile`). Then some preparatory steps are required:

Fetch Hay et al. 2011 model files needed for some simulations using Python:

    ipython
    # import modules
    >>> from urllib.request import urlopen
    >>> import zipfile
    >>> import ssl
    # download
    >>> url = '{}{}'.format('http://senselab.med.yale.edu/ModelDB/eavBinDown.asp',
                        '?o=139653&a=23&mime=application/zip')
    >>> u = urlopen(url, context=ssl._create_unverified_context())
    >>> localFile = open('L5bPCmodelsEH.zip', 'wb')
    >>> localFile.write(u.read())
    >>> localFile.close()
    # unzip:
    >>> myzip = zipfile.ZipFile('L5bPCmodelsEH.zip', 'r')
    >>> myzip.extractall('.')
    >>> myzip.close()
    >>> quit()

Compile NMODL files in folder `./mod`:

    cd mod && nrnivmodl && cd -

Start jobs on the HPC resource via Slurm (<https://slurm.schedmd.com/documentation.html>).
These files must first be modified by the user for the particular HPC resource in question:

    python run_pscan.py
    python run_hay2011_pscan.py

In case all jobs finished without errors, fit some point-neuron network parameters to mimic the spiking activity of the ground-truth network simulation:

    sbatch Fit_LIF_net.job

Finally, if all above steps finished without errors, the figures etc. from the preprint was generated by the provided jupyter notebooks. A jupyter session can be invoked from the terminal issuing:

    jupyter-notebook
    # or
    jupyter-lab  # (recommended)

Note that total file output generated will be several GB, hence these are not provided in this repository. We are aiming to deposit these e.g., on Zenodo.org.

## Docker

Docker (<https://www.docker.com>) provides a solution for packaging all project requirements in a single container.
This can be used for simulations and analysis, and may be supported on the HPC resource (e.g., via Singularity).
Make sure the Docker client is running on the host. Then:

    # build image using Docker:
    docker build -t lfpykernels - < Dockerfile

    # start container mounting local file system, then open a jupyter-lab session:
    docker run --mount type=bind,source="$(pwd)",target=/opt/data -it -p 5000:5000 lfpykernels
    /# cd /opt/data/
    /# jupyter-lab --ip 0.0.0.0 --port=5000 --no-browser --allow-root
    # take note of the URL printed to the terminal, and open it in a browser on the host.

    # oneliner (open URL, then browse to `/opt/data/` and open notebooks):
    docker run --mount type=bind,source="$(pwd)",target=/opt/data -it -p 5000:5000 lfpykernels jupyter-lab --ip 0.0.0.0 --port=5000 --no-browser --allow-root
    # take note of the URL printed to the terminal, and open it in a browser on the host.

    # A working Python/MPI environment should be present in the running container. Hence scripts can be run interactively issuing:
    docker run --mount type=bind,source="$(pwd)",target=/opt/data -it -p 5000:5000 lfpykernels
    /# cd /opt/data/

    # start an interactive Python session
    /# ipython
    >>> import lfpykernels  # etc
    >>> quit()

    # run a simulation with MPI, assuming we have access to 1024 physical CPU cores (also make sure that parameter files have been created by an earlier call to `python run_pscan.py`)
    /# mpiexec -n 1024 python example_network.py c8a7f7ed6495c9fc60060f31d788208f

## Singularity

Instructions to create Singularity (<https://docs.sylabs.io/guides/3.5/user-guide/index.html>#) container named `lfpykernels.sif` using the provided `Dockerfile` recipe,
typically on a HPC resource providing corresponding build tools (cf. <https://apps.fz-juelich.de/jsc/hps/jusuf/cluster/container-runtime.html>).

    module --force purge
    module load Stages/2022  GCCcore/.11.2.0 Apptainer-Tools/2022

    sib upload ./Dockerfile lfpykernels
    sib build --recipe-name lfpykernels --blocking
    # wait ... 
    sib download --recipe-name lfpykernels
    sib remove --recipe-name lfpykernels  # remove recipe
    # test installation
    singularity exec lfpykernels.sif ipython3 -i

Compile NMODL files using the `nrnivmodl` script included in the container:

    cd mod && rm -rf x86_64 && singularity exec ../lfpykernels.sif nrnivmodl && cd ..  # compile NMODL files for the container

Ensure that Slurm jobscripts are configured for singularity, calling the built-in python executable:
    unset DISPLAY
    module --force purge
    module load Stages/2022  GCCcore/.11.2.0 Apptainer-Tools/2022 GCC/11.2.0 ParaStationMPI/5.5.0-1
    srun --mpi=pmi2 singularity exec lfpykernels.sif python3 -u example_network.py adb947bfb931a5a8d09ad078a6d256b0

## Files

Packaged files:

- `run_pscan.py`: main script that submit SLURM jobs for all ball-n-sticks network simulations and hybrid scheme approximations. Adapt as necessary for the HPC resource in question.
- `run_hay2011_pscan.py`: Similar to `run_pscan.py` but for simulation runs using the Hay et al. 2011 L5bPC model.
- `example_network_methods.py`: module with some shared methods required by simulation scripts and analysis notebooks
- `*_network_parameters.py`: shared parameters for the ball-n-sticks or hay2011 based networks
- `*_network.py`: main ball-and-sticks network simulation script
- `*_network_reconstruction.py`: hybrid-scheme implementation of extracellular signal predictions from spike output of `example_network.py`
- `*_network_kernel.py`: hybrid-scheme calculations of kernels
- `plotting.py`: some shared plotting routines
- `FIR_filter.nestml`: NESTML implementation of FIR filter node
- `mod/*.mod`: NMODL language files with quasi-linear and passive-frozen versions of active ion channels of the Hay et al. (2011) L5bPC model
- `BallAndSticks*.hoc`: cell geometry/template files in NEURON's HOC language
- `Fit_LIF_net.py`: script for fitting parameters of spiking point-neuron network
- `Fit_LIF_net.job`: Slurm job script for `Fit_LIF_net.py`. Adapt as necessary for HPC resource.
- `LIF_net.py`: Common classes/methods used by `Fit_LIF_net.py` and `Figure16.ipynb`
- `README.md`: This file

Output:

- `L5bPCmodelsEH/`: Layer 5b pyramidal cell model by Hay et al., (2011)
- `jobs/*.job`: Slurm job scripts generated by `run*pscan.py` script
- `output/*/*`: File output by different simulation scripts
- `logs/*.txt`: Job output to STDOUT/STDERR
- `parameters/*.txt`: parameter dictionary files (for `parameters.ParameterSet`) generated by `run*pscan.py` scripts
- `figures/`: figure output directory for below notebooks
- `PS*.txt`: `parameters.ParameterSpace` files for the different simulations started by `run_pscan.py`
- `hay2011_PS*.txt`: `parameters.ParameterSpace` files for the different simulations started by `run_hay2011_pscan.py`
- `Fit_LIF_net.h5`: file output created by successful runs of `Fit_LIF_net.py` containing parameters and metrics for each generation running the NSGA-II multi-objective optimization algorithm

Notebooks:

- `Figure*.ipynb`: Jupyter notebooks used to generate the corresponding figures post simulation plus some analysis

A 4.6GB .zip file with some precomputed output data corresponding to Hagen E et al. (2022) in PLOS Comput Biol is available for download from <https://doi.org/10.5281/zenodo.6992596>
Unzip it's content in this directory, and one should be able to execute all notebooks.

### Issues

For problems/questions with regards to these codes, please raise an issue at <https://github.com/LFPy/LFPykernels/issues>