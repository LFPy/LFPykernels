# LFPykernels

The ``LFPykernels`` package incorporates forward-model based calculations of causal spike-signal
impulse response functions for finite-sized neuronal network models.


## Build Status

[![DOI](https://zenodo.org/badge/424143558.svg)](https://zenodo.org/badge/latestdoi/424143558)
[![Coverage Status](https://coveralls.io/repos/github/LFPy/LFPykernels/badge.svg?branch=main)](https://coveralls.io/github/LFPy/LFPykernels?branch=main)
[![Documentation Status](https://readthedocs.org/projects/lfpykernels/badge/?version=latest)](https://lfpykernels.readthedocs.io/en/latest/?badge=latest)
[![Lintly flake8 checks](https://github.com/LFPy/lfpykernels/workflows/Lintly%20flake8%20checks/badge.svg)](https://github.com/LFPy/LFPykernels/actions?query=workflow%3A%22Lintly+flake8+checks%22)
[![Python application](https://github.com/LFPy/LFPykernels/workflows/Python%20application/badge.svg)](https://github.com/LFPy/LFPykernels/actions?query=workflow%3A%22Python+application%22)
[![Upload Python Package](https://github.com/LFPy/LFPykernels/workflows/Upload%20Python%20Package/badge.svg)](https://pypi.org/project/LFPykernels)
[![Conda Recipe](https://img.shields.io/badge/recipe-lfpykernels-green.svg)](https://anaconda.org/conda-forge/lfpykernels)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/lfpykernels.svg)](https://anaconda.org/conda-forge/lfpykernels)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/lfpykernels.svg)](https://anaconda.org/conda-forge/lfpykernels)
[![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/lfpykernels.svg)](https://anaconda.org/conda-forge/lfpykernels)
[![License](http://img.shields.io/:license-GPLv3+-green.svg)](http://www.gnu.org/licenses/gpl-3.0.html)


## Citation

If you use this software, please cite it as (change version accordingly):

    Hagen, Espen. (2021). LFPykernels (version/git-SHA/git-tag). Zenodo. https://doi.org/10.5281/zenodo.5720619


BibTex format:

    @software{hagen_espen_2021_5720619,
      author       = {Hagen, Espen},
      title        = {LFPykernels},
      month        = nov,
      year         = 2021,
      note         = {If you use this software, please cite it as below.},
      publisher    = {Zenodo},
      version      = {version/git-SHA/git-tag},
      doi          = {10.5281/zenodo.5720619},
      url          = {https://doi.org/10.5281/zenodo.5720619}
    }


These codes correspond to results shown in the preprint manuscript:

    Brain signal predictions from multi-scale networks using a linearized framework
    Espen Hagen, Steinn H. Magnusson, Torbjørn V. Ness, Geir Halnes, Pooja N. Babu, Charl Linssen, Abigail Morrison, Gaute T. Einevoll
    bioRxiv 2022.02.28.482256; doi: https://doi.org/10.1101/2022.02.28.482256


Bibtex format:

    @article {Hagen2022.02.28.482256,
    	author = {Hagen, Espen and Magnusson, Steinn H. and Ness, Torbjørn V. and Halnes, Geir and Babu, Pooja N. and Linssen, Charl and Morrison, Abigail and Einevoll, Gaute T.},
    	title = {Brain signal predictions from multi-scale networks using a linearized framework},
    	elocation-id = {2022.02.28.482256},
    	year = {2022},
    	doi = {10.1101/2022.02.28.482256},
    	publisher = {Cold Spring Harbor Laboratory},
    	URL = {https://www.biorxiv.org/content/early/2022/03/02/2022.02.28.482256},
    	eprint = {https://www.biorxiv.org/content/early/2022/03/02/2022.02.28.482256.full.pdf},
    	journal = {bioRxiv}
    }

If you use or refer to this work, please cite it as above.
Adaptations or modifications of this work should comply with the provided `LICENSE` file provided with this repository.  


## Features

The ``LFPykernels`` package incorporates forward-model based calculations of causal spike-signal
impulse response functions for finite-sized neuronal network models.
The signals considered are low-frequency extracellular potentials ("local field potential" - LFP)
or current dipole moments (and by extension EEG and MEG like signals) that are
thought to mainly stem from synaptic currents and associated return currents.
The basic idea is that the effect of any spike event in each presynaptic
population on each signal type can be captured by single linearised multicompartment neuron
models representative of each population and simultaneously accounting for known distributions of
cells and synapses in space, distributions of delays, synaptic currents and associated return currents.

A scientific publication describing the present methodology in detail is planned.
The intended use for filter kernels predicted using ``LFPykernels`` is forward-model based signal predictions
from neuronal network simulation frameworks using simplified neuron representations like leaky integrate-and-fire
point neurons or rate-based neurons.
Let ``nu_X(t)`` describe presynaptic population spike rates in units of spikes/dt
and ``H_YX(r, tau)`` predicted spike-signal kernels for the connections between presynaptic populations ``X`` and
postsynaptic populations ``Y`` the full signal may then be computed via the sum over linear convolutions
```
V(r, t) = sum_X sum_Y conv(nu_X, H_YX)(r, t)
```

A more elaborate example combining kernel predictions with a spiking point-neuron network simulation is provided in the example notebook
https://github.com/LFPy/LFPykernels/blob/main/examples/LIF_net_forward_model_predictions.ipynb

For questions, please raise an issue at https://github.com/LFPy/LFPykernels/issues.


## Usage

Example prediction of kernel function ``H`` mapping spike events of a
presynaptic inhibitory population ``X=='I'`` to extracellular potential contributions by a
postsynaptic excitatory population ``Y=='E'`` (see https://github.com/LFPy/LFPykernels/blob/main/examples/README_example.ipynb):

    import matplotlib.pyplot as plt
    import scipy.stats as st
    import numpy as np
    from lfpykernels import GaussCylinderPotential, KernelApprox
    import neuron

    # recompile mod files if needed
    mech_loaded = neuron.load_mechanisms('mod')
    if not mech_loaded:
        os.system('cd mod && nrnivmodl && cd -')
        mech_loaded = neuron.load_mechanisms('mod')
    print(f'mechanisms loaded: {mech_loaded}')

    # misc parameters
    dt = 2**-4  # time resolution (ms)
    t_X = 500  # time of synaptic activations (ms)
    tau = 50  # duration of impulse response function after onset (ms)
    Vrest = -65  # assumed average postsynaptic potential (mV)

    X=['E', 'I']   # presynaptic population names
    N_X = np.array([8192, 1024])  # presynpatic population sizes
    Y = 'E' # postsynaptic population
    N_Y = 8192  # postsynaptic population size
    C_YX = np.array([0.05, 0.05])  # pairwise connection probability between populations X and Y
    nu_X = {'E': 2.5, 'I': 5.0}  # assumed spike rates of each population (spikes/s)
    g_eff = True  # account for changes in passive leak due to persistent synaptic activations

    def set_passive(cell, Vrest):
        """Insert passive leak channel across all sections

        Parameters
        ----------
        cell: object
            LFPy.NetworkCell like object
        Vrest: float
            Steady state potential
        """
        for sec in cell.template.all:
            sec.insert('pas')
            sec.g_pas = 0.0003  # (S/cm2)
            sec.e_pas = Vrest  # (mV)

    # parameters for LFPy.NetworkCell representative of postsynaptic population
    cellParameters={
        'templatefile': 'BallAndSticksTemplate.hoc',
        'templatename': 'BallAndSticksTemplate',
        'custom_fun': [set_passive],
        'custom_fun_args': [{'Vrest': Vrest}],
        'templateargs': None,
        'delete_sections': False,
        'morphology': 'BallAndSticks_E.hoc'}

    populationParameters={
            'radius': 150.0,  # population radius (µm)
            'loc': 0.0,  # average depth of cell bodies (µm)
            'scale': 75.0}  # standard deviation (µm)

    # Predictor for extracellular potentials across depth assuming planar disk source
    # elements convolved with Gaussian along z-axis.
    # See https://lfpykernels.readthedocs.io/en/latest/#class-gausscylinderpotential for details
    probe = GaussCylinderPotential(
        cell=None,
        z=np.linspace(1000., -200., 13),  # depth of contacts (µm)
        sigma=0.3,  # tissue conductivity (S/m)
        R=populationParameters['radius'],  #
        sigma_z=populationParameters['scale'],
        )

    # Create KernelApprox object. See https://lfpykernels.readthedocs.io/en/latest/#class-kernelapprox for details
    kernel = KernelApprox(
        X=X,
        Y=Y,
        N_X=N_X,
        N_Y=N_Y,
        C_YX=C_YX,
        cellParameters=cellParameters,
        populationParameters=populationParameters,
        # function and parameters used to estimate average multapse count:
        multapseFunction=st.truncnorm,
        multapseParameters=[
            {'a': (1 - 2.) / .6, 'b': (10 - 2.) / .6, 'loc': 2.0, 'scale': 0.6},
            {'a': (1 - 5.) / 1.1, 'b': (10 - 5.) / 1.1, 'loc': 5.0, 'scale': 1.1}],
        # function and parameters for delay distribution from connections between a
        # population in X onto population Y:
        delayFunction=st.truncnorm,
        delayParameters=[{'a': -2.2, 'b': np.inf, 'loc': 1.3, 'scale': 0.5},
                         {'a': -1.5, 'b': np.inf, 'loc': 1.2, 'scale': 0.6}],
        # parameters for synapses from connections by populations X onto Y
        synapseParameters=[
            {'weight': 0.00012, 'syntype': 'Exp2Syn', 'tau1': 0.2, 'tau2': 1.8, 'e': 0.0},
            {'weight': 0.002, 'syntype': 'Exp2Syn', 'tau1': 0.1, 'tau2': 9.0, 'e': -80.0}],
        # parameters for spatial synaptic connectivity by populations X onto Y
        synapsePositionArguments=[
            {'section': ['apic', 'dend'],
             'fun': [st.norm],
             'funargs': [{'loc': 50.0, 'scale': 100.0}],
             'funweights': [1.0]},
            {'section': ['soma', 'apic', 'dend'],
             'fun': [st.norm],
             'funargs': [{'loc': -100.0, 'scale': 100.0}],
             'funweights': [1.0]}],
        # parameters for extrinsic synaptic input
        extSynapseParameters={'syntype': 'Exp2Syn', 'weight': 0.0002, 'tau1': 0.2, 'tau2': 1.8, 'e': 0.0},
        nu_ext=40.,  # external activation rate (spikes/s)
        n_ext=450,  # number of extrinsic synapses  
        nu_X=nu_X,
    )

    # make kernel predictions for connection from populations X='I' onto Y='E'
    H = kernel.get_kernel(
        probes=[probe],
        Vrest=Vrest, dt=dt, X='I', t_X=t_X, tau=tau,
        g_eff=g_eff)


## Physical units

Notes on physical units used in `LFPykernels`:

- There are no explicit checks for physical units

- Transmembrane currents are assumed to be in units of (nA)

- Spatial information is assumed to be in units of (µm)

- Voltages are assumed to be in units of (mV)

- Extracellular conductivities are assumed to be in units of (S/m)

- current dipole moments are assumed to be in units of (nA µm)

- Magnetic fields are assumed to be in units of (nA/µm)

- Simulation times are assumed to be in units of (ms) with step size ∆t

- Spike rates are assumed to be in units of (# spikes / ∆t)


## Documentation

The online Documentation of `LFPykernels` can be found here:
https://lfpykernels.readthedocs.io/en/latest


## Dependencies

`LFPykernels` is implemented in Python and is written (and continuously tested) for `Python >= 3.7` (older versions may or may not work).
The main `LFPykernels` module depends on ``LFPy`` (https://github.com/LFPy/LFPy, https://LFPy.readthedocs.io).

Running all unit tests and example files may in addition require `py.test`, `matplotlib`,
`LFPy`.


## Installation

### From development sources (https://github.com/LFPy/LFPykernels)

Install the current development version on https://GitHub.com using `git` (https://git-scm.com):

    $ git clone https://github.com/LFPy/LFPykernels.git
    $ cd LFPykernels
    $ python setup.py install  # --user optional

or using `pip`:

    $ pip install .  # --user optional

For active development, link the repository location

    $ pip install -e .  # --user optional

### Installation of stable releases on PyPI.org (https://www.pypi.org)

Installing stable releases from the Python Package Index (https://www.pypi.org/project/lfpykernels):

    $ pip install lfpykernels  # --user optional

To upgrade the installation using pip:

    $ pip install --upgrade --no-deps lfpykernels


### Installation of stable releases on conda-forge (https://conda-forge.org)

Installing `lfpykernels` from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

    $ conda config --add channels conda-forge

Once the `conda-forge` channel has been enabled, `lfpykernels` can be installed with:

    $ conda install lfpykernels

It is possible to list all of the versions of `lfpykernels` available on your platform with:

    $ conda search lfpykernels --channel conda-forge


## Docker

We provide a Docker (https://www.docker.com) container recipe file with LFPykernels etc.
To get started, install Docker and issue either:

    # build Dockerfile from GitHub
    $ docker build -t lfpykernels https://raw.githubusercontent.com/LFPy/LFPykernels/main/Dockerfile
    $ docker run -it -p 5000:5000 lfpykernels


or

    # build local Dockerfile (obtained by cloning repo, checkout branch etc.)
    $ docker build -t lfpykernels - < Dockerfile
    $ docker run -it -p 5000:5000 lfpykernels


If the docker file should fail for some reason it is possible to store the build log and avoid build caches by issuing

    $ docker build --no-cache --progress=plain -t lfpykernels - < Dockerfile 2>&1 | tee lfpykernels.log


For successful builds, the ``--mount`` option can be used to mount a folder on the host to a target folder as:

    $ docker run --mount type=bind,source="$(pwd)",target=/opt/data -it -p 5000:5000 lfpykernels


which mounts the present working dirctory (``$(pwd)``) to the ``/opt/data`` directory of the container.
Try mounting the ``LFPykernels`` source directory for example (by setting ``source="<path-to-LFPykernels>"``).
Various example files can then be found in the folder ``/opt/data/examples/``
when the container is running.

Jupyter notebook servers running from within the
container can be accessed after invoking them by issuing:

    $ cd /opt/data/examples/
    $ jupyter-notebook --ip 0.0.0.0 --port=5000 --no-browser --allow-root


and opening the resulting URL in a browser on the host computer, similar to:
http://127.0.0.1:5000/?token=dcf8f859f859740fc858c568bdd5b015e0cf15bfc2c5b0c1
