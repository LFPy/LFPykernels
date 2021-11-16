# LFPykernels

The `LFPykernels` package incorporates calculations of causal spike-signal
impulse response functions for finite-sized neuronal network models.


## Build Status

[![DOI](https://zenodo.org/badge/XXXXXXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXXXXXXX)
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


## Features

The `LFPykernels` package incorporates calculations of causal spike-signal
impulse response functions for finite-sized neuronal network models.


## Usage

A basic usage example using a mock 3-segment stick-like neuron,
treating each segment as a point source in a linear, isotropic and homogeneous volume conductor,
computing the extracellular potential in ten different locations
alongside the cell geometry:

    >>> # imports
    >>> import numpy as np
    >>> from lfpykernels import *


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

`LFPykernels` is implemented in Python and is written (and continuously tested) for `Python >= 3.7`.
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

Installing from the Python Package Index (https://www.pypi.org/project/lfpykernels):

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


The ``--mount`` option can be used to mount a folder on the host to a target folder as:

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
