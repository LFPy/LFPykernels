#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''LFPykernels setuptools file

'''

import os
import sys
import shutil
from distutils.spawn import spawn
import setuptools

d = {}
exec(open(os.path.join('lfpykernels', 'version.py')).read(), None, d)
version = d['version']

# try and locate the nrnivmodl or mknrndll script of NEURON in PATH so that the
# NEURON NMODL files LFPy/test/*.mod can be compiled in place and be copied
# as part of the package_data, allowing unit tests to run
if not any(arg in sys.argv for arg in ['sdist', 'upload']):
    if shutil.which('nrnivmodl') is not None:
        os.chdir(os.path.join('lfpykernels', 'tests'))
        for path in ['x86_64', 'arm64', 'aarch64']:
            if os.path.isdir(path):
                shutil.rmtree(path)
        spawn([shutil.which('nrnivmodl')])
        os.chdir(os.path.join('..', '..'))
    elif shutil.which('mknrndll') is not None:
        os.chdir(os.path.join('lfpykernels', 'tests'))
        if os.path.isfile("nrnmech.dll"):
            os.remove("nrnmech.dll")
        spawn([shutil.which('mknrndll')])
        os.chdir(os.path.join('..', '..'))
    else:
        print("nrnivmodl/mknrndll script not found in PATH, thus NMODL " +
              "files could not be compiled. lfpykernels test functions" +
              "will fail")

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='LFPykernels',
    version=version,
    author='LFPy-team',
    author_email='lfpy@users.noreply.github.com',
    description=('Causal spike-signal impulse response functions for ' +
                 'finite-sized neuronal network models'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/LFPy/LFPykernels',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Utilities',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Development Status :: 4 - Beta',
    ],
    python_requires='>=3.7',
    install_requires=['LFPy>=2.2.2'],
    package_data={
        'lfpykernels': [
            os.path.join('tests', '*.py'),
            os.path.join('tests', '*.hoc'),
            os.path.join('tests', '*.mod'),
            ]},
    include_package_data=True,
    extras_require={
        'tests': [
            'pytest',
            'sympy'],
        'docs': [
            'sphinx',
            'numpydoc',
            'sphinx_rtd_theme',
            'recommonmark'],
    },
    dependency_links=[],
    provides=['lfpykernels'],
    zip_safe=False)
