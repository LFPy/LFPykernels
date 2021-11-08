#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''LFPykernels setuptools file

'''

import os
import setuptools

d = {}
exec(open(os.path.join('lfpykernels', 'version.py')).read(), None, d)
version = d['version']


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
        'Programming Language :: Python :: 3.6',
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
    python_requires='>=3.6',
    install_requires=['LFPy>=2.2.2'],
    package_data={
        'lfpykernels': [
            os.path.join(
                    'tests',
                    '*.py')]},
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
