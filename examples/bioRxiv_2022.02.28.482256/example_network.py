#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Demonstrate usage of LFPy.Network with network of ball-and-stick type
morphologies with active HH channels inserted in the somas and passive-leak
channels distributed throughout the apical dendrite. The corresponding
morphology and template specifications are in the files BallAndStick.hoc and
BallAndStickTemplate.hoc.

Execution (w. MPI):

    mpirun -np 2 python example_network.py

Copyright (C) 2017 Computational Neuroscience Group, NMBU.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

"""
# import modules:
import os
import sys
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import scipy.signal as ss
import h5py
from parameters import ParameterSet, ParameterSpace
from mpi4py import MPI
from time import time
import neuron
from LFPy import Network, Synapse, RecExtElectrode, CurrentDipoleMoment, \
    LaminarCurrentSourceDensity
from plotting import draw_lineplot, remove_axis_junk
import example_network_parameters as params
import hashlib
import json


# set up MPI variables:
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

# avoid same sequence of random numbers from numpy and neuron on each RANK,
# e.g., in order to draw unique cell and synapse locations and random synapse
# activation times
GLOBALSEED = 1234
np.random.seed(GLOBALSEED + RANK)

##########################################################################
# Function declarations
##########################################################################


if __name__ == '__main__':
    ##########################################################################
    # Main simulation
    ##########################################################################

    neuron.load_mechanisms('mod')

    #######################################
    # tic tac
    #######################################
    tic = time()

    #######################################
    # Capture command line values
    #######################################
    # load parameter file
    md5 = sys.argv[1]
    pset = ParameterSet(os.path.join('parameters', '{}.txt'.format(md5)))

    # fetch parameter values
    weight_EE = pset['weight_EE']
    weight_IE = pset['weight_IE']
    weight_EI = pset['weight_EI']
    weight_II = pset['weight_II']
    weight_scaling = pset['weight_scaling']
    n_ext = pset['n_ext']

    #####################################
    # save connections only for those parameter combinations that exist in
    # PS1 (and PS2) to save disk space
    ######################################
    def get_save_connections_value(pset, md5):
        PS1 = ParameterSpace('PS1.txt')
        keys = list(PS1.keys())
        for key in keys:
            if key not in list(pset.keys()):
                PS1.pop(key)
        md5s = []
        for ps1pset in PS1.iter_inner():
            # sorted json dictionary
            js = json.dumps(ps1pset, sort_keys=True).encode()
            md5s += [hashlib.md5(js).hexdigest()]

        if md5 in md5s:
            return True
        else:
            return False

    save_connections = get_save_connections_value(pset, md5)

    ##########################################################################
    # Set up shared and population-specific parameters
    ##########################################################################
    # relative path for simulation output:
    OUTPUTPATH = os.path.join('output', md5)

    if RANK == 0:
        # create directory for output:
        if not os.path.isdir(OUTPUTPATH):
            os.mkdir(OUTPUTPATH)
        # remove old simulation output if directory exist
        else:
            for fname in os.listdir(OUTPUTPATH):
                os.unlink(os.path.join(OUTPUTPATH, fname))
    COMM.Barrier()

    # tic tac
    tac = time()
    if RANK == 0:
        with open(os.path.join(OUTPUTPATH, 'tic_tac.txt'), 'w') as f:
            f.write(f'step time_s\nsetup {tac - tic}\n')

    ##########################################################################
    # Parameters dependent on command line input
    ##########################################################################
    # scale connection weights
    weight_EE *= weight_scaling
    weight_IE *= weight_scaling
    weight_EI *= weight_scaling
    weight_II *= weight_scaling

    # synapse max. conductance (function, mean, st.dev., min.):
    weightArguments = [[dict(loc=weight_EE, scale=weight_EE / 10),
                        dict(loc=weight_IE, scale=weight_IE / 10)],
                       [dict(loc=weight_EI, scale=weight_EI / 10),
                        dict(loc=weight_II, scale=weight_II / 10)]]

    # instantiate Network:
    network = Network(OUTPUTPATH=OUTPUTPATH, **params.networkParameters)

    # create E and I populations:
    for j, (name, size, morphology) in enumerate(zip(params.population_names,
                                                     params.population_sizes,
                                                     params.morphologies)):
        popParams = deepcopy(params.populationParameters)
        popParams['cell_args'].update(dict(
            morphology=morphology
        ))
        network.create_population(name=name, POP_SIZE=size,
                                  **popParams)

        # create excitatory background synaptic activity for each cell
        # with Poisson statistics
        for cell in network.populations[name].cells:
            idx = cell.get_rand_idx_area_norm(section='allsec', nidx=n_ext[j])
            for i in idx:
                syn = Synapse(cell=cell, idx=i,
                              **params.extSynapseParameters)
                syn.set_spike_times_w_netstim(
                    interval=params.netstim_interval,
                    seed=np.random.rand() * 2**32 - 1)

    # tic tac
    tic = time()
    if RANK == 0:
        with open(os.path.join(OUTPUTPATH, 'tic_tac.txt'), 'a') as f:
            f.write(f'create {tic - tac}\n')

    # create connectivity matrices and connect populations:
    for i, pre in enumerate(params.population_names):
        for j, post in enumerate(params.population_names):
            # boolean connectivity matrix between pre- and post-synaptic
            # neurons in each population (postsynaptic on this RANK)
            connectivity = network.get_connectivity_rand(
                pre=pre, post=post,
                connprob=params.connectionProbability[i][j]
            )

            # connect network:
            (conncount, syncount) = network.connect(
                pre=pre, post=post,
                connectivity=connectivity,
                syntype=params.synapseModel,
                synparams=params.synapseParameters[i][j],
                weightfun=params.weightFunction,
                weightargs=weightArguments[i][j],
                minweight=params.minweight,
                delayfun=params.delayFunction,
                delayargs=params.delayArguments[i][j],
                mindelay=params.mindelay,
                multapsefun=params.multapseFunction,
                multapseargs=params.multapseArguments[i][j],
                syn_pos_args=params.synapsePositionArguments[i][j],
                save_connections=save_connections,
            )

    # tic tac
    tac = time()
    if RANK == 0:
        with open(os.path.join(OUTPUTPATH, 'tic_tac.txt'), 'a') as f:
            f.write(f'connect {tac - tic}\n')

    # set up extracellular recording device.
    # Here `cell` is set to None as handles to cell geometry is handled
    # internally
    electrode = RecExtElectrode(cell=None, **params.electrodeParameters)

    # set up recording of current dipole moments. Ditto with regards to
    # `cell` being set to None
    current_dipole_moment = CurrentDipoleMoment(cell=None)

    # set up recording of current source density. Ditto with regards to
    # `cell` being set to None
    csd = LaminarCurrentSourceDensity(cell=None, **params.csdParameters)

    # run simulation:
    SPIKES = network.simulate(
        probes=[electrode, current_dipole_moment, csd],
        **params.networkSimulationArguments
    )

    # tic tac
    tic = time()
    if RANK == 0:
        with open(os.path.join(OUTPUTPATH, 'tic_tac.txt'), 'a') as f:
            f.write(f'simulate {tic - tac}\n')

    # save spikes
    if RANK == 0:
        with h5py.File(os.path.join(OUTPUTPATH, 'spikes.h5'), 'w') as f:
            dtype = h5py.special_dtype(vlen=np.dtype('float'))
            for i, name in enumerate(['E', 'I']):
                subgrp = f.create_group(name)
                if len(SPIKES['gids'][i]) > 0:
                    subgrp['gids'] = np.array(SPIKES['gids'][i]).flatten()
                    dset = subgrp.create_dataset('times',
                                                 (len(SPIKES['gids'][i]),),
                                                 dtype=dtype)
                    for j, spt in enumerate(SPIKES['times'][i]):
                        dset[j] = spt
                else:
                    subgrp['gids'] = []
                    subgrp['times'] = []

    # collect somatic potentials across all RANKs to RANK 0:
    if RANK == 0:
        somavs = []
    for i, name in enumerate(params.population_names):
        for j, cell in enumerate(network.populations[name].cells):
            if j == 0:
                somavs_pop = cell.somav
            # else:
            #     somavs_pop = np.vstack((somavs_pop, cell.somav))
        if RANK == 0:
            for j in range(1, SIZE):
                somavs_pop = np.vstack((somavs_pop,
                                        COMM.recv(source=j, tag=15)))
            somavs.append(somavs_pop)
        else:
            COMM.send(somavs_pop, dest=0, tag=15)
        del somavs_pop

    if RANK == 0:
        with h5py.File(os.path.join(OUTPUTPATH, 'somav.h5'), 'w') as f:
            for i, name in enumerate(params.population_names):
                # f[name] = somavs[i]
                f[name] = ss.decimate(somavs[i],
                                      q=int(round(1 // network.dt)),
                                      axis=-1,
                                      zero_phase=True)

    # store lfpykit probe data to file
    if RANK == 0:
        for probe in [electrode, current_dipole_moment, csd]:
            with h5py.File(
                os.path.join(OUTPUTPATH,
                             '{}.h5'.format(probe.__class__.__name__)), 'w'
            ) as f:
                f['data'] = probe.data

    # tic tac
    tac = time()
    if RANK == 0:
        with open(os.path.join(OUTPUTPATH, 'tic_tac.txt'), 'a') as f:
            f.write(f'save {tac - tic}\n')

    ##########################################################################
    # Plot some output on RANK 0
    ##########################################################################

    if RANK == 0:
        # spike raster
        fig, ax = plt.subplots(1, 1)
        for name, spts, gids in zip(
                params.population_names, SPIKES['times'], SPIKES['gids']):
            t = []
            g = []
            for spt, gid in zip(spts, gids):
                t = np.r_[t, spt]
                g = np.r_[g, np.zeros(spt.size) + gid]
            inds = (t >= 500) & (t <= 1000)  # show [500, 1000] ms interval
            ax.plot(t[inds], g[inds], '.', ms=3, label=name)
        ax.legend(loc=1)
        remove_axis_junk(ax, lines=['right', 'top'])
        ax.set_xlabel('t (ms)')
        ax.set_ylabel('gid')
        ax.set_title('spike raster')
        fig.savefig(os.path.join(OUTPUTPATH, 'spike_raster.pdf'),
                    bbox_inches='tight')
        plt.close(fig)

        # somatic potentials
        fig = plt.figure()
        gs = GridSpec(4, 1)
        ax = fig.add_subplot(gs[:2])
        if len(somavs[0] > 8):
            data = ss.decimate(somavs[0][:8], q=16, axis=-1, zero_phase=True)
        else:
            data = ss.decimate(somavs[0], q=16, axis=-1, zero_phase=True),
        draw_lineplot(ax,
                      data,
                      dt=network.dt * 16,
                      T=(500, 1000),
                      scaling_factor=1.,
                      vlimround=16,
                      label='E',
                      scalebar=True,
                      unit='mV',
                      ylabels=False,
                      color='C0',
                      ztransform=True
                      )
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_ylabel('E')
        ax.set_title('somatic potentials')
        ax.set_xlabel('')

        ax = fig.add_subplot(gs[2:])
        if len(somavs[1] > 8):
            data = ss.decimate(somavs[1][:8], q=16, axis=-1, zero_phase=True)
        else:
            data = ss.decimate(somavs[1], q=16, axis=-1, zero_phase=True),
        draw_lineplot(ax,
                      data,
                      dt=network.dt * 16,
                      T=(500, 1000),
                      scaling_factor=1.,
                      vlimround=16,
                      label='I',
                      scalebar=True,
                      unit='mV',
                      ylabels=False,
                      color='C1',
                      ztransform=True
                      )
        ax.set_yticks([])
        ax.set_ylabel('I')

        fig.savefig(os.path.join(OUTPUTPATH, 'soma_potentials.pdf'),
                    bbox_inches='tight')
        plt.close(fig)

        # extracellular potentials, E and I contributions, sum
        fig, axes = plt.subplots(1, 3, figsize=(6.4, 4.8))
        fig.suptitle('extracellular potentials')
        for i, (ax, name, label) in enumerate(zip(axes, ['E', 'I', 'imem'],
                                                  ['E', 'I', 'sum'])):
            draw_lineplot(ax,
                          ss.decimate(electrode.data[name], q=16,
                                      zero_phase=True),
                          dt=network.dt * 16,
                          T=(500, 1000),
                          scaling_factor=1.,
                          vlimround=None,
                          label=label,
                          scalebar=True,
                          unit='mV',
                          ylabels=True if i == 0 else False,
                          color='C{}'.format(i),
                          ztransform=True
                          )
            ax.set_title(label)
        fig.savefig(os.path.join(OUTPUTPATH, 'extracellular_potential.pdf'),
                    bbox_inches='tight')
        plt.close(fig)

        # current-dipole moments, E and I contributions, sum
        fig, axes = plt.subplots(3, 3, figsize=(6.4, 4.8))
        fig.subplots_adjust(wspace=0.45)
        fig.suptitle('current-dipole moments')
        for i, u in enumerate(['x', 'y', 'z']):
            for j, (name, label) in enumerate(zip(['E', 'I', 'imem'],
                                                  ['E', 'I', 'sum'])):
                t = np.arange(current_dipole_moment.data.shape[1]) * network.dt
                inds = (t >= 500) & (t <= 1000)
                axes[i, j].plot(
                    t[inds][::16],
                    ss.decimate(current_dipole_moment.data[name][i, inds],
                                q=16, zero_phase=True),
                    'C{}'.format(j))

                if j == 0:
                    axes[i, j].set_ylabel(r'$\mathbf{p}\cdot\mathbf{e}_{' +
                                          '{}'.format(u) + '}$ (nA$\\mu$m)')
                if i == 0:
                    axes[i, j].set_title(label)
                if i != 2:
                    axes[i, j].set_xticklabels([])
                else:
                    axes[i, j].set_xlabel('t (ms)')
        fig.savefig(os.path.join(OUTPUTPATH, 'current_dipole_moment.pdf'),
                    bbox_inches='tight')
        plt.close(fig)

        # current source density, E and I contributions, sum
        fig, axes = plt.subplots(1, 3, figsize=(6.4, 4.8))
        fig.suptitle('CSD')
        for i, (ax, name, label) in enumerate(zip(axes, ['E', 'I', 'imem'],
                                                  ['E', 'I', 'sum'])):
            draw_lineplot(ax,
                          ss.decimate(csd.data[name], q=16,
                                      zero_phase=True),
                          dt=network.dt * 16,
                          T=(500, 1000),
                          scaling_factor=1.,
                          vlimround=None,
                          label=label,
                          scalebar=True,
                          unit='nA/µm$^3$',
                          ylabels=True if i == 0 else False,
                          color='C{}'.format(i),
                          ztransform=True
                          )
            ax.set_title(label)
        fig.savefig(os.path.join(OUTPUTPATH, 'csd.pdf'),
                    bbox_inches='tight')
        plt.close(fig)

    # population illustration (per RANK)
    if RANK == 0:
        try:
            fig = plt.figure(figsize=(6.4, 4.8 * 2))
            ax = fig.add_subplot(111, projection='3d')
            ax.view_init(elev=5)
            ax.plot(electrode.x, electrode.y, electrode.z, 'ko', zorder=0)
            for i, (name, pop) in enumerate(network.populations.items()):
                for cell in pop.cells:
                    c = 'C0' if name == 'E' else 'C1'
                    for x, y, z, d in zip(cell.x, cell.y, cell.z, cell.d):
                        ax.plot(x, y, z, c, lw=d / 2)
            ax.set_xlabel(r'$x$ ($\mu$m)')
            ax.set_ylabel(r'$y$ ($\mu$m)')
            ax.set_zlabel(r'$z$ ($\mu$m)')
            ax.set_title('network populations')
            fig.savefig(os.path.join(OUTPUTPATH,
                                     'population_RANK_{}.pdf'.format(RANK)),
                        bbox_inches='tight')
        except Exception:
            pass
        plt.close(fig)

    ##########################################################################
    # customary cleanup of object references - the psection() function may not
    # write correct information if NEURON still has object references in memory
    # even if Python references has been deleted. It will also allow the script
    # to be run in successive fashion.
    ##########################################################################
    network.pc.gid_clear()  # allows assigning new gids to threads
    electrode = None
    current_dipole_moment = None
    csd = None
    syn = None
    synapseModel = None
    for population in network.populations.values():
        for cell in population.cells:
            cell.__del__()
            cell = None
        population.cells = None
        population = None
    pop = None
    network = None
    neuron.h('forall delete_section()')
