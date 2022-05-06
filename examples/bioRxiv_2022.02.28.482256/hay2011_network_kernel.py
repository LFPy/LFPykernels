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
import numpy as np
import scipy.signal as ss
import h5py
from parameters import ParameterSet, ParameterSpace
import json
import hashlib
from mpi4py import MPI
from time import time
import neuron
from LFPy import Network, Synapse, RecExtElectrode, CurrentDipoleMoment, \
    LaminarCurrentSourceDensity
from plotting import draw_lineplot
import hay2011_network_parameters as params
import example_network_methods as methods

# set up MPI variables:
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

# avoid same sequence of random numbers from numpy and neuron on each RANK,
# e.g., in order to draw unique cell and synapse locations and random synapse
# activation times
GLOBALSEED = 1234
np.random.seed(GLOBALSEED + RANK)


# class Network parameters:
networkParameters = dict(**params.networkParameters)
networkParameters.update(dict(tstop=600.))

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
    biophys = pset['biophys']
    t_E = np.array([pset['t_E']])
    t_I = np.array([pset['t_I']])
    i_syn = True  # always used
    n_ext = pset['n_ext']
    g_eff = pset['g_eff']
    perseg_Vrest = pset['perseg_Vrest']

    TRANSIENT = 2000

    ##########################################################################
    # Set up shared and population-specific parameters
    ##########################################################################
    # relative path for simulation output:
    PS0 = ParameterSpace('hay2011_PS0.txt')
    pset_0 = ParameterSet(dict(weight_EE=weight_EE,
                               weight_IE=weight_IE,
                               weight_EI=weight_EI,
                               weight_II=weight_II,
                               weight_scaling=weight_scaling,
                               n_ext=PS0['n_ext'].value))
    js_0 = json.dumps(pset_0, sort_keys=True).encode()
    md5_0 = hashlib.md5(js_0).hexdigest()
    OUTPUTPATH_REAL = os.path.join('output', md5_0)

    OUTPUTPATH = os.path.join('output', md5)

    ##########################################################################
    # Parameters dependent on command line input
    ##########################################################################

    with h5py.File(os.path.join(OUTPUTPATH_REAL, 'synapse_connections.h5'
                                ), 'r') as f:
        synapseParameters = []
        for pre in params.population_names:
            synapseParameters.append([])
            for post in params.population_names:
                d = dict()
                for key, value in f['synparams']['{}:{}'.format(pre, post)
                                                 ].items():
                    if key == 'mechanism':
                        d['syntype'] = value[()]
                    else:
                        d[key] = value[()]

                synapseParameters[-1].append(d)

    # create directory for output:
    if RANK == 0:
        if not os.path.isdir(OUTPUTPATH):
            os.mkdir(OUTPUTPATH)
        # remove old simulation output if directory exist
        else:
            for fname in os.listdir(OUTPUTPATH):
                os.unlink(os.path.join(OUTPUTPATH, fname))
    COMM.Barrier()

    # compute population firing rates as this needed in order to incorporate
    # the effective membrane time constant in case of current-based synapses
    if i_syn:
        if RANK == 0:
            nu_X = methods.compute_mean_nu_X(params, OUTPUTPATH_REAL,
                                             TRANSIENT=TRANSIENT)
        else:
            nu_X = None
        nu_X = COMM.bcast(nu_X, root=0)

    # tic tac
    tac = time()
    if RANK == 0:
        with open(os.path.join(OUTPUTPATH, 'tic_tac.txt'), 'w') as f:
            f.write(f'step time_s\nsetup {tac - tic}\n')

    # instantiate Network:
    network = Network(OUTPUTPATH=OUTPUTPATH, **networkParameters)

    # create E and I populations:
    Vrest = dict()
    for j, (name, size) in enumerate(zip(params.population_names,
                                         params.population_sizes)):
        # Extract median soma voltages from actual network simulation and
        # assume this value corresponds to Vrest.
        if RANK == 0:
            if perseg_Vrest:
                with h5py.File(os.path.join(OUTPUTPATH_REAL, 'vmem.h5'
                                            ), 'r') as f:
                    vrest_ = np.median(f[name][()][:, TRANSIENT:], axis=-1)
            else:
                with h5py.File(os.path.join(OUTPUTPATH_REAL, 'somav.h5'
                                            ), 'r') as f:
                    vrest_ = np.median(f[name][()][:, TRANSIENT:])
        else:
            vrest_ = None
        Vrest[name] = COMM.bcast(vrest_, root=0)

        # update cell parameters
        cellParameters = deepcopy(params.cellParameters[name])
        if biophys == 'frozen':
            if name == 'E':
                cellParameters.update({
                    'templatefile': [
                        'L5bPCmodelsEH/models/L5PCbiophys3_frozen.hoc',
                        'L5bPCmodelsEH/models/L5PCtemplate_frozen.hoc'
                    ],
                    'templatename': 'L5PCtemplate_frozen',
                    'custom_fun': [
                        methods.set_V_R,
                        methods.make_cell_uniform
                    ],
                    'custom_fun_args': [dict(Vrest=Vrest[name])] * 2,
                })
            elif name == 'I':
                cellParameters.update({
                    'custom_fun': [
                        methods.set_frozen_hay2011,
                        methods.make_cell_uniform
                    ],
                    'custom_fun_args': [dict(Vrest=Vrest[name])] * 2,
                })
            else:
                raise Exception(f'population {name} not recognized')
        elif biophys == 'lin':
            if name == 'E':
                cellParameters.update({
                    'templatefile': [
                        'L5bPCmodelsEH/models/L5PCbiophys3_lin.hoc',
                        'L5bPCmodelsEH/models/L5PCtemplate_lin.hoc'
                    ],
                    'templatename': 'L5PCtemplate_lin',
                    'custom_fun': [
                        methods.set_V_R,
                        methods.make_cell_uniform],
                    'custom_fun_args': [dict(Vrest=Vrest[name])] * 2,
                })
            elif name == 'I':
                cellParameters.update({
                    'custom_fun': [
                        methods.set_Ih_linearized_hay2011,
                        methods.make_cell_uniform
                    ],
                    'custom_fun_args': [dict(Vrest=Vrest[name])] * 2,
                })
            else:
                raise Exception(f'population {name} not recognized')
        elif biophys == 'pas':
            if name == 'E':
                cellParameters.update({
                    'templatefile': [
                        'L5bPCmodelsEH/models/L5PCbiophys3_pas.hoc',
                        'L5bPCmodelsEH/models/L5PCtemplate_pas.hoc'
                    ],
                    'templatename': 'L5PCtemplate_pas',
                    'custom_fun': [
                        methods.make_cell_uniform],
                    'custom_fun_args': [dict(Vrest=Vrest[name])],
                })
            elif name == 'I':
                cellParameters.update({
                    'custom_fun': [
                        methods.set_pas_hay2011,
                        methods.make_cell_uniform
                    ],
                    'custom_fun_args': [dict(Vrest=Vrest[name])] * 2,
                })
            else:
                raise Exception(f'population {name} not recognized')
        else:
            raise NotImplementedError(f'biophys={biophys} not implemented')

        popParams = deepcopy(params.populationParameters)
        popParams['rotation_args'] = deepcopy(params.rotation_args[name])
        popParams['cell_args'] = cellParameters

        network.create_population(name=name, POP_SIZE=size,
                                  **popParams)

        # in case of current-based synapses, we need lists of segment
        # references for each cell in order to shift g_pas per segment
        if i_syn:
            for cell in network.populations[name].cells:
                cell.allseglist = neuron.h.List()
                for sec in cell.allseclist:
                    for seg in sec:
                        cell.allseglist.append(seg)

        # shared parameters for extrinsic input
        extPar = params.extSynapseParameters.copy()
        # modify parameters for current_based synapses
        if i_syn:
            # extPar['weight'] = - extPar['weight'] * (Vrest[name] - extPar['e'])
            # del extPar['e']  # no longer needed
            try:
                extPar['syntype'] = extPar['syntype'] + 'I'
            except TypeError:  # fix w. python > 3.6
                extPar['syntype'] = extPar['syntype'].decode() + 'I'

        if n_ext[j] > 0:
            # create excitatory background synaptic activity for each cell
            # with Poisson statistics
            for cell in network.populations[name].cells:
                idx = cell.get_rand_idx_area_norm(section='allsec',
                                                  nidx=n_ext[j])
                for i in idx:
                    '''
                    syn = Synapse(cell=cell, idx=i,
                                  **extPar)
                    syn.set_spike_times_w_netstim(
                        interval=params.netstim_interval,
                        seed=np.random.rand() * 2**32 - 1)
                    '''
                    _ = np.random.rand() * 2**32 - 1  # throw away to ensure
                                                      # similar sequences

                    if i_syn & g_eff:
                        # compute and apply shift of seg.g_pas:
                        if extPar['syntype'] == 'Exp2SynI':
                            # compute area under temporal kernel (ms)
                            beta = methods.integrate_beta(extPar['tau1'],
                                                          extPar['tau2'])
                            g_shift = (
                                # uS
                                abs(params.extSynapseParameters['weight'])
                                / cell.area[i]  # um**2
                                * beta  # ms
                                / (params.netstim_interval)  # ms
                            )  # (hekto-S/cm**2)
                            try:
                                cell.allseglist[i].g_pas += (
                                    g_shift * 100.)  # unit: S/cm**2
                            except AttributeError:
                                cell.allseglist[i].gl_hh += (
                                    g_shift * 100.)  # unit: S/cm**2

                        else:
                            errmsg = '{} not supported'.format(d['syntype'])
                            raise NotImplementedError(errmsg)

    # tic tac
    tic = time()
    if RANK == 0:
        with open(os.path.join(OUTPUTPATH, 'tic_tac.txt'), 'a') as f:
            f.write(f'create {tic - tac}\n')

    # recreate synapses of recurrent network using spike times of network
    # for activation times - a.k.a. hybrid scheme
    with h5py.File(os.path.join(OUTPUTPATH_REAL, 'synapse_connections.h5'), 'r'
                   ) as f:
        with h5py.File(os.path.join(OUTPUTPATH_REAL, 'spikes.h5'), 'r') as fs:
            for i, pre in enumerate(params.population_names):
                gids = fs[pre]['gids'][()]
                times = fs[pre]['times'][()]
                for j, post in enumerate(params.population_names):
                    conn = f['{}:{}'.format(pre, post)][()]
                    # synparams = f['synparams']['{}:{}'.format(pre, post)]
                    for cell in network.populations[post].cells:
                        rows = conn['gid'] == cell.gid
                        for c in conn[rows]:
                            if cell.templatename.rfind('L5PCtemplate') >= 0:
                                _section = '['.join(
                                    [cell.templatename] +
                                    c['sec'].decode().split('[')[1:])
                                idxs = cell.get_idx(_section)
                            else:
                                idxs = cell.get_idx(c['sec'].decode())
                            idxx = np.linspace(1 / idxs.size / 2,
                                               1 - 1 / idxs.size / 2,
                                               idxs.size)
                            idx = idxs[np.argmin((c['sec.x'] - idxx)**2)]

                            # make a copy of Synapse kwargs as it is modified
                            # for current based synapses
                            d = synapseParameters[i][j].copy()
                            if i_syn:
                                d['weight'] = - c['weight'] * (Vrest[post]
                                                               - d['e'])
                                del d['e']  # no longer needed
                                try:
                                    d['syntype'] = d['syntype'] + 'I'
                                except TypeError:  # fix w. python > 3.6
                                    d['syntype'] = d['syntype'].decode() + 'I'

                                # compute and apply shift of seg.g_pas:
                                if g_eff:
                                    if d['syntype'] == 'Exp2SynI':
                                        # area under temporal kernel (ms)
                                        beta = methods.integrate_beta(
                                            d['tau1'], d['tau2'])
                                        g_shift = (abs(c['weight'])
                                                   / cell.area[idx]
                                                   * beta
                                                   * nu_X[pre])  # deci-S/cm**2
                                        try:
                                            cell.allseglist[idx].g_pas += (
                                                g_shift * 0.1)  # unit: S/cm**2
                                        except AttributeError:
                                            cell.allseglist[idx].gl_hh += (
                                                g_shift * 0.1)  # unit: S/cm**2
                                    else:
                                        errmsg = '{} not supported'.format(
                                            d['syntype'])
                                        raise NotImplementedError(errmsg)
                            else:
                                d['weight'] = c['weight']
                                if isinstance(d['syntype'], bytes):
                                    d['syntype'] = d['syntype'].decode()

                            # create Synapse
                            syn = Synapse(cell=cell,
                                          idx=idx,
                                          **d)

                            # set activation times
                            syn.set_spike_times(
                                c['delay'] + t_E
                                if pre == 'E'
                                else
                                c['delay'] + t_I)

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
        # probes = [],
        probes=[electrode, current_dipole_moment, csd],
        **params.networkSimulationArguments
    )

    # tic tac
    tic = time()
    if RANK == 0:
        with open(os.path.join(OUTPUTPATH, 'tic_tac.txt'), 'a') as f:
            f.write(f'simulate {tic - tac}\n')

    if RANK == 0:
        # save spikes
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
        # extracellular potentials, E and I contributions, sum
        fig, axes = plt.subplots(1, 3, figsize=(6.4, 4.8))
        fig.suptitle('extracellular potentials')
        for i, (ax, name, label) in enumerate(zip(axes, ['E', 'I', 'imem'],
                                                  ['E', 'I', 'sum'])):
            draw_lineplot(ax,
                          ss.decimate(electrode.data[name], q=16,
                                      zero_phase=True),
                          dt=network.dt * 16,
                          T=(100, 500),
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
                inds = (t >= 100) & (t <= 500)
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

    ##########################################################################
    # customary cleanup of object references - the psection() function may not
    # write correct information if NEURON still has object references in memory
    # even if Python references has been deleted. It will also allow the script
    # to be run in successive fashion.
    ##########################################################################
    network.pc.gid_clear()  # allows assigning new gids to threads
    electrode = None
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
