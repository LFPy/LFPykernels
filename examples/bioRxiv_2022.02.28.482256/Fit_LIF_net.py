#!/usr/bin/env python
# coding: utf-8
from pymoo.factory import get_decision_making, get_reference_directions
from pymoo.model.callback import Callback
from pymoo.optimize import minimize
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.model.problem import Problem
import nest
import hashlib
import json
from example_network_parameters import (networkParameters, population_names,
                                        population_sizes)
import example_network_methods as methods
import scipy.signal as ss
import h5py
import matplotlib.pyplot as plt
import numpy as np
from parameters import ParameterSpace, ParameterSet
import os
import LIF_net


# print only NEST errors
nest.set_verbosity('M_ERROR')


# Load parameter spaces
PS0 = ParameterSpace('PS0.txt')
PS1 = ParameterSpace('PS1.txt')

# simulation parameters
TRANSIENT = 2000
dt = networkParameters['dt']
tstop = networkParameters['tstop']


# scipy.signal.psd settings for power spectra
Fs = 1000 / dt
NFFT = 2048
noverlap = 1536
detrend = 'constant'
cutoff = 200

# figure out which ground truth data to compare with
for pset in PS1.iter_inner():
    weight_EE = pset['weight_EE']
    weight_IE = pset['weight_IE']
    weight_EI = pset['weight_EI']
    weight_II = pset['weight_II']
    weight_scaling = pset['weight_scaling']
    pset_0 = ParameterSet(dict(weight_EE=weight_EE,
                               weight_IE=weight_IE,
                               weight_EI=weight_EI,
                               weight_II=weight_II,
                               weight_scaling=weight_scaling,
                               n_ext=PS0['n_ext'].value))
    js_0 = json.dumps(pset_0, sort_keys=True).encode()
    md5_0 = hashlib.md5(js_0).hexdigest()
    OUTPUTPATH_REAL = os.path.join('output', md5_0)

    break
print(f'comparing with ground truth dataset: {OUTPUTPATH_REAL}')


# mean firing rates, rates, rate spectras of "real" network populations
mean_nu_X = dict()
nu_X = dict()
psd_X = dict()
with h5py.File(os.path.join(OUTPUTPATH_REAL, 'spikes.h5'), 'r') as f:
    for i, X in enumerate(population_names):
        times = np.concatenate(f[X]['times'])
        mean_nu_X[X] = LIF_net.get_mean_spike_rate(
            times, TRANSIENT=TRANSIENT, tstop=tstop) / population_sizes[i]
        bins, nu_X[X] = LIF_net.get_spike_rate(
            times, TRANSIENT=TRANSIENT, tstop=tstop, dt=dt)
        freqs, psd_X[X] = LIF_net.get_psd(
            nu_X[X], Fs=Fs, NFFT=NFFT, noverlap=noverlap, detrend=detrend,
            cutoff=cutoff)


def get_lif_data(x=[1.5, 1.5, -20., -15., 27.6, 3.0, 3.0, 3.0, 3.0, 269, 100]):
    params = dict(
        X=['E', 'I'],
        N_X=[8192, 1024],
        C_m_X=[x[9], x[10]],  # [269., 100.],
        tau_m_X=[10., 10.],
        E_L_X=[-65., -65.],
        C_YX=[[0.5, 0.5], [0.5, 0.5]],
        J_YX=[[x[0], x[1]], [x[2], x[3]]],
        delay_YX=[[x[5], x[6]], [x[7], x[8]]],
        tau_syn_YX=[[0.5, 0.5], [0.5, 0.5]],
        n_ext=PS0['n_ext'].value,
        nu_ext=40.,
        # J_ext=27.6,
        J_ext=x[4],
        model='iaf_psc_exp',
        dt=2**-4,
        local_num_threads=32,
    )
    net = LIF_net.Network(**params)
    net.simulate(tstop=tstop + dt)

    psd_X = []
    nu_X = []
    mean_nu_X = []
    for i, X in enumerate(population_names):
        times = nest.GetStatus(net.spike_recorders[X])[0]['events']['times']
        times = times[times >= TRANSIENT]

        mean_nu_X += [LIF_net.get_mean_spike_rate(times,
                                          TRANSIENT=TRANSIENT,
                                          tstop=tstop) / population_sizes[i]]
        _, lif_nu_X = LIF_net.get_spike_rate(
            times, TRANSIENT=TRANSIENT, tstop=tstop, dt=dt)
        nu_X += [lif_nu_X]
        _, lif_psd_X = LIF_net.get_psd(
            lif_nu_X, Fs=Fs, NFFT=NFFT, noverlap=noverlap, detrend=detrend,
            cutoff=cutoff)
        psd_X += [lif_psd_X]

    return mean_nu_X, nu_X, psd_X


class LifNet(Problem):
    """inherited pymoo Problem class required for minimizer"""

    def __init__(self):
        super().__init__(
            n_var=11,
            n_obj=4,
            n_constr=0,
            xl=[1.1, 1.5, -25., -14., 28, 1.0, 1.0, 1.0, 1.0, 270, 100],
            xu=[1.8, 2.1, -18., -8., 32, 4.0, 4.0, 4.0, 4.0, 310, 120])

    def _evaluate(self, x, out, *args, **kwargs):
        f0 = []
        f1 = []
        f2 = []
        f3 = []
        for xi in x:
            lif_mean_nu_X, lif_nu_X, lif_psd_X = get_lif_data(xi)

            # RMSE between mean rates
            f0 += [np.sqrt((mean_nu_X['E'] - lif_mean_nu_X[0])**2)]
            f1 += [np.sqrt((mean_nu_X['I'] - lif_mean_nu_X[1])**2)]

            # RMSE between rate spectras
            f2 += [np.sqrt((psd_X['E'] - lif_psd_X[0])**2).sum()]
            f3 += [np.sqrt((psd_X['I'] - lif_psd_X[1])**2).sum()]

        out['F'] = np.column_stack([f0, f1, f2, f3])


class MyCallback(Callback):

    def __init__(self) -> None:
        super().__init__()
        self.data["F"] = []
        self.data["X"] = []

    def notify(self, algorithm):
        self.data["F"].append(algorithm.pop.get("F"))
        self.data["X"].append(algorithm.pop.get("X"))


# minimize using NSGA-II multiobjective evolutionary algorithm
problem = LifNet()
algorithm = NSGA2(pop_size=100)
res = minimize(problem,
               algorithm,
               ('n_gen', 20),
               seed=1,
               save_history=True,
               callback=MyCallback(),
               verbose=True)


# Find a suitable parameter combination performing well with all features, from
# https://pymoo.org/decision_making/index.html#Pseudo-Weights
ref_dirs = get_reference_directions("das-dennis", 4, n_partitions=12)
F = res.F

weights = np.array([0.25, 0.25, 0.25, 0.25])
a, pseudo_weights = get_decision_making(
    "pseudo-weights", weights).do(res.F, return_pseudo_weights=True)


# Store data to HDF5
with h5py.File('Fit_LIF_net.h5', 'w') as f:
    f["F"] = np.r_[res.algorithm.callback.data["F"]]
    f["X"] = np.r_[res.algorithm.callback.data["X"]]
    f["F_opt"] = res.F[[a]]
    f["X_opt"] = res.X[[a]][0]


# simulate network using updated params
lif_mean_nu_X, lif_nu_X, lif_psd_X = get_lif_data(x=res.X[[a]][0])


# compare power spectra with updated params
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
for i, X in enumerate(population_names):
    ax.semilogy(freqs, psd_X[X], label=f'{X} mc {mean_nu_X[X]} s-1')
    ax.semilogy(freqs, lif_psd_X[i], label=f'{X} lif {lif_mean_nu_X[i]} s-1')
ax.set_xlabel('$f$ (Hz)', labelpad=0)
ax.set_ylabel(r'PSD$_\nu$ (s$^{-2}$/Hz)')
ax.legend()
fig.savefig(os.path.join('figures', 'Fit_LIF_net_PSD.pdf'),
            bbox_inches='tight')


# compare firing rate auto- and cross-correlation functions:
maxlag = 100  # ms
lag = (np.arange(nu_X['E'].size) - nu_X['E'].size // 2) * dt
lag_inds = (lag >= -maxlag) & (lag <= maxlag)
fig, axes = plt.subplots(2, 2, figsize=(16, 8), sharex=True, sharey=True)
for i, X in enumerate(population_names):
    for j, Y in enumerate(population_names):
        xcorr_mc = np.correlate(methods.zscore(nu_X[X]), methods.zscore(
            nu_X[Y]), 'same') / methods.zscore(nu_X[X]).size
        xcorr_lif = np.correlate(methods.zscore(lif_nu_X[i]), methods.zscore(
            lif_nu_X[j]), 'same') / methods.zscore(lif_nu_X[i]).size

        axes[i, j].plot(lag[lag_inds], xcorr_mc[lag_inds], label='mc')
        axes[i, j].plot(lag[lag_inds], xcorr_lif[lag_inds], label='lif')

        axes[i, j].set_title(f'$C_{{{X}{Y}}}$')

        if j == 0:
            axes[i, j].set_ylabel(r'$C_{XY}$ (-)', labelpad=0)
        if i == 1:
            axes[i, j].set_xlabel(r'$\tau$ (ms)', labelpad=0)

        axes[i, j].axis(axes[i, j].axis('tight'))
        axes[i, j].set_ylim(ymax=0.25)
        axes[i, j].legend()
fig.savefig(os.path.join('figures', 'Fit_LIF_net_CXY.pdf'),
            bbox_inches='tight')
