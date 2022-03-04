#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pymoo.factory import get_visualization
from pymoo.factory import get_decision_making, get_reference_directions
from pymoo.visualization.scatter import Scatter
from pymoo.model.callback import Callback
from pymoo.optimize import minimize
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.model.problem import Problem
from nest import raster_plot
import nest
import pandas as pd
import hashlib
import json
from copy import deepcopy
import scipy.stats as st
# import example_network_parameters as params
from example_network_parameters import (networkParameters, population_names,
                                        population_sizes)
import scipy.signal as ss
import h5py
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np
from parameters import ParameterSpace, ParameterSet
import os
from time import time
# get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


# In[ ]:


nest.set_verbosity('M_ERROR')


# In[ ]:


PS0 = ParameterSpace('PS0.txt')
PS1 = ParameterSpace('PS1.txt')
# PS2 = ParameterSpace('PS2.txt')


# In[ ]:


TRANSIENT = 2000
dt = networkParameters['dt']
tstop = networkParameters['tstop']
tau = 100  # time lag relative to spike for kernel predictions


# In[ ]:


# plt.mlab.psd/csd settings
Fs = 1000 / dt
NFFT = 1024 * 2
noverlap = 768 * 2
detrend = 'constant'


# In[ ]:


# figure out which real LFP to compare with
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


# In[ ]:


def get_spike_rate(times):
    bins = (np.arange(TRANSIENT / dt, tstop / dt + 1)
            * dt - dt / 2)
    hist, _ = np.histogram(times, bins=bins)
    return bins, hist.astype(float) / dt


# In[ ]:


def get_mean_spike_rate(times):
    times = times[times >= TRANSIENT]
    return times.size / (tstop - TRANSIENT) * 1000


# In[ ]:


def get_psd(nu, cutoff=200):
    freqs, Pxx = ss.welch(nu, fs=Fs, nperseg=NFFT,
                          noverlap=noverlap, detrend=detrend)
    return freqs[freqs <= cutoff], Pxx[freqs <= cutoff]


# In[ ]:


def zscore(x):
    return (x - x.mean()) / x.std()


# In[ ]:


# mean firing rates, rates, rate spectras of "real" network populations
mean_nu_X = dict()
nu_X = dict()
psd_X = dict()
with h5py.File(os.path.join(OUTPUTPATH_REAL, 'spikes.h5'), 'r') as f:
    for i, X in enumerate(population_names):
        times = np.concatenate(f[X]['times'])
        mean_nu_X[X] = get_mean_spike_rate(times) / population_sizes[i]
        bins, nu_X[X] = get_spike_rate(times)
        freqs, psd_X[X] = get_psd(nu_X[X])

mean_nu_X


# In[ ]:


# plot firing rate time series of "real" network
inds = (bins[:-1] >= TRANSIENT) & (bins[:-1] <= TRANSIENT + 100)
fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(16, 9))
with h5py.File(os.path.join(OUTPUTPATH_REAL, 'spikes.h5'), 'r') as f:
    for i, (X, N_X) in enumerate(zip(population_names, population_sizes)):
        axes[i].step(bins[:-1][inds], nu_X[X][inds])
        axes[i].set_title(r'$\nu_{%s}(t)$' % X)
        axes[i].set_ylabel(r'$\nu_X$ (# spikes/s)')
axes[1].set_xlabel('$t$ (ms)')


# In[ ]:


# firing rate spectra
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
for X in population_names:
    ax.plot(freqs, psd_X[X], label=X)
ax.set_xlabel('$f$ (Hz)', labelpad=0)
ax.set_ylabel(r'PSD$_\nu$ (s$^{-2}$/Hz)')
ax.legend()


# In[ ]:


# create NEST LIF network
class Network(object):
    """Class implementing a LIF network"""

    def __init__(self,
                 X=['E', 'I'],
                 N_X=[8192, 1024],
                 C_m_X=[269., 100.],
                 tau_m_X=[10., 10.],
                 E_L_X=[-65., -65.],
                 C_YX=[[0.5, 0.5], [0.5, 0.5]],
                 J_YX=[[1.5, 1.5], [-20., -15.]],
                 delay_YX=[[3., 3.], [3.0, 3.0]],
                 tau_syn_YX=[[0.5, 0.5], [0.5, 0.5]],
                 n_ext=[450, 160],
                 nu_ext=40.,
                 J_ext=27.6,
                 model='iaf_psc_exp',
                 dt=2**-4,
                 verbose=False,
                 **kwargs
                 ):
        self.X = X
        self.N_X = N_X
        self.C_m_X = C_m_X
        self.tau_m_X = tau_m_X
        self.E_L_X = E_L_X
        self.C_YX = C_YX
        self.J_YX = J_YX
        self.delay_YX = delay_YX
        self.tau_syn_YX = tau_syn_YX
        self.n_ext = n_ext
        self.nu_ext = nu_ext
        self.J_ext = J_ext
        self.model = model
        self.dt = dt
        self.verbose = verbose

        self._create()
        self._connect()

    def _create(self):
        """Create network nodes and connections"""

        nest.ResetKernel()
        nest.SetKernelStatus(
            dict(
                local_num_threads=32,
                resolution=self.dt,
                tics_per_ms=1000 / dt))

        if self.verbose:
            print('creating...')
        # neurons
        self.neurons = {}
        for (
                X,
                N,
                C_m,
                tau_m,
                E_L,
                (tau_syn_ex,
                 tau_syn_in)) in zip(
                self.X,
                self.N_X,
                self.C_m_X,
                self.tau_m_X,
                self.E_L_X,
                self.tau_syn_YX):
            params = dict(
                C_m=C_m,
                tau_m=tau_m,
                E_L=E_L,
                V_reset=E_L,
                tau_syn_ex=tau_syn_ex,
                tau_syn_in=tau_syn_in
            )
            self.neurons[X] = nest.Create(self.model, N, params)

        # poisson generators
        self.poisson = {}
        for X, n_ext in zip(self.X, self.n_ext):
            self.poisson[X] = nest.Create(
                'poisson_generator', 1, dict(
                    rate=self.nu_ext * n_ext))

        # spike recorders
        self.spike_recorders = {}
        for X in self.X:
            self.spike_recorders[X] = nest.Create('spike_recorder', 1)

    def _connect(self):
        if self.verbose:
            print('connecting...')
        for i, X in enumerate(self.X):
            # recurrent connections
            for j, Y in enumerate(self.X):
                if self.verbose:
                    print(f'connecting {X} to {Y}...')
                conn_spec = dict(
                    rule='pairwise_bernoulli',
                    p=self.C_YX[i][j],
                )
                syn_spec = dict(
                    synapse_model='static_synapse',
                    weight=nest.math.redraw(
                        nest.random.normal(
                            mean=self.J_YX[i][j],
                            std=abs(self.J_YX[i][j]) * 0.1,
                        ),
                        min=0. if self.J_YX[i][j] >= 0 else np.NINF,
                        max=np.Inf if self.J_YX[i][j] >= 0 else 0.,
                    ),

                    delay=nest.math.redraw(
                        nest.random.normal(
                            mean=self.delay_YX[i][j],
                            std=self.delay_YX[i][j] * 0.5,
                        ),
                        min=0.3,
                        max=np.Inf,
                    )
                )

                nest.Connect(
                    self.neurons[X],
                    self.neurons[Y],
                    conn_spec,
                    syn_spec)

            # poisson generators
            if self.verbose:
                print(f'connecting poisson_generator[{X}] to {X}...')
            nest.Connect(
                self.poisson[X],
                self.neurons[X],
                'all_to_all',
                dict(
                    weight=self.J_ext))

            # recorders
            if self.verbose:
                print(f'connecting spike_recorder[{X}] to {X}...')
            nest.Connect(self.neurons[X], self.spike_recorders[X])

    def simulate(self, tstop=6000):
        """Instantiate and run simulation"""
        if self.verbose:
            print('simulating...')
        nest.Simulate(tstop)
        if self.verbose:
            print('done!')


# In[ ]:


'''
# test Network class
tic = time()
net = Network()
net.simulate(tstop=tstop + dt)
toc = time()
print('simulation time:', toc-tic)
'''

'''
# In[ ]:


# mean firing rates of "real" network populations
lif_mean_nu_X = dict()  # mean spike rates
lif_nu_X = dict()  # binned firing rate
lif_psd_X = dict()
for i, X in enumerate(population_names):
    times = nest.GetStatus(net.spike_recorders[X])[0]['events']['times']
    times = times[times >= TRANSIENT]

    lif_mean_nu_X[X] = get_mean_spike_rate(times)
    _, lif_nu_X[X] = get_spike_rate(times)
    _, lif_psd_X[X] = get_psd(lif_nu_X[X])

lif_mean_nu_X


# In[ ]:


# power spectra
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
for X in population_names:
    ax.plot(freqs, psd_X[X], label=X)
    ax.plot(freqs, lif_psd_X[X], label=X)
ax.set_xlabel('$f$ (Hz)', labelpad=0)
ax.set_ylabel(r'PSD$_\nu$ (s$^{-2}$/Hz)')
ax.legend()


# In[ ]:
'''


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
    )
    net = Network(**params)
    net.simulate(tstop=tstop + dt)

    psd_X = []
    nu_X = []
    mean_nu_X = []
    for i, X in enumerate(population_names):
        times = nest.GetStatus(net.spike_recorders[X])[0]['events']['times']
        times = times[times >= TRANSIENT]

        mean_nu_X += [get_mean_spike_rate(times) / population_sizes[i]]
        _, lif_nu_X = get_spike_rate(times)
        nu_X += [lif_nu_X]
        _, lif_psd_X = get_psd(lif_nu_X)
        psd_X += [lif_psd_X]

    return mean_nu_X, nu_X, psd_X


'''
# In[ ]:


lif_mean_nu_X, lif_nu_X, lif_psd_X = get_lif_data(
    x=[1.5, 1.5, -20., -15., 27.6, 3.0, 3.0, 3.0, 3.0])


# In[ ]:
'''


class LifNet(Problem):
    """inherited pymoo Problem class required for minimizer"""

    def __init__(self):
        super().__init__(n_var=11,
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


# In[ ]:


class MyCallback(Callback):

    def __init__(self) -> None:
        super().__init__()
        self.data["F"] = []
        self.data["X"] = []

    def notify(self, algorithm):
        self.data["F"].append(algorithm.pop.get("F"))
        self.data["X"].append(algorithm.pop.get("X"))


# In[ ]:


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


# In[ ]:


# Find a suitable parameter combination performing well with all features, from
# https://pymoo.org/decision_making/index.html#Pseudo-Weights
ref_dirs = get_reference_directions("das-dennis", 4, n_partitions=12)
F = res.F

weights = np.array([0.25, 0.25, 0.25, 0.25])
a, pseudo_weights = get_decision_making(
    "pseudo-weights", weights).do(res.F, return_pseudo_weights=True)
'''
plot = get_visualization("petal", bounds=(0, 0.5), reverse=True)
plot.add(res.F[[a]])
plot.show()
'''

# In[ ]:

'''
res.F[[a]], res.X[[a]][0]
'''

# In[ ]:

'''
# plot location in feature space
plot = Scatter(figsize=(12, 12))
for F in res.algorithm.callback.data["F"]:
    plot.add(F, color="gray")
plot.add(res.F, color="black")
plot.add(res.F[[a]], color="red")
plot.show()
'''

# In[ ]:


# Store data to HDF5
with h5py.File('Fit_LIF_net.h5', 'w') as f:
    f["F"] = np.r_[res.algorithm.callback.data["F"]]
    f["X"] = np.r_[res.algorithm.callback.data["X"]]
    f["F_opt"] = res.F[[a]]
    f["X_opt"] = res.X[[a]][0]


# In[ ]:


# pymoo messes up MPL rcParams
# get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


# simulate network using updated params
lif_mean_nu_X, lif_nu_X, lif_psd_X = get_lif_data(x=res.X[[a]][0])


# In[ ]:


# compare power spectra with updated params
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
for i, X in enumerate(population_names):
    ax.semilogy(freqs, psd_X[X], label=f'{X} mc {mean_nu_X[X]} s-1')
    ax.semilogy(freqs, lif_psd_X[i], label=f'{X} lif {lif_mean_nu_X[i]} s-1')
ax.set_xlabel('$f$ (Hz)', labelpad=0)
ax.set_ylabel(r'PSD$_\nu$ (s$^{-2}$/Hz)')
ax.legend()
fig.savefig('Fit_LIF_net_PSD.pdf', bbox_inches='tight')


# In[ ]:


# compare firing rate auto- and cross-correlation functions:
maxlag = 100  # ms
lag = (np.arange(nu_X['E'].size) - nu_X['E'].size // 2) * dt
lag_inds = (lag >= -maxlag) & (lag <= maxlag)
fig, axes = plt.subplots(2, 2, figsize=(16, 8), sharex=True, sharey=True)
for i, X in enumerate(population_names):
    for j, Y in enumerate(population_names):
        xcorr_mc = np.correlate(zscore(nu_X[X]),
                                zscore(nu_X[Y]), 'same') / zscore(nu_X[X]).size
        xcorr_lif = np.correlate(zscore(lif_nu_X[i]), zscore(
            lif_nu_X[j]), 'same') / zscore(lif_nu_X[i]).size

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
fig.savefig('Fit_LIF_net_CXY.pdf', bbox_inches='tight')


# In[ ]:
