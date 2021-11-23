#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Class definitions for lfpykernels

Copyright (C) 2021 Computational Neuroscience Group, NMBU.

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
import numpy as np
import scipy.stats as st
import scipy.integrate as si
from copy import deepcopy
from LFPy import TemplateCell, Synapse
import lfpykit
import neuron
from warnings import warn


def integrate_beta(tau_1, tau_2):
    """
    Return the integral of the beta function from 0 to infty

    Parameters
    ----------
    tau_1: float
        rise time constant
    tau_2: float
        decay time constant. tau_2 > tau_1

    Returns
    -------
    float
    """
    tp = (tau_1 * tau_2) / (tau_2 - tau_1) * np.log(tau_2 / tau_1)
    u = np.exp(-tp / tau_2) - np.exp(-tp / tau_1)
    return (tau_2 - tau_1) / u


class KernelApprox(object):
    '''Class for computing linear spike-to-signal filter kernels resulting
    from presynaptic spiking activity and resulting postsynaptic currents

    Parameters
    ----------
    X: list of str
        presynaptic populations
    Y: str
        postsynaptic population
    N_X: array of int
        presynaptic population sizes
    N_Y: int
        postsynaptic population size
    C_YX: array of float
        pairwise connection probabilities betwen populations X and Y
    multapseFunction: callable
        ``scipy.stats.rv_discrete`` or ``scipy.stats.rv_continuous`` like
        function for determining mean number of synapse
        instances per connection between populations X and Y.
        Default is ``scipy.stats.truncnorm``.
    multapseParameters: list of dict
        kwargs for ``multapseFunction``
    cellParameters: dict
        kwargs for ``LFPy.TemplateCell`` class for cell representative of the
        entire postsynaptic population
    populationParameters: dict
        keys: ``radius``, ``loc``, ``scale`` with float values representing
        radius in xy-plane and mean and standard deviation of cell positions
        along the z-axis
    delayFunction: callable
        ``scipy.stats.rv_continuous`` like callable with pdf method
        for delays between presynaptic populations ``X`` and postsynaptic
        population ``Y``.
        Default is ``scipy.stats.truncnorm``.
    delayParameters: list of dict
        kwargs for callable ``delayFunction``
    synapseParameters: list of dict
        kwargs for ``LFPy.Synapse``, assuming conductance based synapse which
        will be linearized to current based synapse for connections between
        populations X and Y
    synapsePositionArguments: list of dict
        kwargs for ``KernelApprox.get_rand_idx_area_and_distribution_prob``
        method for connections between populations X and Y
    extSynapseParameters: dict
        shared parameters for extrinsic synapses distributed homogeneously
        across morphology
    nu_ext: float
        activation rate of extrinsic synapses (1/s)
    n_ext: float
        number of extrinsic synapses
    nu_X: dict of floats
        presynaptic population rates (1/s)
    '''

    def __init__(
            self,
            X=['E'],
            Y='E',
            N_X=np.array([1024]),
            N_Y=1024,
            C_YX=np.array([0.1]),
            cellParameters=dict(),
            populationParameters=dict(radius=100, loc=0, scale=50),
            rotationParameters=dict(x=0., y=0.),
            multapseFunction=st.truncnorm,
            multapseParameters=[dict(a=-0.2, b=1.6, loc=2, scale=5)],
            delayFunction=st.truncnorm,
            delayParameters=[dict(a=-4.0, b=np.inf, loc=1.5, scale=0.3)],
            synapseParameters=[dict(weight=0.001, syntype='Exp2Syn',
                                    tau1=0.2, tau2=1.8, e=0.)],
            synapsePositionArguments=[dict(section=['soma', 'apic'],
                                           fun=[st.norm, st.norm],
                                           funargs=[dict(loc=0., scale=100.),
                                                    dict(loc=750., scale=100.)
                                                    ],
                                           funweights=[0.5, 1.])],
            extSynapseParameters=dict(syntype='Exp2Syn', weight=0.0005,
                                      tau1=0.2, tau2=1.8, e=0.0),
            nu_ext=40.,
            n_ext=128.,
            nu_X=dict(E=1.)
    ):
        # set attributes
        self.X = X
        self.Y = Y
        self.N_X = N_X
        self.N_Y = N_Y
        self.C_YX = C_YX
        self.multapseFunction = multapseFunction
        self.multapseParameters = multapseParameters
        self.cellParameters = cellParameters
        self.populationParameters = populationParameters
        self.rotationParameters = rotationParameters
        self.delayFunction = delayFunction
        self.delayParameters = delayParameters
        self.synapseParameters = synapseParameters
        self.synapsePositionArguments = synapsePositionArguments
        self.extSynapseParameters = extSynapseParameters
        self.nu_ext = nu_ext
        self.n_ext = n_ext
        self.nu_X = nu_X

    def get_delay(self, X, dt, tau):
        '''Get normalized transfer function for conduction delay distribution
        for connections between population X and Y

        Parameters
        ----------
        X: str
            presynaptic population name
        dt: float
            time resolution
        tau: float
            time lag

        Returns
        -------
        h_delta: ndarray
            shape (2 * tau // dt + 1) array with transfer function for delay
            distribution
        '''
        t = np.linspace(-tau, tau, int(2 * tau // dt + 1))
        [i] = np.where(np.array(self.X) == X)[0]
        h_delay = self.delayFunction(**self.delayParameters[i]).pdf(t)
        return h_delay / h_delay.sum()

    def draw_rand_pos(self, SIZE, radius, loc, scale, cap=None):
        """
        Draw ``SIZE`` random locations within radius ``radius`` in xy-plane,
        at mean depth ``loc`` and standard deviation ``scale`` along z-axis.

        Parameters
        ----------
        SIZE: int
            Population size
        radius: float
            Population radius (µm)
        loc: float
            expected mean depth (µm)
        scale: float
            expected standard deviation of depth (µm)
        cap: None, float or length to list of floats
            if float, cap distribution between [loc-cap, loc+cap),
            if list, cap distribution between [loc-cap[0], loc+cap[1]]

        Returns
        -------
        pos: ndarray
            shape (SIZE, 3) ndarray with randomly chosen locations
        """
        pos = np.zeros((SIZE, 3))

        for i in range(SIZE):
            pos[i, 0] = (np.random.rand() - 0.5) * radius * 2
            pos[i, 1] = (np.random.rand() - 0.5) * radius * 2
            while np.sqrt(pos[i, 0]**2 + pos[i, 1]**2) >= radius:
                pos[i, 0] = (np.random.rand() - 0.5) * radius * 2
                pos[i, 1] = (np.random.rand() - 0.5) * radius * 2
        pos[:, 2] = np.random.normal(loc=loc, scale=scale, size=SIZE)
        if cap is not None:
            if type(cap) in [float, np.float32, np.float64]:
                while not np.all((pos[:, 2] >= loc - cap)
                                 & (pos[:, 2] < loc + cap)):
                    ii = (pos[:, 2] < loc - cap) ^ (pos[:, 2] > loc + cap)
                    pos[ii, 2] = np.random.normal(loc=loc, scale=scale,
                                                  size=ii.sum())
            elif isinstance(cap, list):
                assert len(cap) == 2, \
                    'cap = {} is not a length 2 list'.format(float)
                while not np.all((pos[:, 2] >= loc - cap[0])
                                 & (pos[:, 2] < loc + cap[1])):
                    ii = (pos[:, 2] < loc - cap[0]
                          ) ^ (pos[:, 2] > loc + cap[1])
                    pos[ii, 2] = np.random.normal(loc=loc, scale=scale,
                                                  size=ii.sum())
            else:
                raise Exception('cap = {} is not None'.format(float),
                                'a float or length 2 list of floats')

        return pos

    def get_rand_idx_area_and_distribution_prob(self, cell, section='allsec',
                                                z_min=-1E6, z_max=1E6,
                                                fun=st.norm,
                                                funargs=dict(loc=0, scale=100),
                                                funweights=None):
        """
        Return probability
        normalized to the membrane area of each segment multiplied by
        the value of the probability density function of ``fun``, a function
        in the ``scipy.stats`` module with corresponding function arguments
        in ``funargs`` on the interval [z_min, z_max]

        Parameters
        ----------
        section: str
            string matching a section name
        z_min: float
            lower depth interval
        z_max: float
            upper depth interval
        fun: function or str, or iterable of function or str
            if function a scipy.stats method, if str, must be method in
            scipy.stats module with the same name (like ``norm``),
            if iterable (list, tuple, numpy.array) of function or str some
            probability distribution in scipy.stats module
        funargs: dict or iterable
            iterable (list, tuple, numpy.array) of dict, arguments to fun.pdf
            method (e.g., w. keys ``loc`` and ``scale``)
        funweights: None or iterable
            iterable (list, tuple, numpy.array) of floats, scaling of each
            individual fun (i.e., introduces layer specificity)
        """
        poss_idx = cell.get_idx(section=section, z_min=z_min, z_max=z_max)
        if poss_idx.size == 0:
            print('No possible segment idx match query - returning '
                  'empty array')
            return np.array([])
        else:
            p = np.zeros_like(cell.area)
            p[poss_idx] = cell.area[poss_idx]

            # scale with density function
            if type(fun) in [list, tuple, np.ndarray]:
                assert type(funargs) in [list, tuple, np.ndarray]
                assert type(funweights) in [list, tuple, np.ndarray]
                assert len(fun) == len(funargs) & len(fun) == len(funweights)
                mod = np.zeros(poss_idx.shape)
                for f, args, scl in zip(fun, funargs, funweights):
                    if isinstance(f, str) and f in dir(st):
                        f = getattr(st, f)
                    df = f(**args)
                    mod += df.pdf(x=cell.z[poss_idx].mean(axis=-1)) * scl
                p[poss_idx] = p[poss_idx] * mod
            else:
                if isinstance(fun, str) and fun in dir(st):
                    fun = getattr(st, fun)
                df = fun(**funargs)
                p *= df.pdf(x=cell.z[poss_idx].mean(axis=-1))

            # normalize
            p /= p.sum()
            return p

    def _get_multapsedist(self, i):
        """
        Parameters
        ----------
        i: int
            presynaptic population index

        Returns
        -------
        scipy.stats._distn_infrastructure.rv_sample
            instance of class ``scipy.stats.rv_discrete``
        """
        if hasattr(self.multapseFunction, 'pdf'):
            # assume we're dealing with a scipy.stats.rv_continuous
            # like method. Then evaluate pdf at positive integer
            # values and feed as custom scipy.stats.rv_discrete
            # distribution
            d = self.multapseFunction(**self.multapseParameters[i])
            # number of multapses must be on interval [1, 100]
            # (cap at 100 as more sounds completely unreasonable)
            xk = np.arange(1, 100)
            pk = d.pdf(xk)
            pk /= pk.sum()
            multapsedist = st.rv_discrete(values=(xk, pk))
            # this aint pretty:
            mssg = (
                'multapseFunction: '
                + self.multapseFunction(**self.multapseParameters[i]
                                        ).__str__()
                + f'w. multapseargs: {self.multapseParameters[i]} resulted '
                + f'in {multapsedist.mean()} synapses'
            )
            assert multapsedist.mean() >= 0, mssg
        elif hasattr(self.multapseFunction, 'pmf'):
            # assume we're dealing with a scipy.stats.rv_discrete
            # like method that can be used to generate random
            # variates directly
            multapsedist = self.multapseFunction(**self.multapseParameters[i])
            mssg = (
                'multapseFunction: '
                + self.multapseFunction(**self.multapseParameters[i]
                                        ).__str__()
                + f'w. multapseargs: {self.multapseParameters[i]} resulted '
                + f'in {multapsedist.mean()} synapses'
            )
            assert multapsedist.mean() >= 0, mssg
        else:
            raise NotImplementedError(
                'multapseFunction must be like scipy.stats.rv_discrete '
                + ' or scipy.stats.rv_continuous')

        return multapsedist

    def get_kernel(self, probes, Vrest=-65, dt=2**-4,
                   X='E', t_X=200, tau=50,
                   g_eff=True):
        '''Compute linear spike-to-signal filter kernel mapping presynaptic
        population firing rates/spike trains to signal measurement, e.g., LFP.

        Parameters
        ----------
        probes: list of objects
            list of ``LFPykit.models`` like instances
            (should be instantiated with cell=None).
        Vrest: float
            Mean/Expectation value of postsynaptic membrane voltage used
            for linearization of synapse conductances
        dt: float
            temporal resolution (ms)
        X: str
            presynaptic population for kernel, must be element in
            ``<KernelApprox instance>.X``
        t_X: float
            time of presynaptic event (ms)
        tau: float
            half-duration of filter kernel -- full duration is (2 * tau + dt)
        g_eff: bool
            if True (default), account for contributions by synaptic
            conductances to the effective membrane time constant from
            presynaptic populations X and extrinsic connections.

        Returns
        -------
        H_YX: dict of ndarray
            shape (n_channels, 2 * tau // dt + 1) linear response kernel
        '''
        # get conduction delay transfer function for connections from X to Y
        h_delta = self.get_delay(X, dt, tau)

        # assess index of presynaptic population in X
        (X_i, ) = np.where(np.array(self.X) == X)[0]

        # estimate number of connections as in Potjans&Diesmann2014
        # K_YX = np.log(1. - C_YX) / np.log(1. - 1. / (N_X * N_Y))
        # accurate for small K/(N_X N_Y):
        K_YX = self.C_YX * self.N_X * self.N_Y

        # account for one-sided truncated distribution of synapses per
        # connection with expectation 'loc' and standard deviation 'scale'
        for i in range(len(self.X)):
            multapsedist = self._get_multapsedist(i)

            # total number of connections
            K_YX[i] = K_YX[i] * multapsedist.mean()

        # per neuron indegree
        k_YX_in = K_YX / self.N_Y

        # per neuron outdegree
        k_YX_out = K_YX / self.N_X

        # class NetworkCell parameters:
        cellParameters = deepcopy(self.cellParameters)
        cellParameters.update(dict(
            dt=dt,
            tstop=t_X + tau,
            delete_sections=True
        )
        )

        # Create Cell object representative of whole population
        cell = TemplateCell(**cellParameters)

        # set cell rotation
        cell.set_rotation(**self.rotationParameters)

        # need lists of segment references for each cell in order to shift
        # g_pas per segment
        cell.allseglist = neuron.h.List()
        for sec in cell.allseclist:
            for seg in sec:
                cell.allseglist.append(seg)

        if g_eff:
            # perturb passive leak conductance value due to the missing
            # extrinsic synapses assumed to be homogeneously distributed per
            # surface area.
            # Per-compartment indegree expectation
            rho_ext = self.n_ext * cell.area / cell.area.sum()

            # compute and apply shift of seg.g_pas:
            extPar = self.extSynapseParameters
            if extPar['syntype'] == 'Exp2Syn':
                # compute area under temporal kernel (ms)
                beta = integrate_beta(extPar['tau1'], extPar['tau2'])
                for i in range(cell.totnsegs):
                    g_shift = (
                        rho_ext[i]  # (dimensionless)
                        * abs(extPar['weight'])  # uS
                        / cell.area[i]  # um**2
                        * beta  # ms
                        * self.nu_ext  # s**-1
                    )  # (deci-S/cm**2)
                    cell.allseglist[i].g_pas += (
                        g_shift * 0.1)  # unit: S/cm**2
            else:
                errmsg = '{} not supported'.format(extPar['syntype'])
                raise NotImplementedError(errmsg)

        # iterate over all presynaptic populations in order to offset g_L
        # correctly
        for iii in range(len(self.X)):
            # compute probabilities of connections per compartment which will
            # be used to distribute synaptic input currents.
            # Account for distribution of cells along z-axis:
            # Somatic placement is Gaussian w. z_mean +/- z_std
            # Synapse placement is Gaussian w. u_mean +/- u_std
            # The convolution of two Gaussians is a Gaussian with
            # mean (z_mean + u_mean) and st.dev sqrt(z_std^2 + u_std^2)
            syn_pos = deepcopy(self.synapsePositionArguments[iii])

            for h, funarg in enumerate(syn_pos['funargs']):
                # NOTE: ignoring shifting synapse placements by the mean
                # somatic depth, which may be implemented as:
                # syn_pos['funargs'][h]['loc'] = \
                #     funarg['loc'] + self.populationParameters['loc']
                syn_pos['funargs'][h]['scale'] = \
                    np.sqrt(funarg['scale']**2 +
                            self.populationParameters['scale']**2)

            # per-compartment connection probability
            p = self.get_rand_idx_area_and_distribution_prob(
                cell,
                **syn_pos)

            # synapses per compartment for each connection.
            # Distinguish between input densities and outdegrees as
            # indegrees affects g_pas while outdegrees scale synaptic weights.
            rho_YX_in = p * k_YX_in[iii]
            rho_YX_out = p * k_YX_out[iii]

            # perturb passive leak conductance value due to the linearized
            # synapses
            if g_eff:
                d = self.synapseParameters[iii].copy()
                if d['syntype'] == 'Exp2Syn':
                    # compute area under temporal kernel (ms)
                    beta = integrate_beta(d['tau1'], d['tau2'])
                    for i in range(cell.totnsegs):
                        g_shift = (
                            rho_YX_in[i]  # (dimensionless)
                            * abs(d['weight'])  # uS
                            / cell.area[i]  # um**2
                            * beta  # ms
                            * self.nu_X[self.X[iii]]  # s**-1
                        )  # (deci-S/cm**2)
                        cell.allseglist[i].g_pas += (
                            g_shift * 0.1)  # unit: S/cm**2
                else:
                    errmsg = '{} not supported'.format(extPar['syntype'])
                    raise NotImplementedError(errmsg)

            if iii == X_i:
                # modify synapse parameters to account for current-based
                # synapses linearized around Vrest
                d = self.synapseParameters[iii].copy()
                d['weight'] = - d['weight'] * (Vrest - d['e'])
                del d['e']  # no longer needed
                d['syntype'] = d['syntype'] + 'I'

                # create synapses activated by spike time of presynaptic
                # population X
                # setting weight scaled by synapses per compartment
                for idx in range(cell.totnsegs):
                    di = d.copy()
                    di['weight'] = di['weight'] * rho_YX_out[idx]
                    syn = Synapse(cell, idx=idx, **di)
                    syn.set_spike_times(np.array([t_X]))

        # simulate and record transmembrane currents
        cell.simulate(rec_imem=True)

        # compute and extract kernels
        inds = (cell.tvec >= (t_X - tau)) & (cell.tvec <= (t_X + tau))
        H_YX = dict()
        for probe in probes:
            probe.cell = cell
            if probe.__class__.__name__ in ['PointSourcePotential',
                                            'LineSourcePotential']:
                warn('results are non-deterministic for ' +
                     f'probe {probe.__class__.__name__}.' +
                     'Support may be deprecated in the future.')
                # draw offsets from the distribution of cells in space
                # with homogeneous density
                offsets = self.draw_rand_pos(
                    SIZE=2000,
                    radius=self.populationParameters['radius'],
                    loc=self.populationParameters['loc'],
                    scale=self.populationParameters['scale'])

                M = None
                x_0, y_0, z_0 = probe.x, probe.y, probe.z  # save
                for d in offsets:
                    probe.x = probe.x + d[0]
                    probe.y = probe.y + d[1]
                    probe.z = probe.z + d[2]

                    if M is None:
                        M = probe.get_transformation_matrix()
                    else:
                        M = np.dstack((M, probe.get_transformation_matrix()))

                    probe.x, probe.y, probe.z = x_0, y_0, z_0  # restore

                # spatially averaged transformation matrix
                M = M.mean(axis=-1)
            else:
                M = probe.get_transformation_matrix()

            data = M @ cell.imem

            # apply delay distribution function
            for h, d in enumerate(data):
                data[h, ] = np.convolve(d, h_delta, 'same')

            # subtract kernel offsets at tau == 0:
            data = data - data[:, cell.tvec == t_X]

            # force negative time lags to be zero
            data[:, cell.tvec < t_X] = 0

            # extract kernel
            H_YX[probe.__class__.__name__] = data[:, inds]

            # unset cell
            probe.cell = None

        return H_YX


class GaussCylinderPotential(lfpykit.LinearModel):
    """Compute electric potential of electric sources that are treated as
    inhomogeneous current source density cylinders that are Gaussian along the
    vertical z-axis and constant within a fixed radius in the radial directions
    (xy-plane).

    Parameters
    ----------
    cell: object
        ``CellGeometry`` object or similar
    z: ndarray
        contact point locations
    sigma: float
        conductivity
    R: float
        disk radius
    sigma_z: float > 0
        standard deviation of spatial filter
    """

    def __init__(self, cell, z, sigma=0.3, R=100,
                 sigma_z=50.):
        super().__init__(cell=cell)

        # check input
        assert isinstance(z, np.ndarray), \
            'z must be of type numpy.ndarray'
        assert z.ndim == 1, \
            'z must be of shape (n_coords, )'
        assert isinstance(sigma, float) and sigma > 0, \
            'sigma must be a float number greater than zero'

        # set attributes
        self.z = z
        self.sigma = sigma
        self.R = R
        self.sigma_z = sigma_z

    def _f(self, z_e, z_i):
        return 1 / (2 * self.sigma) * (
            np.sqrt((z_e - z_i)**2 + self.R**2) - abs(z_e - z_i))

    def _g(self, z):
        return np.exp(-(z / self.sigma_z)**2 / 2) / (
            self.sigma_z * np.sqrt(2 * np.pi))

    def _func(self, z, z_e=0, z_i=0):
        return self._f(z_e - z, z_i) * self._g(z)

    def get_transformation_matrix(self):
        '''
        Get linear response matrix

        Returns
        -------
        response_matrix: ndarray
            shape (n_coords, n_seg) ndarray

        Raises
        ------
        AttributeError
            if ``cell is None``
        '''
        if self.cell is None:
            raise AttributeError(
                '{}.cell is None'.format(self.__class__.__name__))
        M = np.empty((self.z.size, self.cell.totnsegs))
        for j, z_e in enumerate(self.z):
            for i in range(self.cell.totnsegs):
                M[j, i], _ = si.quad(self._func, -np.inf, np.inf,
                                     args=(z_e, self.cell.z[i].mean()),
                                     limit=1000)
        M /= (np.pi * self.R**2)  # Distribute current evenly across surface

        return M


class KernelApproxCurrentDipoleMoment(lfpykit.CurrentDipoleMoment):
    """Modified ``lfpykit.CurrentDipoleMoment`` like class that ignores
    contributions to the current dipole moment in the the x- and y-directions
    due to rotational symmetry around the z-axis.

    Parameters
    ----------
    cell: object
        ``CellGeometry`` object or similar
    """

    def __init__(self, cell):
        super().__init__(cell=cell)

    def get_transformation_matrix(self):
        '''
        Get linear response matrix

        Returns
        -------
        response_matrix: ndarray
            shape (3, n_seg) ndarray

        Raises
        ------
        AttributeError
            if ``cell is None``
        '''
        if self.cell is None:
            raise AttributeError(
                '{}.cell is None'.format(self.__class__.__name__))
        return np.stack([np.zeros(self.cell.totnsegs),
                         np.zeros(self.cell.totnsegs),
                         self.cell.z.mean(axis=-1)])
