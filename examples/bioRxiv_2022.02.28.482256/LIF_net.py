#!/usr/bin/env python
# coding: utf-8
import numpy as np
import scipy.signal as ss
import nest


def get_spike_rate(times, TRANSIENT=2000, tstop=10000, dt=2**-4):
    bins = (np.arange(TRANSIENT / dt, tstop / dt + 1)
            * dt - dt / 2)
    hist, _ = np.histogram(times, bins=bins)
    return bins, hist.astype(float) / dt


def get_mean_spike_rate(times, TRANSIENT=2000, tstop=10000):
    times = times[times >= TRANSIENT]
    return times.size / (tstop - TRANSIENT) * 1000


def get_psd(nu, Fs=16000, NFFT=2048, noverlap=1536,
            detrend='constant', cutoff=200):
    '''Truncated power spectrum'''
    freqs, Pxx = ss.welch(nu, fs=Fs, nperseg=NFFT,
                          noverlap=noverlap, detrend=detrend)
    return freqs[freqs <= cutoff], Pxx[freqs <= cutoff]


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
                 local_num_threads=1,
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
        self.local_num_threads = local_num_threads

        self._create()
        self._connect()

    def _create(self):
        """Create network nodes and connections"""

        nest.ResetKernel()
        nest.SetKernelStatus(
            dict(
                local_num_threads=self.local_num_threads,
                resolution=self.dt,
                tics_per_ms=1000 /
                self.dt))

        if self.verbose:
            print('creating...')
        # neurons
        self.neurons = {}
        for (X, N, C_m, tau_m, E_L, (tau_syn_ex, tau_syn_in)
             ) in zip(self.X, self.N_X, self.C_m_X,
                      self.tau_m_X, self.E_L_X, self.tau_syn_YX):
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

    def create_fir_filters(self, H_YX, tau=100):
        """Create sets of FIR filters for each signal probe

        Parameters
        ----------
        kernels: nested dict
        """
        keys = list(H_YX[f'{self.X[0]}:{self.X[0]}'].keys())
        n_filters_per_kernel = []
        for key in keys:
            n_filters_per_kernel += [
                H_YX[f'{self.X[0]}:{self.X[0]}'][key].shape[0]]

        # create container
        if not hasattr(self, 'fir_filters'):
            self.fir_filters = dict()

        # create sets of FIR_filter nodes for each type of measurement and
        # presynaptic population (ignoring Poisson generators)
        for k in keys:
            self.fir_filters[k] = dict()
            # presynaptic
            for i, X in enumerate(self.X):
                # self.fir_filters[k][X] = []
                # postsynaptic
                for j, Y in enumerate(self.X):
                    # Discarding the elements for time lags tau <= 0 for the
                    # FIR filter
                    tau_inds = np.arange(
                        0, int(2 * tau / self.dt) + 1) > int(tau / self.dt)

                    H = H_YX[f'{Y}:{X}'][k][:, tau_inds]

                    # create one fir_filter node for each row
                    fir_filter_params = []
                    for h in H:
                        fir_filter_params += [dict(
                            N=h.size,  # filter order
                            h=h,  # filter coefficients
                        )]
                    self.fir_filters[k][f'{Y}:{X}'] = nest.Create(
                        'fir_filter_nestml', H.shape[0], fir_filter_params)

        # create recording multimeters, one per probe, sample every dt
        self.multimeters = dict()
        for k in keys:
            self.multimeters[k] = dict()
            for X in self.X:
                for Y in self.X:
                    self.multimeters[k][f'{Y}:{X}'] = nest.Create(
                        'multimeter', 1, {'interval': self.dt, 'label': f'{k}_{Y}:{X}'})
                    self.multimeters[k][f'{Y}:{X}'].set({"record_from": ["y"]})

        # connect FIR filter nodes to the corresponding presynaptic populations
        # at minimum delay
        for k in keys:
            for i, X in enumerate(self.X):
                for Y in self.X:
                    nest.Connect(self.neurons[f'{X}'],
                                 self.fir_filters[k][f'{Y}:{X}'],
                                 syn_spec=dict(delay=self.dt))

        # connect multimeters
        for k in keys:
            for i, X in enumerate(self.X):
                for Y in self.X:
                    nest.Connect(self.multimeters[k][f'{Y}:{X}'],
                                 self.fir_filters[k][f'{Y}:{X}'])

    def simulate(self, tstop=6000):
        """Instantiate and run simulation"""
        if self.verbose:
            print('simulating...')
        nest.Simulate(tstop)
        if self.verbose:
            print('done!')
