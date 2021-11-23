#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Copyright (C) 2021 Computational Neuroscience Group, NMBU.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

"""

import os
import sys
import posixpath
import unittest
import numpy as np
import lfpykit
import lfpykernels
import scipy.stats as st
import LFPy
import neuron

# for tests to run load mechanisms
if "win32" in sys.platform:
    pth = os.path.join(lfpykernels.__path__[0], 'tests', 'nrnmech.dll')
    pth = pth.replace(os.sep, posixpath.sep)
    if pth not in neuron.nrn_dll_loaded:
        neuron.h.nrn_load_dll(pth)
        neuron.nrn_dll_loaded.append(pth)
else:
    neuron.load_mechanisms(os.path.join(lfpykernels.__path__[0], 'tests'))


def get_cell(n_seg=4):
    cell = lfpykit.CellGeometry(x=np.array([[0.] * 2] * n_seg),
                                y=np.array([[0.] * 2] * n_seg),
                                z=np.array([[1. * x, 1. * (x + 1)]
                                            for x in range(n_seg)]),
                                d=np.array([1.] * n_seg))
    return cell


class TestSuite(unittest.TestCase):
    """
    test methods and modules
    """

    def test_TestSuite_00(self):
        '''test TestSuite'''
        self.assertTrue(True)

    def test_KernelApproxCurrentDipoleMoment_00(self):
        '''cell with radial orientation'''
        cell = get_cell(n_seg=3)
        cdm = lfpykernels.KernelApproxCurrentDipoleMoment(cell)
        M = cdm.get_transformation_matrix()

        imem = np.array([[-1., 1.],
                         [0., 0.],
                         [1., -1.]])

        P = M @ imem

        P_gt = np.array([[0., 0.], [0., 0.], [2., -2.]])

        self.assertTrue(np.all(P_gt == P))

    def test_KernelApproxCurrentDipoleMoment_01(self):
        '''cell with tangential orientation'''
        cell = get_cell(n_seg=3)
        cell.x = cell.z
        cell.z = np.zeros_like(cell.z)
        cdm = lfpykernels.KernelApproxCurrentDipoleMoment(cell)
        M = cdm.get_transformation_matrix()

        imem = np.array([[-1., 1.],
                         [0., 0.],
                         [1., -1.]])

        P = M @ imem

        P_gt = np.array([[0., 0.], [0., 0.], [0., 0.]])

        self.assertTrue(np.all(P_gt == P))

    def test_GaussCylinderPotential_00(self):
        '''check that forward solution is similar to that of point source
        in limit of small disk radius, small standard deviation and large
        source to measurement site separation'''
        cell = get_cell(n_seg=3)

        z = np.array([-1000., 1000.])
        gauss = lfpykernels.GaussCylinderPotential(
            cell=cell, z=z, sigma=0.3, R=.1, sigma_z=.1
        )
        point = lfpykit.PointSourcePotential(
            cell=cell, x=np.zeros(2), y=np.zeros(2), z=z, sigma=0.3
        )

        np.testing.assert_allclose(point.get_transformation_matrix(),
                                   gauss.get_transformation_matrix())

    def test_KernelApprox_00(self):
        '''test that the basic methods of the KernelApprox class works'''

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

        # parameters
        Vrest = -65
        dt = 2**-4
        tau = 100

        cellParameters = dict(
            templatefile=os.path.join(lfpykernels.__path__[0], 'tests',
                                      'BallAndSticksTemplate.hoc'),
            templatename='BallAndSticksTemplate',
            custom_fun=[set_passive],
            custom_fun_args=[{'Vrest': Vrest}],
            templateargs=None,
            delete_sections=False,
            morphology=os.path.join(lfpykernels.__path__[0], 'tests',
                                    'BallAndSticks_E.hoc')
        )

        populationParameters = dict(radius=100, loc=0, scale=50)

        params = dict(
            X=['E'],
            Y='E',
            N_X=np.array([1024]),
            N_Y=1024,
            C_YX=np.array([0.1]),
            cellParameters=cellParameters,
            populationParameters=populationParameters,
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
        )

        # instantiate object
        kernel = lfpykernels.KernelApprox(
            **params
        )

        # check delay distribution
        delay = kernel.get_delay('E', dt=dt, tau=tau)
        assert delay.shape == (int(tau // dt) * 2 + 1, )
        assert delay.min() == 0
        assert np.all(delay[:int(tau // dt) + 1] == 0)
        np.testing.assert_almost_equal(delay.sum(), 1)

        # check get_rand_idx_area_and_distribution_prob method with
        # a uniform distribution
        cell = LFPy.TemplateCell(**cellParameters)
        p = kernel.get_rand_idx_area_and_distribution_prob(
            cell,
            section='allsec',
            fun=st.uniform,
            funargs=dict(loc=-1E3, scale=2E3))
        np.testing.assert_allclose(p, cell.area / cell.area.sum())

        # check kernel predictions using two different forward models
        z = np.linspace(-500, 500, 11)
        gauss = lfpykernels.GaussCylinderPotential(
            cell=None, z=z, sigma=0.3, R=100, sigma_z=100
        )

        cdm = lfpykernels.KernelApproxCurrentDipoleMoment(cell=None)

        H_EE = kernel.get_kernel(probes=(gauss, cdm), Vrest=Vrest, X='E',
                                 dt=dt, tau=tau, t_X=200, g_eff=True)

        assert (H_EE['GaussCylinderPotential'].shape
                == (z.size, int(2 * tau // dt) + 1))
        assert np.all(H_EE['GaussCylinderPotential'][:, :int(tau // dt) + 1]
                      == 0)
        assert (H_EE['KernelApproxCurrentDipoleMoment'].shape
                == (3, int(2 * tau // dt) + 1))
        assert np.all(
            H_EE['KernelApproxCurrentDipoleMoment'][:, :int(tau // dt) + 1]
            == 0)
