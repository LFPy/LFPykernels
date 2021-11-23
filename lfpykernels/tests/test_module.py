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

import unittest
import numpy as np
import lfpykit
import lfpykernels


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

    def test_KernelApproxCurrentDipoleMoment(self):
        cell = get_cell(n_seg=3)
        cdm = lfpykernels.KernelApproxCurrentDipoleMoment(cell)
        M = cdm.get_transformation_matrix()

        imem = np.array([[-1., 1.],
                         [0., 0.],
                         [1., -1.]])

        P = M @ imem

        P_gt = np.array([[0., 0.], [0., 0.], [2., -2.]])

        self.assertTrue(np.all(P_gt == P))

    def test_GaussCylinderPotential(self):
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
