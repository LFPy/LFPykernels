#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Initialization of LFPykernels

Copyright (C) 2021 Computational Neuroscience Group, NMBU.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

:Classes:
  * KernelApprox:
        Class for computing linear spike-to-signal filter kernels resulting
        from presynaptic spiking activity and resulting postsynaptic currents
  * GaussCylinderPotential:
        Compute electric potential of electric sources that are treated as
        inhomogeneous current source density cylinders that are Gaussian along
        the vertical z-axis and constant within a fixed radius in the radial
        directions (xy-plane).
  * KernelApproxCurrentDipoleMoment:
        Modified ``lfpykit.CurrentDipoleMoment`` class that ignores contributions
        to the current dipole moment in the the x- and y-directions
        due to rotational symmetry around the z-axis.


:Modules:
  * kernel_approx
  * version
"""

from .version import version as __version__
from .kernel_approximator import (KernelApprox,
                                  GaussCylinderPotential,
                                  KernelApproxCurrentDipoleMoment)
