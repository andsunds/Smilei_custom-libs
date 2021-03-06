### NIST_2018_Cu_ionization_energy_eV -- the ionization energies of
### copper in eV, along with their uncertainties.
###
### Data from the NIST Standard Reference Databas v.5.6.1
### https://dx.doi.org/10.18434/T4W30F

######################################################################
# Copyright 2019-2020 ANDRÉAS SUNDSTRÖM
#
# This file is part of `NIST_ionization_energies`.
#
# `NIST_ionization_energies` is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# (version 3) as published by the Free Software Foundation.
#
# `NIST_ionization_energies` is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with `NIST_ionization_energies`. If not, see
# <http://www.gnu.org/licenses/>.
######################################################################

import numpy as np


Cu={'levels': np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]),
    'energies': np.array([7.72638, 20.29239, 36.841, 57.38, 79.8, 103.0, 139.0, 166.0, 198.0, 232.2, 265.33, 367.0, 401.0, 436.0, 483.1, 518.7, 552.8, 632.5, 670.608, 1690.5, 1800, 1918, 2044, 2179.4, 2307.3, 2479.1, 2586.954, 11062.4312, 11567.613]),
    'uncertainty':np.array([4e-06, 6e-05, 0.012, 0.05, 0.7, 1.0, 1.2, 2.1, 2.2, 0.5, 0.25, 0.9, 0.3, 0.6, 0.9, 1.2, 2.1, 0.6, 0.02, 0.9, 3, 4, 6, 1.5, 1.0, 2.2, 0.007, 0.0016, 0.004])
}
