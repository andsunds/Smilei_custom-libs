### NIST_2018_Li_ionization_energy_eV -- the ionization energies of
### lithium in eV, along with their uncertainties.
###
### Data from the NIST Standard Reference Databas v.5.7
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

Li={'levels': np.array([1, 2, 3]),
    'energies': np.array([5.39171495, 75.6400964, 122.4543581]),
    'uncertainty': np.array([4e-8, 1.3e-6, 8e-7]),
}
