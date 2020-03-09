# Evolutionary dynamics code
This respository contains code used to simulate evolution of
cotranslational folding.

The main program, `folding_evolution`, is located at
`src/app_folding_evolution/folding_evolution.cpp`. It is code for an
MPI program. Supporting libraries are in `src`, although some are not
used, and some code within `src` are from other projects.

To compile, an installation of HDF5 and MPI are needed. Furthermore,
the software suite latPack, https://github.com/proteins247/latPack,
must be installed.

# Repository branches
The default branch, `master`, has not been updated in some
time. Several variants of the code are located in different branches:

- `new_fitness_eval`: Evolution of lattice protein sequences under a
  fitness function that evaluates cotranslational folding simulations.
- `no_translation`: Same as above, but protein folding simulations are
  based on *in vitro* folding performance.
- `stability_only`: Fitness evaluation in which protein starts in
  folded state.
- `translation_rate`: Simulation mode where mutations are synonymous
  changes that only change the translation rate.

# Copyright
copyright 2020 by Victor Zhao

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
