# Evolutionary dynamics code
This respository contains code used to simulate evolution of
cotranslational folding.

The main program, `folding_evolution`, is located at
`src/app_folding_evolution/folding_evolution.cpp`. It is code for an
MPI program. Supporting libraries are in `src`, although some are not
used. These libraries were written by others in the Shakhnovich
research group. I added a few functions to `gencode`. Some other
evolutionary simulation code written by others for other research
projects are in `src`.

To compile, an installation of HDF5 and MPI are needed. Furthermore,
the software suite latPack, https://github.com/proteins247/latPack,
must be installed. `install_folding_evolution.sh` is how I compile and
install my `folding_evolution` executable as a module on the Harvard
cluster; it should be modified before use elsewhere.

# Repository branches
Several variants of the code are located in different branches:

- `new_fitness_eval`: Same as master branch. Evolution of lattice
  protein sequences under a fitness function that evaluates
  cotranslational folding simulations.
- `no_translation`: Protein folding simulations for fitness evaluation
  are based on *in vitro* folding performance.
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
