Last modified on 1/19/2017 by Rostam


The following code generates sequences of two cubic lattice proteins
undergoing monoclonal evolution of monomers or dimers. 


Directory name	Contents
./commondata	lattice protein conformations - representative subset of 10000
					of all 103346 structures
				physical values of amino acids such as contact energies from Miyazawa and
					Jernigan 1996 
./src			functions for folding and binding
./RNG			random number generator
./example		contains code to do protein design from a starting sequence in seq.txt 


To run example:
1) cd example/
2) ./compile.sh
3) ./protein-design.exe 1.0 1 1 1 2 0

The path of the output file is dim2/v0.dat 

Check out ./protein-design.exe --help for information on the input arguments
