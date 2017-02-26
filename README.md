# Lib Gamefft 

### Author: Alioune Schurz

This program enables to identify the main structures with the given targeted molecular weight and seed scaffolds.
The core optimization problem solved within this procedure (referred as CSCCP) is solved in this distribution by implementing a Dynamic programming algorithm (DP):

A GPU accelerated version of the DP, named GAME, developed by Alioune Schurz is implemented in this library. 


## Dependencies 

### A. for basic Functionalities:

- numpy >= 1.7.1
- openbabel-python >= 1.3
- jinja2

### B. for GPU support:

- cuda 5
- pycuda
don't forget to configure the path and to mount the device in /dev (c.f. CUDA installation guide)

### C. for cluster support:

make sure you install hadoop-1.0.3

## Installation

Please make sure that the libgamefft library is installed on your system. (Refer to the section A)

### A. installing the libgamefft

This will compile and install the library in /usr/lib/
From the distribution's root run:

cd lib/libgamefft
make
sudo make install


### B. installing the Python distribution

It is very simple, just do:

1. Install libgamefft (cf. section A)
2. From the distribution's root run:  python setup.py install 
3. Test the distribution by running test.py 

## Command line interface

When installing the distribution, a script called csccp-solver-cli is installed. 

### A. Solving CSCCPs

Thanks to it you can solve CSCCP using different methods. Go at GAME/bin and type:

$> csccp-solver-cli -s GAME/examples/normal/example-data -m 168.195105 -v idp -l cc -r 3 -c 0 --dec 5

The program of csccp-solver-cli (GAME) will use the seed scaffolds included in the sepcified foldder after "-s" (the file in this case is s0000000001) 
to identified the structures having a targeted weight 168.195105 (-m 168.195105) in configuration 0 (-c 0) and ourput the top 3 (-r 3)  most probable
molecules.  The csccp-solver-cli used the DP algorithm (-v idp), implemented in C++ (-v cc).

PS. In general, r is set to 3~10, dec is set to 1~5
    In the file of *.cIdx, assume that the number in line 3 is x, we can set the index of configuration (-c) from 0 to x-1.

You should get this output in this case:

------------------csccp info----------------
scaffold: s0000000001
n: 2
R: 3
mass_peak (w0): 168.195105
mass_peak_min:168.111007448
mass_peak_max: 168.279202553
scaffold_weight_rel2config:106.12194
scaffold_probability_rel2config: 0.026296875
min_possible_weight: 152.120115
max_possible_weight: 498.210115
wmin: 61.9890674475
wmax: 62.1572625525
configuration: [7, 3]
number of sidechains (K): [5, 6]
number of compounds: 30
W,P,sidechain_smiles won't be displayed...
------------------end------------------------
Starting Iterative Dynamic Programming
Finished with 1 iterations: RR=30<=30. Reason: len(filtered_results)=4 < RR
	Probability	Weight		Smile
____________________________________________________________________________________
1.	1.972266e-04	1.681951e+02	C1C/C(=C\CC)C(=O)CC1CO
2.	1.150488e-04	1.681951e+02	C1C/C(=C\CCO)C(=O)CC1C
3.	6.574219e-05	1.681951e+02	C1C/C(=C\O)C(=O)CC1C(C)C


### B. Query information on scaffolds

To query information on a scaffold file you can type:

$> csccp-solver-cli -s data/full_set/normal/s0000000001 -i

You should get this output:

name: s0000000001
popularity: 20
scaffold_weight_plus_hydrogens: 110.1537
min_possible_weight: 152.120115
max_possible_weight: 498.210115
number of configurations: 8
max number of compounds: 0:30 1:36 2:6 3:90 4:60 5:180 6:5 7:20


Other informations
==================
The identified(validated) main structures in our four tesing natural products were provided in the GAME/examples/Datasets/structures

All seed scaffolds (*.cIdx files) in our collected database that can be used as input of the "csccp-solver-cli" program were provided in the GAME/examples/core_index
The index(filename) of seed scaffolds of the four datasets were listed in the GAME/examples/Datasets/scaffolds. (the parameter of -c in csccp-solver-cli)
The targeted molecular weights of the four datasets were also provided in GAME/examples/MW/. (the parameter of -m in csccp-solver-cli)