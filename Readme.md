This repository contins all the code to reproduce the figures in the paper.
This repository requires:
- Python version 3.6 or higher
- Mathematica and wolframscript installed
- jupyter notebooks
Files are divided by content
- Section 3 contains a Mathematica notebook
- Section 4 contains:
	- lib: a folder with interal routines for python scripts
	- macroscopic-cavity: a folder containing
		- cavity.wls, which produces the files *.txt in this folder
		- macroscopic_cavity.ipynb, which converts the text files into heatmaps
	-pruning-aOC; percolation on synthetic network:
		-type_I: a folder for type I networks, it contains:
			- perc_typeI.py: which run the percolation problem and create a .pkl file with output.
			- type_1.ipynb which combines the simulation output with the theory from micrscopic-cavity folder
			
		-typeII: a folder for type II networks, it contains:
			-perc_typeII.py: as above
			type_II.ipynb: as above
			powers.wls: theory fot the macroscopic cavity of type II networks
	-single-instance: microscopic cavity dynamics
		-percolation.py: generate the output in the file .pkl
		- non-linear.ipyn compare the final state of the simulated dynamics with macroscopic cavity as computed in folder macroscopic-cavity
	-knockout-cascade:
		-knockout_cascade.py: the output is saved in the folder kncokout-cascade_multiple
		-knckout-cascade.ipynb: which loads from the folder  kncokout-cascade_multiple
		
		


pruning-aOC  This folder for a family of synthetic network  a finite fraction of genes are removed, and the  size of the "AND out component" is evauated.
