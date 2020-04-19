# Percolation of gene regulatory network,

This repository contins all the code to reproduce the figures in the paper.
This repository requires:

- Python version 3.6 or higher
- Mathematica and wolframscript installed
- jupyter notebooks

Files are divided by content. Inside each folder, use the following order to execute files:

1. *.wls
2. *.py
3. *ipynb

Also Section_4/macroscopic-cavity should be run first.

---
## Content of the repository
- Section 3 contains a Mathematica notebook
- Section 4 contains:
	- lib: a folder with interal routines for python scripts
	- macroscopic-cavity: a folder containing
		- cavity.wls, which produces the files *.txt in this folder
		- macroscopic_cavity.ipynb, which converts the text files into heatmaps
	- pruning-aOC; percolation on synthetic network:
		- type_I: a folder for type I networks, it contains:
			- perc_typeI.py: which run the percolation problem and create a .pkl file with output.
			- type_1.ipynb which combines the simulation output with the theory from micrscopic-cavity folder
			
		- typeII: a folder for type II networks, it contains:
			- perc_typeII.py: as above
			- type_II.ipynb: as above
			- powers.wls: theory fot the macroscopic cavity of type II networks
			
	- single-instance: microscopic cavity dynamics
		- percolation.py: generate the output in the file .pkl
		- non-linear.ipyn compare the final state of the simulated dynamics with macroscopic cavity (from folder macroscopic-cavity)
	- knockout-cascade:
		- knockout_cascade.py: the output is saved in the folder kncokout-cascade_multiple
		- knckout-cascade.ipynb: which loads from the folder  kncokout-cascade_multiple
		
		

## Extensions for bipartite graphs in networkx 
This repository relies heavily on networkx. Some method for direct bipartite graph generation have been extendend. In particular:

- random graph 
- configurational model

Moreover in- and out-component are computed for directed graph.
Go in the `Section-4/lib` to know more