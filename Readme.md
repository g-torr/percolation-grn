# Percolation of gene regulatory network,

This repository contins all the code to reproduce the figures in the paper.
This repository requires:

- Python version 3.6 or higher
- Mathematica and wolframscript installed
- jupyter notebooks

Each folder performs a task. Inside each folder, use the following order to execute files:

1. *.wls
2. *.py
3. *ipynb

Also Section_4/macroscopic-cavity should be run first.

---
## Content of the repository
- Section 3 contains a Mathematica notebook
- Section 4 contains:
	- lib: a folder with interal routines for python scripts
	- macroscopic-cavity: this folder compute the macroscopic cavity theory for type 1 networks.
	- pruning-aOC; percolation  on synthetic network:	
	- single-instance: microscopic cavity dynamics
	- knockout-cascade: single gene removal


## Extensions for bipartite graphs in networkx 
This repository relies heavily on networkx. Some method for direct bipartite graph generation have been extendend. In particular:

- random graph 
- configurational model

Moreover in- and out-component are computed for directed graph.
Go in the `Section-4/lib` to know more
