# Percolation on the gene regulatory network,
## Giuseppe Torrisi et al J. Stat. Mech. (2020) 083501 https://doi.org/10.1088/1742-5468/aba7b0.

See this [twitter thread](https://twitter.com/pepp_torr/status/1260539098277859328) for a short summary of what is in the paper.

This repository contins all the code to reproduce the figures of the paper.
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
- lib: a folder with interal routines for python scripts
- Section 3 contains a Mathematica notebook
- Section 4 contains:

	- macroscopic-cavity: this folder compute the macroscopic cavity theory for type 1 networks.
	- pruning-aOC; percolation  on synthetic network:	
	- single-instance: microscopic cavity dynamics
	- knockout-cascade: single gene removal
- Section 5 contains:
	strongly connected component evaluation.


## Extensions for bipartite graphs in networkx 
This repository relies heavily on networkx. Some method for direct bipartite graph generation have been extendend. In particular:

- directed biparite random graph generation 
- configurational model for directed bipartite graph

Moreover in- and out-component are computed for directed graph.
Go in the `lib` folder  to know more
