The detailed structure of the Section 4:
There are the following folders:

  - macroscopic-cavity:   this folder compute the macroscopic cavity theory for type 1 networks,
  - pruning-aOC: percolation  on synthetic network,
  - single-instance: microscopic cavity dynamics,
  - knockout-cascade: single gene removal
  
####  Inner structure of the folders
	- macroscopic-cavity:   this folder compute the macroscopic cavity theory for type 1 networks.
		- cavity.wls, which produces the files *.txt in this folder
		- macroscopic_cavity.ipynb, which converts the text files into heatmaps
		
	- pruning-aOC;   percolation  on synthetic network:
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
	- knockout-cascade: single gene removal
		- knockout_cascade.py: the output is saved in the folder kncokout-cascade_multiple
		- knckout-cascade.ipynb: which loads from the folder  knockout-cascade_multiple
		
		
