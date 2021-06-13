## dual-method-tiler

Generate2D_PA.m: Raw mathematica code containing various functions related to dual method

Run_Generate_LI_QC.nb: Script that can generate a tiling. You'll have to copy the contents of the file into a mathematica notebook, and save the notebook in the folder that contains GridMethod2D_PA.m

Inputs
- gammasum: the LI class label
- seed: seed for random number generator, which generates the phases subject to constraint that phases sum up to gammasum
- accuracy (play around with this to give different degrees of approximant... I've provided four values that would give the lowest four approximants)
- kmin and kmax: determine the number of gridlines in a grid (kmax - kmin + 1). The value needs to be sufficiently large to get the tiles for the unit cell. I think -7 to 7 was sufficient for the accuracy = 0.01 approximant.