Package for Markov-Chain Monte-Carlo tomography and Bayesian earthquake
location.  The idea is pretty simple: try lots of different velocity models
and see which ones explain the data better than others.  Now, each point 
in the velocity model then has a posterior PDF.  Finally, generate a 
set of models from the posterior PDFs, relocate an indepenent earthquake
catalog in each of these sampled models, and ideally out comes earthquake
locations with velocity model uncertainty and a tomography model with some
notion of resolution.

Requirements:

(i)   MPI Fortran and C compilers
(ii)  HDF5
(iii) MKL or ESSL for random number generation

LICENSE:
Apache 2.0

AUTHOR:
Ben Baker
