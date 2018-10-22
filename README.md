# XXZ - spin chain simulation 
A simulator for different types of Heisenberg XXZ spin chains with spin 1/2-particles written in python. For an easy introduction to simulating a spin chain with the code use the introduction.ipynb notebook. If you want to understand the simulation in more detail you can look through the .py files.

# Introduction
This is the simulation written for my bachelor thesis. General physics description can be found in bsc.pdf (To be inserted!). A technical description of the code will be provided in the future. 
This code is intended to serve as a introduction to spin chain simulations and might be not the fastest code to simulate the dynamics of such a quantum system. However it is easy to extend the code with needed functions or add it to a parallel cluster setup. 

# Code structure
The main code is defined in four python-files ( np_helper.py, permutations.py, qm_helper.py, rungekutta.py), which gives the basic functionality. The simulations are written with Jupyter Notebooks. However all notebooks, analysing the dynamics of a system, share the same general setup, which is introduced in a introduction.ipynb. Other example notebooks will be uploaded after some code cleanup. In the end all notebooks used to produce the results from bsc.pdf will be published in this repository.




