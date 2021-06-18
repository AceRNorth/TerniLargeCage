# TerniLargeCage

The scripts in this repository were written to study the dynamics of caged populations of Anopheles mosquitoes to accompany the paper "The full suppression of Anopheles gambiae populations in indoor large cages by gene drive technology" by Hammond, Pollegioni, Persampieri et al. (2021).
The script TerniSim.cpp encodes a simulation of the cage population and is accompanied by the header file headSim.h.
RunSimulation is a compiled version of TerniSim.cpp (compiled using the Intel compiler v14.0.2). The simulation takes in input parameters representing demographic/genetic processes, and the set up of the cage (number of days it is followed for etc), and outputs time-series of the frequencies of the different mosquito genotypes.
The script GetSample.wl is a Mathematica package that generates a random parameter vector (from a prior distribution), runs the simulation with these parameters, compares the simulated data to real data from the Cage experiment described in the paper, and finally computes a number of measures of distance from the simulated to the real data. This is used for ABC fitting of the parameters.
The script GenerateSamples.nb is a Mathematica Notebook file that can be used to repeatedly run the GetSample.wl script, in order to conduct ABC parameter estimation.
The Mathematica scripts were written using Mathematica Version 12.2.0.0 on a Linux system.
