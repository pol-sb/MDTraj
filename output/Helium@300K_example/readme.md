# Helium Simulation

This directory contains the plot of a Helium gas simulation. 

## Initial Parameters
We consider the following parameters for the simulation.

*	Density = 0.005(reduced units)
*	Nc = 6
*	Initial structure = simple cubic
*	Initial velocities = Bimodal distribution at temperature.
*	Temperature = 300 K
*	Time step = 0.001
*	Number of melting steps =10000
*	Number of main simulation steps =100000
*	Collect data every (steps) =100


The density has been calculated by considering a sigma equals to 2.5 Angstroms. The density of Helium at ambient conditions is about 0.18 Kg/m^3 , its mass is approximated 4 g/mol.

  
## Results

We observe with image **ene-allplot.png** and **ene-totplot.png** the conservation of the total energy. Since it is a gas with a extremly low density, the potential interaction is almost imperceptible and the kinetic energy dominates (see figures **ene-kinplot.png**). 

The average temperature oscillates around the 300K (see figure **tempplot.png**). The pressure also oscillates between -0.6 and 0.3, with an average of 0 (figure **pressureplot.png**). 

The radial function distribution has the typical shape of a the one corresponding to a gas. It has a first pick and immediately stabilizes at 1. 
