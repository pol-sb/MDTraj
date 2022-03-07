<p align="center">
  <a href="" rel="noopener">
 <img width=200px height=200px src="./imgs/cs_latt.png" alt="Project logo"></a>
</p>

<h3 align="center">Molecular Dynamics Simulation of a Van der Waals Gas</h3>

<div align="center">

[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![GitHub Issues](https://img.shields.io/github/issues/Eines-Informatiques-Avancades/Project-I.svg)](https://github.com/Eines-Informatiques-Avancades/Project-I/issues)
[![GitHub Pull Requests](https://img.shields.io/github/issues-pr/Eines-Informatiques-Avancades/Project-I.svg)](https://github.com/Eines-Informatiques-Avancades/Project-I/pulls)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)

</div>

---

<p align="center"> The main of this project is to create a simple parallel Molecular Dynamics simulation code. It is the final project of Advanced computation tools course of the Atomistic and Multiscale Computational Modelling in Physics, Chemistry and Biochemistry Master at the Universitat de Barcelona.
    <br> 
</p>

## Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Prerequisites](#prerequisites)
- [Installing](#installing)
- [Input parameters](#parameters)
- [Output files and plots](#output)
- [Deployment](#deployment)
- [Usage](#usage)
- [Built Using](#built_using)
- [TODO](#todo)
- [Contributing](../CONTRIBUTING.md)
- [Authors](#authors)

## About <a name = "about"></a>
In this project we aim to develope a simple parelling Molecular Dynamics simulation code. We implement three possible initial structures (sc, fcc, diamond) and two initial configurations (bimodal distribution or to 0). We implement two different integration algorithms (the velocity verlet with and without thermostat andersen). The potential used is the Lennard-Jones.

## Getting Started <a name = "getting_started"></a>

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See [deployment](#deployment) for notes on how to deploy the project on a live system.

### Prerequisites  <a name = "prerequisites"></a>

The core of this program works in FORTRAN 90, so a FORTRAN compiler must be installed before trying to build the code. We higly reccoment to install [gfortran](https://gcc.gnu.org/wiki/GFortran). It can be installed using the following commands:
```
sudo apt-get update
sudo apt-get install gfortran
```

Python is used for the results plotting and representation. A python version higher or equal than `python 3.6` is needed, and additionally the following libraries are needed:

- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)

To install these libraries, the python package manager `pip` can be used. The required versions are available from the included [requirements.txt](src/requirements.txt) file and can be easily installed using the following command while in the same directory:

```
python -m pip install -r requirements.txt
```
`pip` can normally be installed from your distribution package manager.

### Installing. <a name = "installing"></a>

Donwload the zip file and uncompress in your working directory, you can use:


```
 unzip Project-I-master.zip 
```

Move to main directory and tip:

```
make
```

To run a simulation you must modify the parameter.h file, in the input directory, se the input parameter section. Since it is necessary to recompile if you want to add the changes on the input parameters, we recommend to use the command

```
make
```
or
```
make all
``` 
Nonetheles, it can be done in three steps by the following terminal commands:
```
make compile
make run
make plot
```



## Input parameters. <a name = "parameters"></a>

**If not specified, all the units are in reduced units**

The number of unit cells that are simulated is choosed in the **nc** parameters, 

The density have units of (particles / reduced units of distances), take special care in avoid density greater than 0.6 .

The dimension parameters can't be changed in this version.

Can be generated three type of lattice: Simple cubic (=1), face cubic centered (=2), diamond (=3). Also it will be implemented a read from file subroutine in the next version.

Temperature is in kelvin units.

The initial velocities of the particles can be setted up to 0 or to a bimodal distribution. To choose use the **vel_opt** parameter (= 1, bimodal) (= 0, zero initial velocities).

To control de time step change the **dt** parameter (in ps units) less than 0.01 ps make in the simulation unstable. The number of steps that will be computed is selected in the **ntimes** parameter. 

The output information frequency is controlled by the **everyt** parameter. AVOID HIGH RATES OF PRINTIN!!! Printing is a limiting stage in this software.

The **rc** parameter represents the cut-off of the forces calculation, at higher cut-off better precision but higher times of calculation.

At last, you must choose the parameters for the force-field ([Lennard-Jones](https://es.wikipedia.org/wiki/Potencial_de_Lennard-Jones) type). **σ** is the distance to the zero potential point in hte potential and **ε** is the depth of the potential well.


## Output files and plots <a name = "output"></a>


Containing the initial structure:

  * **init_conf_sc.xyz**: The initial simple cubic structure is stored.
  * **init_conf_fcc.xyz**: The initial face centered cubic structure is stored.
  * **init_conf_diamond.xyz**: The initial diamond structure is stored.

Containing the thermodynamics parameters:

  * **temp.dat**: It contains the temperatures of the temperature for some time-steps.
  * **energy.dat**: It contains the energy of the temperature for some time-steps
  * **pressure.dat**: It contains the pressure of the temperature for some time-steps
  * **rdf.dat**: It contains the data of the radial distribution function.

The plots:
 * **ene-allplot.png**: Evolution of the kinetic, potential and total energy.
 * **ene-kinplot.png**: Evolution of the kinetic energy.
 * **ene-potplot.png**: Evolution of the potential energy.
 * **ene-totplot.png**: Evolution of the tiotal energy.
 * **presseplot.png**: Evolution of the pressure.
 * **rdfplot.png**:  Radial distribution function.
 * **tempplot.png**: Evolution of the temperature. 

## Running the tests <a name = "tests"></a>

Tests will be implemented in the next version
### Break down into end to end tests


## Usage <a name="usage"></a>

Add notes about how to use the system.

## Deployment <a name = "deployment"></a>

Add additional notes about how to deploy this on a live system.

## Built Using <a name = "built_using"></a>

- [Fortran](https://fortran-lang.org/) - Fortran
- [Express](https://www.python.org/) - Python
- Modules, etc...?

## TODO <a name = "todo"></a>
**Parallelization:**
- Initialization:
  - sc: Raul
  - fcc: Raul
  - diamond: Marc
  - bimodal: Marc
- Thermostat:
  - kinetic: Lucas
  - andersen thermostat: Raul
- Integration:
  - euler: Pol
  - velocity verlet: Pol
  - velocity verlet with thermostat: Marc
- Forces:
  - Force + Lj: Lucas
- Statistic:
  - g(r): Pol     

## Authors <a name = "authors"></a>
- [@LucasFernandezStolpa](https://github.com/LucasFernandezStolpa) - Coordinator
- [@perasperadastra](https://github.com/perasperadastra)
- [@pol-sb](https://github.com/pol-sb)
- [@Mtunica](https://github.com/Mtunica)

See also the list of [contributors](https://github.com/Eines-Informatiques-Avancades/Project-I/contributors) who participated in this project.

