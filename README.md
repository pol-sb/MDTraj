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

<p align="center"> The main of this project is to create a simple parelling Molecular Dynamics simulation code. The present code is the final project of the Advanced computation tools for the Atomistic and Multiscale Computational Modelling in Physics, Chemistry and Biochemistry at the Universitat de Barcelona.
    <br> 
</p>

## Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Prerequisites](#prerequisites)
- [Installing](#installing)
- [Input parameters](#parameters)
- [Deployment](#deployment)
- [Usage](#usage)
- [Built Using](#built_using)
- [TODO](../TODO.md)
- [Contributing](../CONTRIBUTING.md)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## About <a name = "about"></a>
In this project we aim to develope a simple parelling Molecular Dynamics simulation code...

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
make build
make
```

To run a simulation you must modify the parameter.h file, in the input directory, se the input parameter section.

## Input parameters. <a name = "parameters"></a>

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


## Authors <a name = "authors"></a>

- [@perasperadastra](https://github.com/perasperadastra) -> Initialization and integration
- [@LucasFernandezStolpa](https://github.com/LucasFernandezStolpa) -> Boundary conditions and forces
- [@pol-sb](https://github.com/pol-sb) -> Statistics and visualization of results
- [@Mtunica](https://github.com/Mtunica) -> Main and makefile

See also the list of [contributors](https://github.com/Eines-Informatiques-Avancades/Project-I/contributors) who participated in this project.

## Acknowledgements <a name = "acknowledgement"></a>

- Hat tip to anyone whose code was used
- Inspiration
- References
