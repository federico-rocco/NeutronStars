# NeutronStars

A python package for the study of neutron stars.

## Description

The package can solve a system of differential equations in order to determine the Mass and Radius of a neutron star with a given central pressure. Both the non-relativistic and full 
relativistic case are available using respectively Newton and Tolman-Oppenheimer-Volkoff structures equations.
With the implementation of different pressure-energy density equations of state different kinds of stars can be simulated:
1. a polytropic function in the form <img src="https://latex.codecogs.com/svg.image?P=K\cdot\varepsilon^{\gamma}" title="P=K\cdot\varepsilon^{\gamma}" />, managed by the class `Polytropic`. This eos is used for the non-relativistic case of a pure neutron star.
2. implicit form of the equations <img src="https://latex.codecogs.com/svg.image?P=P(x)" title="P=P(x)" /> and <img src="https://latex.codecogs.com/svg.image?\varepsilon=\varepsilon(x)" title="\varepsilon=\varepsilon(x)" />. This eos is used for a pure neutron star in the full relativistic case.
3. piecewise polytropic eos, where in each range a polytropic eos is implemented. This eos is used for non-pure stars like Read neutron stars.

## Prerequisites

1. astropy
2. matplotlib
3. numpy
4. scipy
5. tqdm

## Installation

You might want to create a virtual environment (`venv`) where install the package:

1. `python3 -m venv venv` (Linux) or `python -m venv venv` (Win)

If so, activate it:

2. `. venv/bin/activate` (Linux) or `venv\Scripts\activate` (Win)

If you want to intall the package in your base environment, simply ignore the first 2 steps.

To download and install the package you can simply use

3. `pip install git+https://github.com/federico-rocco/NeutronStar`

You can also download the zip file with the code and unzip it in a chosen folder.

Then:

4. `cd NeutronStars-main`

The installation steps are:

5. `python -m pip install -r requirements.txt`

to install the prerequisites, then

6. `python setup.py install`

## Usage

The folder 'examples' contains examples on how to use the package. They show how to use the different classes:
1. `NonRelativisticNS.py` uses the eos class `Polytropic`
2. `RelativisticNS.py` uses `Implicit`
3. `ReadNS.py` uses `Piecewise`.

All these scripts produce:
- the profile m(r) and p(r) of the studied star with a particular central value (either pressure or density)
- the mass-vs-radius profile of the star with varying central value
- the mass-vs-central-value and the radius-vs-central-value profile of the star with varying central value

Except that for ReadNS, where the error is too big to be of some insterest, the solution of both Newton and Tolman-Oppenheimer-Volkoff equations is showed. In this way, the impact of the relativistic corrections is shown.

To run tests:

1. `cd nst/examples`

and then run the script of interest, like

2. `python NonRelativisticNS.py`

The produced plots will be available in the folder `NSOutput`.

## Bibliography

"Neutron stars for undergraduates", R. R. Silbar, S. Reddy, 2004.

"Compact Stars for Undergraduates", I. Sagert et al, 2004.

"Constraints on a phenomenologically parametrized neutron-star equation of state", J. Read et al, 2009.
