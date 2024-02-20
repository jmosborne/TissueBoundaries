# TissueBoundaries

## Code for D.P.J. Germano, A. Zanca, S.T. Johnston, J.A. Flegg & J.M. Osborne  "Free and Interfacial Boundaries in Individual-Based Models of Multicellular Biological systems"  https://doi.org/10.1007/s11538-023-01214-8

Published in the Buletin of Mathematical Biology in 2023. 

This project contains the code necesary for running the simulations presented in "An adaptive numerical method for multi-cellular simulations of organ development and disease"

Before looking at this, you may wish to look at some of the basic user tutorials for Chaste https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials.

## Getting the code and installing dependencies 

Before running these examples you will need to install Chaste's dependencies (https://chaste.github.io/getting-started/) and the source code for the latest release (release 2024.1) https://github.com/Chaste/Chaste.
The easiest way to do this is using an Ubuntu machine (or an Ubuntu virtual machine) as discussed on https://chaste.github.io/docs/installguides/ubuntu-package/. 
Note that Chaste is only fully supported on Linux/Unix systems, so users of Windows or MacOS we reccomend using Docker (https://github.com/Chaste/chaste-docker).

Once the chaste dependencies and installed and the source code is downloaded go to the `projects` folder and use the command 
`git clone https://github.com/jmosborne/TissueBoundaries.git`  
to download the code for this project.

Now the project should be installed, and everything should compile and run correctly. 
You can now run the tests or simulations, or create your own test suites.

## Documentation
There are two folders - `src` and `test`.
 1. The `src` folder contains the classes necessary to run the simulation. These define the additional classes and  methods not in the core chaste code.
 2. The `test` folder contains:
  * `TestInternalVoid.hpp` - Test suite to run the Void Closure example from Figures 4 and 5.
  * `TestCompetingPopulations.hpp ` - Test suite to run the TissueCollision example from Figures 8 and 9.
  * `TestGrowingMonolayer.hpp` - Test suite to run the Tissue Growth example from Figures 6 and 7.
  
## Running tests
You can then run tests and simulations with (note this assumes the file structure used in the Chaste Docker),
```
cd /lib
cmake ./src
make TestInternalVoid

```

**NB**: the paper was developed with the release 2024.1, it will not work with with release 2021.1 or under.

For further information on using Chaste, see the extensive guide material (https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides).
You may also wish to look at some of the basic user tutorials (https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials/). Note these links will be updated to github website soon.
