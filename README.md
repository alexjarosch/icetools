# icetools
icetools provides a development environment for numerical ice flow models/simulations. Nonlinear rheology and Stokes approach are implemented. "Ready to use" programs/scripts are available based on different numerical libraries, from educational to HPC.

This project has been migrated from sourceforge.net to github on February 26th, 2019.

## Requirements

icetools requires a functioning installation of the [FEniCS project](https://fenicsproject.org/).
I recommend you use my singularity container with FEniCS 1.6 installed.
This setup is tested and works. Other versions of FEniCS are currently not supported.

### Install singularity 2.6

Installation instructions for your system can be found [here](https://www.sylabs.io/guides/2.6/user-guide/installation.html).
My container currently supports singularity version 2.6.

### Building the container

First you need to clone my repo to a directory of your liking.
```bash
git clone https://github.com/alexjarosch/icetools.git
```
Now you can build the singularity container with
```bash
cd icetools/FEniCS_container
sudo singularity build fenics_icetools.simg build_script
```
This will take a while, go get some coffee.
