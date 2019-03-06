# icetools
icetools provides a development environment for numerical ice flow models/simulations. Nonlinear rheology and Stokes approach are implemented. "Ready to use" programs/scripts are available based on different numerical libraries, from educational to HPC.

This project has been migrated from sourceforge.net to github on February 26th, 2019.

## Requirements

icetools requires a functioning installation of the [FEniCS project](https://fenicsproject.org/).
I recommend you use my singularity container with FEniCS 2016.2 installed.
This setup is tested and works. Other versions of FEniCS are currently not supported.

### Install singularity 2.6

Installation instructions for your system can be found [here](https://www.sylabs.io/guides/2.6/user-guide/installation.html).
My container currently supports singularity version 2.6.

### Building the container

First you need to clone my repo to a directory of your liking.
```shell
git clone https://github.com/alexjarosch/icetools.git
```
Now you can build the singularity container with
```shell
cd icetools/FEniCS_container
sudo singularity build fenics_icetools.simg build_script
```
This will take a while, go get some coffee.
When the container is done, the file `fenics_icetools.simg` will be created, which is the FEniCS_container

## Running the examples

### 2D Case

To run the 2D case, first move to the `icetools` directory (maybe you need a `cd ..`).
Now you can start a singularity container shell:
```shell
singularity shell FEniCS_container/fenics_icetools.simg
```
which will bring you inside the container. There you can run the 2D case:
```shell
python icetools_2d_demo.py
```
After the run is completed, you can visualize with [paraview](https://www.paraview.org/), as icetools is producing VTK files as output.
The result for the velocity (stored in `velocity_2D.pvd`) should look like this:
![2D Results](https://github.com/alexjarosch/icetools/raw/master/figs/2d_result.jpeg "Velocity in 2D example")
