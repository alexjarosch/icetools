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

### Building the FEniCS_container

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
When the container is done, the file `fenics_icetools.simg` will be created, which is the FEniCS_container.

### Using the existing FEniCS_container

If you not like to build the container or face troubles doing so, you can download the latest icetools release, which contains a pre-built container based on the current `build_script`.
This container can also be used with the current development version of icetools.

## Running the examples

### 2D Case

To run the 2D case, first move to the `icetools` directory (maybe you need a `cd ..`).
Now you can start the 2D case right within the singularity container:
```shell
singularity exec FEniCS_container/fenics_icetools.simg python icetools_2d_demo.py
```
which will run `python icetools_2d_demo.py` inside the container and place output in your local directory. Consult the singularity documentation to learn more on working with this type of container environment.

After the run is completed, you can visualize with [paraview](https://www.paraview.org/), as icetools is producing VTK files as output.

The result for the velocity (stored in `velocity_2D.pvd`) should look like this:
![2D Results](https://github.com/alexjarosch/icetools/raw/master/figs/2d_result.jpeg "Velocity in 2D example")

### 3D Case

The 3D case is executed similarly to the 2D case. As this case has a bigger mesh, you have the choice of running it as a sequential code or in parallel. Both of the below runs will take quite some time, depending on how many mesh nodes you use.

Using just a single core you can run:
```shell
singularity exec FEniCS_container/fenics_icetools.simg python icetools_3d_demo.py
```
However if you would like to use 6 cores you can utilize mpirun:
```shell
singularity exec FEniCS_container/fenics_icetools.simg mpirun -np 6 python icetools_3d_demo.py
```
After the run is completed, you can agian visualize with [paraview](https://www.paraview.org/) and the 3D velocity should look like this:
![3D Results](https://github.com/alexjarosch/icetools/raw/master/figs/3d_result.jpeg "Velocity in 3D example")
