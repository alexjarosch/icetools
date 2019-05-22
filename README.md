# icetools
icetools provides a development environment for numerical ice flow models/simulations. Nonlinear rheology and Stokes approach are implemented. "Ready to use" programs/scripts are available based on different numerical libraries, from educational to HPC.

This project has been migrated from sourceforge.net to github on February 26th, 2019.

### Initial development

In 2008 the first version of icetools was released. This version was based on the
[GetFEM++](http://getfem.org/) library and its Matlab (R) interface.
A corresponding publication was describing the initial tests and implementation:

A. H. Jarosch, “Icetools: A full Stokes finite element model for glaciers”, Computers & Geosciences, Vol. 34, iss. 8, p. 1005-1014, 2008, doi: [10.1016/j.cageo.2007.06.012](http://dx.doi.org/10.1016/j.cageo.2007.06.012).

Later on around 2010-2011 the code was migrated to utilize the FEniCS project.

### Publications/Projects utilizing icetools

Here a brief and incomplete list of publications that have utilized icetools:
* A. Wirbel, A. H. Jarosch, and L. Nicholson, “Modelling debris transport within glaciers by advection in a full-Stokes ice flow model,” The Cryosphere, Vol. 12, Iss. 1, p. 189–204, 2018, doi: [10.5194/tc-12-189-2018](http://dx.doi.org/10.5194/tc-12-189-2018).
* A. H. Jarosch and M. T. Gudmundsson, “A numerical model for meltwater channel evolution in glaciers,” The Cryosphere, Vol. 6, Iss. 2, p. 493–503, 2012, doi: [10.5194/tc-6-493-2012](http://dx.doi.org/10.5194/tc-6-493-2012).
* A. H. Jarosch and M. T. Gudmundsson, “Numerical studies of ice flow over subglacial geothermal heat sources at Grímsvötn, Iceland, using Full Stokes equations,” Journal of Geophysical Research-Earth Surface, vol. 112, p. F02008, 2007, doi: [10.1029/2006JF000540](http://dx.doi.org/10.1029/2006JF000540).
* A. H. Jarosch, “Icetools: A full Stokes finite element model for glaciers”, Computers & Geosciences, Vol. 34, iss. 8, p. 1005-1014, 2008, doi: [10.1016/j.cageo.2007.06.012](http://dx.doi.org/10.1016/j.cageo.2007.06.012).

Software projects utilize icetools:
* [debadvect](https://github.com/awirbel/debadvect).

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

If you don't like to build the container or face troubles doing so, you can download a pre-built container based on the current `build_script` from the icetools 1.3.0 release [here](https://github.com/alexjarosch/icetools/releases).
This container, called `fenics_icetools.simg`, can also be used with the current development version of icetools.

## Running the examples

### 2D Case

To run the 2D case, first move to the `icetools` directory (maybe you need a `cd ..`).
Now you can start the 2D case right within the singularity container:
```shell
singularity exec FEniCS_container/fenics_icetools.simg python3 icetools_2d_demo.py
```
which will run `python icetools_2d_demo.py` inside the container and place output in your local directory. Consult the singularity documentation to learn more on working with this type of container environment.

After the run is completed, you can visualize with [paraview](https://www.paraview.org/), as icetools is producing VTK files as output.

The result for the velocity (stored in `velocity_2D.pvd`) should look like this:
![2D Results](https://github.com/alexjarosch/icetools/raw/master/figs/2d_result.jpeg "Velocity in 2D example")

### 3D Case

The 3D case is executed similarly to the 2D case. As this case has a bigger mesh, you have the choice of running it as a sequential code or in parallel. Both of the below runs will take quite some time, depending on how many mesh nodes you use.

Using just a single core you can run:
```shell
singularity exec FEniCS_container/fenics_icetools.simg python3 icetools_3d_demo.py
```
However if you would like to use 6 cores you can utilize mpirun:
```shell
singularity exec FEniCS_container/fenics_icetools.simg mpirun -np 6 python3 icetools_3d_demo.py
```
After the run is completed, you can agian visualize with [paraview](https://www.paraview.org/) and the 3D velocity should look like this:
![3D Results](https://github.com/alexjarosch/icetools/raw/master/figs/3d_result.jpeg "Velocity in 3D example")
