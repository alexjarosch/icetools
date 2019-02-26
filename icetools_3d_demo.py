"""A 3D demo of the icetools code"""

__author__ = "Alexander H. Jarosch (research@alexj.at)"
__date__ = "2010-06-09 -- 2011-09-12"
__copyright__ = "Copyright (C) 2011 Alexander H. Jarosch"
__license__  = "GNU GPL Version 3"
__version__ = "0.9"

"""This file is part of icetools.

icetools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

icetools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with icetools.  If not, see <http://www.gnu.org/licenses/>."""


from dolfin import *
import time
import numpy
import matplotlib.pyplot as plt

# Start the timer to benchmark the run.
tstart = time.time()

# Define some parameters and load the mesh. If you are unfamiliar 
# with ice dynamics and do not understand the code below, go read about
# it, e.g. Jarosch A. H. (2008): Icetools: a full Stokes finite element model
# for glaciers. Computers & Geosciences, 34, 1005-1014, 
# doi:10.1016/j.cageo.2007.06.012
# or 
# Cuffey K. M. and Paterson W. S. B. (2010): The Physics of Glaciers, 4th ed. 
# Chapter 8

cs_length = 40.0        # side length of ice column, used in BC definition
cs_height = 120.0       # height of ice column to be used in BC definition           

max_it = 50             # number of maximum viscosity iterations
max_u_err_def = 0.005   # max velocity error between solutions of u in m/year
max_u_err = max_u_err_def/(365*24*3600) # max error in m/sec
visc_limit = 1e15       # viscosity limit

# Define the body force f, i.e. gravity. Here I use a rotated coordinate
# system similar to the Cuffey & Paterson p. 295 Fig. 8.5a, just in 3D 
alpha = 10.0            # inclination of glacier
alphar = alpha*pi/180
g = -9.81               # gravitational constant
rho = 917.0             # ice density
f_x0 = sin(alphar)*g*rho
f_x2 = cos(alphar)*g*rho
# A note here on how Fenics defines coordinates
# x[0] == x
# x[1] == y
# x[2] == z
f = Constant((f_x0, 0, f_x2))

Aglen = 2.4e-24         # Glen flow parameter for temperate ice (Cuffey & Paterson,2010 p. 73)
nglen = 3.0             # Glen's n
nu = Constant(8e13)     # initial viscosity of ice, before viscosity iteration

# Load the mesh from a gmsh generated file
mesh = Mesh("column3D.xml.gz")
# or select Fenics generated mesh. To use this mesh, uncomment the line below.
# Note that nx, ny, nz define the number of mesh points in each dimension.
# nx = 4
# ny = nx
# nz = 3*nx
# mesh = Box(0, 0, 0, qs_length, qs_length, qs_height, nx, ny, nz)

# Define the sub domains for the boundary conditions
def NoslipBoundary(x, on_boundary):
    return x[2] < DOLFIN_EPS and on_boundary

# Define the periodic boundary on the vertical faces in X direction
class PeriodicBoundary_x(SubDomain):

    def inside(self, x, on_boundary):
        return x[0] == 0 and on_boundary

    def map(self, x, y):
        y[0] = x[0] - cs_length
        y[1] = x[1]
        y[2] = x[2]
# Define the periodic boundary on the vertical faces in Y direction
class PeriodicBoundary_y(SubDomain):

    def inside(self, x, on_boundary):
        return x[1] == 0 and on_boundary

    def map(self, x, y):
        y[1] = x[1] - cs_length
        y[0] = x[0]
        y[2] = x[2]   

# Define function spaces
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# Apply a no-slip boundary condition for velocity
noslip = Constant((0,0,0))
bc0 = DirichletBC(W.sub(0), noslip, NoslipBoundary)
# Apply the periodic boundary condition in X
pbc_x = PeriodicBoundary_x()
bc1 = PeriodicBC(W.sub(0), pbc_x)
# Apply the periodic boundary condition in Y
pbc_y = PeriodicBoundary_y()
bc2 = PeriodicBC(W.sub(0), pbc_y)

# Collect boundary conditions
bcs = [bc0, bc1, bc2]

# Define variational problem
(v_i, q_i) = TestFunctions(W)
(u_i, p_i) = TrialFunctions(W)
# Define the Stokes equation
a = (inner(grad(v_i), nu*grad(u_i)) - div(v_i)*p_i + q_i*div(u_i))*dx
L = inner(v_i, f)*dx

# Compute solution
U = Function(W)
solve(a == L, U, bcs)

# Split the mixed solution using a shallow copy
(u, p) = U.split()

# Save the linear solution in VTK format
ufile_pvd = File("velocity_linear.pvd")
ufile_pvd << u
pfile_pvd = File("pressure_linear.pvd")
pfile_pvd << p

# Now let's compute the analytical solution for parallel flow
# in accordance with Cuffey and Paterson (2010) p. 310. eq. 8.34
# I have written here the equation a bit different to implement it easier.
H = cs_height
Aterm1 = (2*Aglen)/(nglen+1)*(rho*g*sin(alphar))**nglen;
ue_x0 = '%g * (pow(%g,%g+1) - pow(%g-x[2],%g+1))' % (Aterm1, H, nglen, H, nglen)
ue_x1 = '0.0'

# Evaluate the equation on the mesh to create a 3D field of the solution
ue_E = Expression((ue_x0,ue_x1,ue_x1))
u_e = interpolate(ue_E, V)

# Write the solution to a file 
uefile_pvd = File("velocity_exact.pvd")
uefile_pvd << u_e

# Calculate the L2 error norm between the numerical and analytical solution
eps_a = errornorm(u, u_e, "L2")
print 'L2 norm error u-uA:', eps_a*365*24*3600, ' in m/year'

# Calculate the strain invariant and viscosity
V_s = FunctionSpace(mesh, "CG", 2)
v_s = TestFunction(V_s)
w_s = TrialFunction(V_s)

# Define the viscosity function, see e.g. Jarosch (2008)
def nu_ice(u):
    epsxx = u[0].dx(0)
    epsyy = u[1].dx(1)
    epszz = u[2].dx(2)
    epsxy = 0.5*(u[0].dx(1) + u[1].dx(0))
    epsxz = 0.5*(u[0].dx(2) + u[2].dx(0))
    epsyz = 0.5*(u[1].dx(2) + u[2].dx(1))
    eps_dot = sqrt(0.5*(epsxx**2 + epsyy**2 + epszz**2 + 2*epsxy**2 + 2*epsxz**2 + 2*epsyz**2))
    
    return 0.5 * pow(Aglen, (-1/nglen)) * pow(eps_dot, ((1-nglen)/nglen))
# Calculate the non-linear viscosity as a variational form  
a_s = inner(w_s, v_s)*dx
L_s = inner(nu_ice(u), v_s)*dx
nu = Function(V_s)
# Solve for viscosity
solve(a_s == L_s, nu)
# Limit the viscosity due to some crossover stress (e.g. Jarosch 2008)
nu_lim = nu.vector().array()
# find all nu values larger than the limit and replace them with the limit
nu_lim_idx = nu_lim > visc_limit
nu_lim[nu_lim_idx] = visc_limit
nu.vector()[:] = nu_lim

# Also create some numpy arrays to store errors
err_array = numpy.array( [] )
ueps_array = numpy.array( [] )

it = 0      # iteration counter
u_eps = 1   # set to one since we do not have two u iterations to calculate it
# Iterate on the viscosity until either the predefined error is reached or the
# maximum iteration number
while u_eps > max_u_err and it < max_it:
    # Basically this is a classical fixed point iteration.
    # copy the old u to u_k
    u_k = u
    
    # Solve the Stokes problel
    # Redefine the Stokes problem to include the new type of viscosity
    a = (inner(grad(v_i), nu*grad(u_i)) - div(v_i)*p_i + q_i*div(u_i))*dx
    U = Function(W)
    solve(a == L, U, bcs)
    (u, p) = U.split()

    # Solve the Viscosity problem
    # Redefine also the viscosity problem to use new velocity
    L_s = inner(nu_ice(u), v_s)*dx
    nu = Function(V_s)
    solve(a_s == L_s, nu)
    # Safe the new viscosity into an array and limit it once again
    nu_lim = nu.vector().array()
    nu_lim_idx = nu_lim > visc_limit
    nu_lim[nu_lim_idx] = visc_limit
    nu.vector()[:] = nu_lim
    
    # Calculate the L2 error norms
    u_eps = errornorm(u, u_k, "L2")
    eps_a = errornorm(u, u_e, "L2")
    
    # Plug them into an array
    err_array = numpy.append(err_array, eps_a)
    ueps_array = numpy.append(ueps_array, u_eps)
    
    print '########################################'
    print '### L2 u - uA:', eps_a*365*24*3600, ' in m/year'
    print '### L2 diff u:', u_eps*365*24*3600, ' in m/year'
    print '########################################'
    it = it + 1
    
    
t_needed = (time.time() - tstart)

print '\n'
print '====================== R E S U L T ============================================='
if it == max_it:
    print 'WARNING! Max number of iterations taken, error between u and uA still:',  eps_a*365*24*3600, ' in m/year'
    print 'diff U error is: ', u_eps*365*24*3600, ' and should be: ', max_u_err*365*24*3600, ' in mm/year'
    print 'probably your grid size is too large in some regions to get this kind of accuracy'
    print '\n'
    
print "it took: " + str(t_needed) + " seconds to solve the problem ;)"       
print 'All done ;)'
print 'Iterations taken: ', it, ' diff U is: ', u_eps*365*24*3600, ' mm/year'
print '================================================================================'

# Save the final solution
uendfile_pvd = File("velocity_end.pvd")
uendfile_pvd << u

# Plot the errors
plt.figure(1)
plt.subplot(211)
plt.plot(numpy.log(err_array))

plt.xlabel('Iterations')
plt.ylabel('log L2norm u - uA')
plt.title('Difference between analytical and numerical solution')

plt.subplot(212)
plt.plot(numpy.log(ueps_array))
plt.xlabel('Iterations')
plt.ylabel('log L2norm uk+1 - uk')

plt.show()
