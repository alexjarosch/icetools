"""A 2D demo of the icetools code"""

__author__ = "Alexander H. Jarosch (research@alexj.at)"
__date__ = "2010-06-09 -- 2015-05-15"
__copyright__ = "Copyright (C) 2015 Alexander H. Jarosch"
__license__  = "GNU GPL Version 3"
__version__ = "1.1"

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

parameters["form_compiler"]["quadrature_degree"] = 2

# test comment

# define some general constants
g = -9.81               # gravitational constant
rho = 917.0             # fluid density
Aglen = 2.4e-24         # Glen flow parameter for temperate ice (Cuffey & Paterson,2010 p. 73)
nglen = 3.0             # Glen's n
glen_fact = 0.5 * Aglen**(-1.0/nglen)
 
# Define the body force f, i.e. gravity, driving the flow
alpha = 10.0            # inclination of glacier
alphar = alpha*pi/180
f_x0 = sin(alphar)*g*rho
f_x1 = cos(alphar)*g*rho
f = Constant((f_x0, f_x1))

# generate a mesh
m_height = 100
m_width = 200
mesh = RectangleMesh(0,0,m_width,m_height,20,10,'right')

# Define the sub domains for the boundary conditions
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS and on_boundary

# Define the periodic boundary on the vertical faces in X direction
class PeriodicBoundary(SubDomain):

    def inside(self, x, on_boundary):
        return x[0] == 0 and on_boundary

    def map(self, x, y):
        y[0] = x[0] - m_width
        y[1] = x[1]

# defining pbc as the periodic boundary
pbc = PeriodicBoundary()

# Define function spaces
V = VectorFunctionSpace(mesh, "CG", 2, constrained_domain=pbc)
Q = FunctionSpace(mesh, "CG", 1, constrained_domain=pbc)
W = V * Q

w = Function(W)
(u, p) = split(w)
(v, q) = TestFunctions(W)

# Apply a no-slip boundary condition for velocity
noslip = Constant((0.0,0.0))
bc0 = DirichletBC(W.sub(0), noslip, NoslipBoundary())
# Collect boundary conditions
bcs = [bc0]

# comment: to solve the linear problem: u & p needs to be Trialfunctions first before the variational form is defined.

nu = 8e13 # linear viscosity for for first linear guess

epsilon = sym(grad(u))
F = (2*nu*inner(epsilon, grad(v)) - div(u)*q - div(v)*p)*dx - inner(v, f)*dx

dF = derivative(F, w)
pde = NonlinearVariationalProblem(F, w, bcs, dF)
solver = NonlinearVariationalSolver(pde)
solver.parameters["symmetric"] = True
solver.parameters["newton_solver"]["maximum_iterations"] = 5
solver.parameters["newton_solver"]["linear_solver"] = 'superlu_dist'
solver.solve()

# define the nonlinear stokes equation directly:
def visc(u):
    
    eps_dot = sqrt(0.5*inner(sym(grad(u)),sym(grad(u)))) # second invariant of strain
    nu_out = glen_fact * eps_dot**((1.0 - nglen)/nglen)
#     return nu_out
    return Min(nu_out, 2e15)    # would introduce a viscosity limit

nu = visc(u)

epsilon = sym(grad(u))
F = (2*nu*inner(epsilon, grad(v)) - div(u)*q - div(v)*p)*dx - inner(v, f)*dx

dF = derivative(F, w)
pde = NonlinearVariationalProblem(F, w, bcs, dF)
solver = NonlinearVariationalSolver(pde)
solver.parameters["symmetric"] = True
solver.parameters["newton_solver"]["maximum_iterations"] = 80
solver.parameters["newton_solver"]["error_on_nonconvergence"]  = False
solver.parameters["newton_solver"]["relaxation_parameter"] = 0.6
solver.parameters["newton_solver"]["relative_tolerance"] = 1E-4
solver.parameters["newton_solver"]["linear_solver"] = 'superlu_dist'
solver.solve()

(u, p) = w.split(deepcopy=True)

# compare to analytical solution
H = m_height
Aterm1 = (2*Aglen)/(nglen+1)*(rho*g*sin(alphar))**nglen;
ue_x0 = '%g * (pow(%g,%g+1) - pow(%g-x[1],%g+1))' % (Aterm1, H, nglen, H, nglen)
ue_x1 = '0.0'
# Evaluate the equation on the mesh to create a 3D field of the solution
ue_E = Expression((ue_x0,ue_x1))
u_e = interpolate(ue_E, V)

nu_out = project(nu, Q)

# Save the final solution
uendfile_pvd = File("velocity_2D.pvd")
uendfile_pvd << u
u_eendfile_pvd = File("velocity_2D_analytical.pvd")
u_eendfile_pvd << u_e
pendfile_pvd = File("pressure_2D.pvd")
pendfile_pvd << p
nuendfile_pvd = File("nu_2D.pvd")
nuendfile_pvd << nu_out

plot((u[0]-u_e[0])*365*24*3600., title="Velocity error in x")
plot(u[1]*365*24*3600., title="Velocity error in y")
plot(p, title="Pressure")
plot(nu_out, title="Viscosity")
interactive(True)

