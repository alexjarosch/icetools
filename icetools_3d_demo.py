from dolfin import *

parameters["form_compiler"]["quadrature_degree"] = 2

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
f_x2 = cos(alphar)*g*rho
f = Constant((f_x0, 0, f_x2))

# generate a mesh
Lx = 200.
Ly = 200.
Lz = 100.
Nx = 20
Ny = 20
Nz = 10
mesh = BoxMesh(0., 0., 0., Lx, Ly, Lz, Nx, Ny, Nz)

class PeriodicBoundary(SubDomain):

    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two slave edges
        return bool((near(x[0], 0) or near(x[1], 0)) and 
            (not ((near(x[0], Lx) and near(x[1], 0)) or 
                  (near(x[0], 0) and near(x[1], Ly)))) and on_boundary)

    def map(self, x, y):
        if near(x[0], Lx) and near(x[1], Ly):
            y[0] = x[0] - Lx
            y[1] = x[1] - Ly
            y[2] = x[2]
        elif near(x[0], Lx):
            y[0] = x[0] - Lx
            y[1] = x[1]
            y[2] = x[2]
        elif near(x[1], Ly):
            y[0] = x[0]
            y[1] = x[1] - Ly
            y[2] = x[2]
        else:
            y[0] = -1000
            y[1] = -1000
            y[2] = -1000

# Define the sub domains for the boundary conditions
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] < DOLFIN_EPS and on_boundary

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
noslip = Constant((0.0,0.0,0.0))
bc0 = DirichletBC(W.sub(0), noslip, NoslipBoundary())
# Collect boundary conditions
bcs = [bc0]

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
H = Lz
Aterm1 = (2*Aglen)/(nglen+1)*(rho*g*sin(alphar))**nglen;
ue_x0 = '%g * (pow(%g,%g+1) - pow(%g-x[2],%g+1))' % (Aterm1, H, nglen, H, nglen)
ue_x1 = '0.0'
ue_x2 = '0.0'
# Evaluate the equation on the mesh to create a 3D field of the solution
ue_E = Expression((ue_x0,ue_x1,ue_x2))
u_e = interpolate(ue_E, V)

nu_out = project(nu, Q)

u_err = project(u - u_e, V)

# Save the final solution
uendfile_pvd = File("velocity_3D.pvd")
uendfile_pvd << u
u_eendfile_pvd = File("velocity_3D_analytical.pvd")
u_eendfile_pvd << u_e
u_errendfile_pvd = File("velocity_3D_err.pvd")
u_errendfile_pvd << u_err
pendfile_pvd = File("pressure_3D.pvd")
pendfile_pvd << p
nuendfile_pvd = File("nu_3D.pvd")
nuendfile_pvd << nu_out


