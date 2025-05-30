from dolfin import *
import matplotlib.pyplot as plt

# Parameters
n = 20
eps = Constant(0.001)

# Domain and mesh
mesh = UnitSquareMesh(n, n)
plot(mesh, title="mesh")
plt.savefig("gr1.png", format="png", dpi=200)
plt.clf()

# Define function space
V = FunctionSpace(mesh, "CG", 2)

# Subdomain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] + x[1] > 1.0 - DOLFIN_EPS and on_boundary
# Define boundary condition
bc = DirichletBC(V, Constant(0.0), DirichletBoundary())

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("x[0]", degree=2)
a = eps*inner(grad(u), grad(v))*dx + u*v*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)
print ("max|u(x)| = ", u.vector().norm('linf') )

# Plot  
plot(u, title="Solution")
plt.savefig("gr2.png", format="png", dpi=200)
plt.clf()

