from dolfin import *

n = 5 
mesh = UnitSquareMesh(n,n)
V = FunctionSpace(mesh, 'Lagrange', 1)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0) # this means that the heat is distributed evenly
#f = Expression('x[0]*x[1]')
g = Constant(50.0) # what does this mean 
u_ = Function(V)  # a function in a finite element function space is a function that is a linear combination of the basisfunctions which make up or in math terms make up that space. 


# we need some boundary perhaps
def boundary(x):
	return x[0]< DOLFIN_EPS or x[0] > 1- DOLFIN_EPS 


u0 = Constant(0.0)
bc = DirichletBC(V,u0,boundary)



a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds
solve(a == L, u_, bc)
plot(u_, interactive = True, title = 'my solution')
