from dolfin import * 

n = 50 
mesh = UnitSquareMesh(n,n)
V = FunctionSpace(mesh, 'CG',1)
w = 0.3
f3 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= x[0]*x[0] && x[1] >= 0.5 )? 0.0:1.0', w=w)	
k = interpolate(f3, V)
plot(k, interactive = True)
