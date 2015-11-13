from dolfin import * 

n = 100
mesh = UnitSquareMesh(n,n)
facets = FacetsFunction("size_t", mesh)  # a vector of size n
facets.array()
print facets.array()

