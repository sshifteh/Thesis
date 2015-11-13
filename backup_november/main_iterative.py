

from dolfin import * 
from indicator_function_expression import ind_func_expr
from K_update import K_update
from capping import capping
from stokes_solver import stokes_solver 
from WSS import WSS
#from WSS_bdry import WSS_bdry  
#from indicator_funciton_conditions2 import ind_func_expr

from indicator_function_v3 import indicator_function

n = 5
w=0.3
mesh = UnitSquareMesh(n,n, 'left/right')
Pspace = FunctionSpace(mesh, 'CG',1)
Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
Pspace = FunctionSpace(mesh, 'Lagrange', 1)
Wspace = MixedFunctionSpace([Vspace, Pspace]) 
DG1 =FunctionSpace(mesh, 'DG', 1) # ShearStressSpace
DG0 = FunctionSpace(mesh, 'DG', 0)	  # DPspace CHANGE this name to DG0 much more intuitive than discont. pressure space

# iterate 50 times
# return the K_func 
# store it to xdmf file to open it paraview 

# step 1: Take K old 
K_1 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] < 0.7 )? 0.0:1.0', w=w)


# perhaps using round in the capping might solve the problem 


def iterative(K_1):
	# step 2: Send it to Stokes solver 
	U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_1, n=n )

	# step 3: calculate the WSS over it and its bdry 
	wss = WSS(U=U, DG1=DG1)
	#wss_bdry = WSS_bdry(K_Func=K_Func, DPspace=DPspace, WSS=wss)

	# Step 4: This is the diffucult step. we need to create an indicator function now 
	indicator_function(WSS=wss, K_Func=K_Func, DG0=DG0, DG1=DG1)	


	#ind_f = indicator_function(WSS=wss, K_Func=K_Func, DG0=DG0, DG1=DG1 )
	# works!! finally 

	# step 5: add together old and new K functions 
	# ind_f is a function in DG0 space and old K function is inPspace
	# have to solve this first 

	#K_new = K_update(ind_func = ind_f, K_Func = K_Func)
	#plot(K_new, interactive = False, title = 'K new ')
	# doesnt look that great, lets try capping 

	#step 6: capping 
	#K_capped = capping(K_new)
	#plot(K_capped, interactive = False, title = 'K capped')
	# actually looks really good, 

	# step 7: run fluid trough it 
	#stokes_solver(w=w, mesh=mesh, Vspace=Vspace, Pspace=Pspace, Wspace=Wspace, K_array=K_capped, n=n )
	#return K_capped 
	


# why the wss_bdry is getting bigger and bigger ? 

iterative(K_1=K_1)

#for i in range(5):
	#K_updated = iterative(K_1 = K_1)
	#K_1= K_updated # this didnt work 
	


