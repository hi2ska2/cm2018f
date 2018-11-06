import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
######## Initial Condition Set #########
N_y = 9
N_z = 5
Configuration = np.zeros((N_y,N_z))
'''
	model = Configuration.T 
	If N_y=3, N_z = 4,
	Configuration.T is
		1	2	3
		4	5	6 
		7	8	9 
		10	11	12
'''
_boundary = sys.argv[1]
if(_boundary=='1'): # top-right =1, top-left = 0, bottom = 0
	Configuration[N_y-1,0],Configuration[N_y-2,0] = 1,1
elif(_boundary=='2'): # top-right = 0, top-left = 1, bottom = 0
	Configuration[0,0],Configuration[1,0] = 1,1
elif(_boundary=='3'): # top-right = 0, top-left = 1, bottom = 0
	Configuration[:,N_z-1] = np.ones(N_y)
elif(_boundary=='4'): # top-right = 0, top-left = 1, bottom = 0
	Configuration[0,0],Configuration[1,0] = 1,1
	Configuration[N_y-1,0],Configuration[N_y-2,0] = 1,1
	Configuration[:,N_z-1] = np.ones(N_y)
else:
	print("Please type argument after your terminal commend")
	exit()
	
print(Configuration.T)
_matrix_part = np.zeros((N_y*N_z,N_y,N_z))
i = 0
for z in range(N_z):
	for y in range(N_y):
		if(y==0 or y==N_y-1 or z==0 or z==N_z-1):
			if((z==0 and y==0) or (z==0 and y==1) or (z==0 and y==N_y-1) or (z==0 and y==N_y-2) or z==N_z-1):
				_matrix_part[i,y,z]=1
				print(y,z,"is Dirichlet boundary")
			elif(z==0):
				_matrix_part[i,y,z+1],_matrix_part[i,y-1,z],_matrix_part[i,y+1,z]= 1,0.5,0.5
				_matrix_part[i,y,z]=-2
				print(y,z,"is Neumann boundary")
			elif(y==0):
				_matrix_part[i,y,z+1],_matrix_part[i,y,z-1],_matrix_part[i,y+1,z]= 0.5,0.5,1
				_matrix_part[i,y,z]=-2
				print(y,z,"is Neumann boundary")
			elif(y==N_y-1):
				_matrix_part[i,y,z+1],_matrix_part[i,y,z-1],_matrix_part[i,y-1,z]= 0.5,0.5,1
				_matrix_part[i,y,z]=-2
				print(y,z,"is Neumann boundary")
			else:
				print("error")
				exit()
		else:
			_matrix_part[i,y,z-1],_matrix_part[i,y,z+1],_matrix_part[i,y-1,z],_matrix_part[i,y+1,z]= 1,1,1,1
			_matrix_part[i,y,z]=-4
		i=i+1

Matrix = np.reshape(np.swapaxes(_matrix_part,1,2),(-1,N_y*N_z))
Vector = Configuration.T.flatten()
######################################

## result
phi = np.matmul(np.linalg.inv(Matrix),Vector)
print(np.reshape(phi,(N_z,N_y)))

## Plot
fig = plt.figure()
fontsize=12

ax = fig.gca(projection='3d')
Y_grid = np.linspace(1, N_y, N_y)
Z_grid = np.linspace(1, N_z, N_z)
Y_grid, Z_grid = np.meshgrid(Y_grid, Z_grid)

Phi_plot = np.reshape(phi,(N_z,N_y))
sulf=ax.plot_surface(Z_grid, Y_grid, Phi_plot, color=cm.get_cmap('plasma_r', 50), linewidth=2, antialiased=False, alpha=0.7)
ax.set_xlabel("Z",fontsize=fontsize)
ax.set_ylabel("Y",fontsize=fontsize)
ax.set_zlabel("$\mathrm{\Phi}$",fontsize=fontsize)
ax.set_zlim(0.0, 1.0)
plt.colorbar( sulf, extend='both', shrink=0.7,boundaries=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])

plt.show()