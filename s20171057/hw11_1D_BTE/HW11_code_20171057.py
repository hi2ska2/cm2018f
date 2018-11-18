import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
######## Initial Condition Set #########
q_electron = 1.602192E-19		#C
m_electron = 9.109534e-34
k_B = 1.380662E-23				#J/K
Temp = 300						#K
V_D = 0.001						#V
tau = 1e-14						#sec

L_1 = 10e-9
L_2 = 10e-9
L_3 = 10e-9
del_x = 0.1e-9

total_width = L_1 + L_2 + L_3
N_grid_1 = int(L_1/del_x)
N_grid_2 = int(L_2/del_x)
N_grid_3 = int(L_3/del_x)
total_N_grid = N_grid_1 + N_grid_2 + N_grid_3 + 1
print("Total grid in the width = ",total_N_grid)
print("The index of first interface = ",N_grid_1+1)					# V[N_grid_1]
print("The index of second interface = ",N_grid_1+N_grid_2+1)		# V[N_grid_1+N_grid_2]

H_grid = np.logspace(-6,-3,num=30)					#V
#H_grid = np.logspace(0.00001,0.1,num=30)					#V
V_vector = np.block([
		np.zeros(N_grid_1),
		np.linspace(0,V_D, num=N_grid_2+1),
		V_D*np.ones(N_grid_3)
		])

f_0_plot = []
f_1_plot = []
for H in H_grid:		
	f_source = np.sqrt(2*np.pi)/(1+np.exp(q_electron*H/(k_B*Temp)))
	f_drain = np.sqrt(2*np.pi)/(1+np.exp(q_electron*(H+V_D)/(k_B*Temp)))

	C_left = H + 0.5*(V_vector[0:-2]+V_vector[1:-1])
	C_right = H + 0.5*(V_vector[1:-1]+V_vector[2:])
	Matrix = np.zeros((total_N_grid,total_N_grid))
	Matrix[0,0],Matrix[total_N_grid-1,total_N_grid-1] = 1,1
	Matrix[1:-1,0:-2] += C_left*np.eye(total_N_grid-2)
	Matrix[1:-1,2:] += C_right*np.eye(total_N_grid-2)
	Matrix[1:-1,1:-1] +=  -(C_left+C_right)*np.eye(total_N_grid-2) 

	Vector = np.zeros(total_N_grid)
	Vector[0] = f_source
	Vector[-1] = f_drain
	## result
	f_0 = np.matmul(np.linalg.inv(Matrix),Vector)
	f_0_plot.append(f_0)
	f_1 = -tau*np.sqrt(2*q_electron*(H+V_vector)/m_electron)*np.gradient(f_0,del_x)
	f_1_plot.append(f_1)

## Plot
fig = plt.figure()
fontsize=12

x_grid = np.linspace(0, total_width, total_N_grid)
ax = fig.gca(projection='3d')
x_grid, H_grid = np.meshgrid(x_grid, H_grid)

sulf=ax.plot_surface(x_grid*1e9, np.log10(H_grid) , np.array(f_1_plot), cmap=cm.plasma, linewidth=2, antialiased=False, alpha=0.7)
ax.set_xlabel("X (nm)",fontsize=fontsize)
ax.set_ylabel("log10(H) ",fontsize=fontsize)
ax.set_zlabel("$\mathrm{f_1}$",fontsize=fontsize)
#ax.set_zlim(0.0, 1.0)
#plt.colorbar( sulf, extend='both', shrink=0.7,boundaries=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])

plt.show()