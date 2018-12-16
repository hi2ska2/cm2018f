import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
sys.path.append('../')
from Constants import *
_jaco=np.load('./jaco.npy')
######## Voltage Set #########
V_gate = float(sys.argv[1])
V_source = 0.0
V_drain = 0.0
################################
'''    model discription
							^ z
			 				  > y 
			Oxide
--------------------------
  Si_1  |  Si_2  |  Si_3  
--------------------------
			Oxide

'''

phi = np.zeros((N_y,N_z))
_res = np.zeros((N_y,N_z))
for z in range(N_z):
	for y in range(N_y):
		if(y==0 or y==N_y-1 or z==0 or z==N_z-1):	#Boundary condition
			if((z==0 or z==N_z-1) and (y>=N_Si_1_y and y<N_y-N_Si_3_y)):	# Gate
				phi[y,z] = 0.33374 + V_gate
			elif(y==0 and (z>=N_Oxi_z and z<N_z-N_Oxi_z)):	# Source
				phi[y,z] = Volt_thermal*np.log(N_don_1/n_intrin) + V_source
			elif(y==N_y-1 and (z>=N_Oxi_z and z<N_z-N_Oxi_z)):	# Drain
				phi[y,z] = Volt_thermal*np.log(N_don_3/n_intrin) + V_drain
		else: 
			phi[y,z] = Volt_thermal*np.log(N_don_yvec[y]/n_intrin) + V_source

## Newton step 
switch = True
n_elec = np.zeros((N_y,N_z))
while(switch):
	n_elec[:,N_Oxi_z:-N_Oxi_z] = n_intrin*np.exp(phi[:,N_Oxi_z:-N_Oxi_z]/Volt_thermal)
	_Jaco = _jaco.copy()
	## make Jacobian
	i=0
	for z in range(N_z):
		for y in range(N_y):
			if((z==0 or z==N_z-1) and (y>=N_Si_1_y and y<N_y-N_Si_3_y)):	# Gate
				pass
			elif(y==0 and (z>=N_Oxi_z and z<N_z-N_Oxi_z)):	# Source
				pass
			elif(y==N_y-1 and (z>=N_Oxi_z and z<N_z-N_Oxi_z)):	# Drain
				pass
			else: 
				_res[y,z] = np.sum(_jaco[i]*phi)
				if(z>=N_Oxi_z+1 and z<=N_z-N_Oxi_z-2):		# Si region
					_res[y,z] -= coeff*(-N_don_yvec[y]+n_elec[y,z])
					_Jaco[i,y,z] -= coeff*n_elec[y,z]/Volt_thermal
				elif(z==N_Oxi_z or z == N_z-N_Oxi_z-1):			# interface
					_res[y,z] -= coeff*(-N_don_yvec[y]+n_elec[y,z])*0.5
					_Jaco[i,y,z] -= 0.5*coeff*n_elec[y,z]/Volt_thermal
			i=i+1
	Jaco = np.reshape(np.swapaxes(_Jaco,1,2),(-1,N_y*N_z))
	Res = _res.T.flatten()
	update = np.linalg.solve(Jaco,-Res)
	update = update.reshape(N_z,N_y).T
	phi = phi + update
	norm_update = np.linalg.norm(update,np.inf)
	print(norm_update)
	if(norm_update<1e-5):
		switch = False

n_elec = np.zeros((N_y,N_z))
n_elec[:,N_Oxi_z+1:-N_Oxi_z-1] = n_intrin*np.exp(phi[:,N_Oxi_z+1:-N_Oxi_z-1]/Volt_thermal)
n_elec[:,N_Oxi_z]=n_elec[:,N_Oxi_z]/2
n_elec[:,-N_Oxi_z-1]=n_elec[:,-N_Oxi_z-1]/2
np.savetxt('Step1_Phi_result_'+str(V_gate)+'.out',phi,fmt="%e")
np.savetxt('Step1_ne_result_'+str(V_gate)+'.out',n_elec,fmt="%e")
## Plot_Phi
fig = plt.figure()
fontsize=12

ax = fig.gca(projection='3d')
Y_grid = np.linspace(0, total_width, N_y)
Z_grid = np.linspace(0, total_hight, N_z)
Y_grid, Z_grid = np.meshgrid(Y_grid, Z_grid)

sulf=ax.plot_surface(Y_grid, Z_grid, phi.T, linewidth=2, antialiased=False, alpha=0.7)
ax.set_xlabel("Width (nm)",fontsize=fontsize)
ax.set_ylabel("High (nm)",fontsize=fontsize)
ax.set_zlabel("$\mathrm{\Phi}$",fontsize=fontsize)

fig.savefig('Step1_Phi_result_'+str(V_gate)+'.png')

## Plot_n_elec
fig = plt.figure()
fontsize=12

ax = fig.gca(projection='3d')
Y_grid = np.linspace(0, total_width, N_y)
Z_grid = np.linspace(0, total_hight, N_z)
Y_grid, Z_grid = np.meshgrid(Y_grid, Z_grid)

sulf=ax.plot_surface(Y_grid, Z_grid, n_elec.T, linewidth=2, antialiased=False, alpha=0.7)
ax.set_xlabel("Width (nm)",fontsize=fontsize)
ax.set_ylabel("High (nm)",fontsize=fontsize)
ax.set_zlabel("$\mathrm{n_e}$",fontsize=fontsize)

fig.savefig('Step1_ne_result_'+str(V_gate)+'.png')