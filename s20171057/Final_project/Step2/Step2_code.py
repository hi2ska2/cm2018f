import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
sys.path.append('../')
from Constants import *
_jaco=np.load('./jaco.npy')
######## Structure Set #########
V_gate = float(sys.argv[1])
V_source = 0.0
V_drain_grid = np.linspace(0,0.1,2)
################################
'''    model discription
							^ z
			 				  > y 
			Oxide(Boundary)
--------------------------
  Si_1  |  Si_2  |  Si_3  		<= I just consider this area for continueity eq.
--------------------------
			Oxide(Boundary)

'''
mobility = 1430e-4		# m^2/V
Bernoulli = lambda x : 0 if x==0 else x/(np.exp(x)-1)
CurrentDensity = lambda n_plus, n_minus, phi_Vt, delta : (n_plus*Bernoulli(phi_Vt)-n_minus*Bernoulli(-phi_Vt))*Volt_thermal/delta
current_result = []
for V_drain in V_drain_grid:
	## Poisson Solver to get Boundary Conditions for Drift-Diffusion solver
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
				phi[y,z] = Volt_thermal*np.log(N_don_yvec[y]/n_intrin) + V_drain

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
	n_elec[:,N_Oxi_z:-N_Oxi_z] = n_intrin*np.exp(phi[:,N_Oxi_z:-N_Oxi_z]/Volt_thermal)
	n_elec[:,N_Oxi_z]=n_elec[:,N_Oxi_z]
	n_elec[:,-N_Oxi_z-1]=n_elec[:,-N_Oxi_z-1]
	phi_ref = phi.copy()

	## Drift-Diffusion Solver 
	'''    model discription
								^ z
								> y 
	----------- BC ------------
	Si_1  |  Si_2  |  Si_3  		<= I just consider Si area for continueity eq.
	----------- BC ------------
	'''
	_res_phi = np.zeros((N_y,N_Si_z+1))
	_res_ne = np.zeros((N_y,N_Si_z+1))
	_Jaco_phi_phi = _jaco[N_Oxi_z*N_y:N_z*N_y-N_Oxi_z*N_y,:,N_Oxi_z:N_z-N_Oxi_z].copy()
	_Jaco_phi_ne = np.zeros((N_y*(N_Si_z+1),N_y,N_Si_z+1))
	_Jaco_ne_ne_const = np.zeros((N_y*(N_Si_z+1),N_y,N_Si_z+1))
	_Jaco_ne_phi = np.zeros((N_y*(N_Si_z+1),N_y,N_Si_z+1))
	rate = 1e-20
	i=0
	for z in range(N_Si_z+1):
		for y in range(N_y):
			if(z==0 or z==N_Si_z):	# Ox-Si interfaces
				_Jaco_ne_ne_const[i,y,z] = 1
				_Jaco_phi_phi[i] = 0
				_Jaco_phi_phi[i,y,z] = 1
			elif(y==0):	# Source
				phi[y,z] = Volt_thermal*np.log(N_don_1/n_intrin) + V_source
				_Jaco_ne_ne_const[i,y,z] = 1
				_Jaco_phi_phi[i] = 0
				_Jaco_phi_phi[i,y,z] = 1
			elif(y==N_y-1):	# Drain
				phi[y,z] = Volt_thermal*np.log(N_don_3/n_intrin) + V_drain
				_Jaco_ne_ne_const[i,y,z] = 1
				_Jaco_phi_phi[i] = 0
				_Jaco_phi_phi[i,y,z] = 1
			else:
				_Jaco_ne_ne_const[i,y,z] = 2*Volt_thermal/del_width + 2*Volt_thermal/del_hight
				_Jaco_ne_ne_const[i,y-1,z] = -1*Volt_thermal/del_width
				_Jaco_ne_ne_const[i,y+1,z] = -1*Volt_thermal/del_width
				_Jaco_ne_ne_const[i,y,z-1] = -1*Volt_thermal/del_hight
				_Jaco_ne_ne_const[i,y,z+1] = -1*Volt_thermal/del_hight
				_Jaco_phi_ne[i,y,z] = -coeff
			i=i+1
	## Newton step
	Jaco = np.zeros((2*N_y*(N_Si_z+1),2*N_y*(N_Si_z+1))) 
	Res = np.zeros(2*N_y*(N_Si_z+1))
	switch = True
	#while(switch):
	for i in range(3):
		phi_Si = phi[:,N_Oxi_z:N_z-N_Oxi_z]
		ne_Si = n_elec[:,N_Oxi_z:N_z-N_Oxi_z]
		_Jaco_ne_ne = _Jaco_ne_ne_const.copy()
		## make Jacobian
		i=0
		for z in range(N_Si_z+1):
			for y in range(N_y):
				if(z==0 or z==N_Si_z):	# Ox-Si interfaces
					_res_phi[y,z] = 0
					_res_ne[y,z] = 0
				elif(y==0):	# Source
					_res_phi[y,z] = phi[y,z] - Volt_thermal*np.log(N_don_1/n_intrin) - V_source
					_res_ne[y,z] = ne_Si[y,z]-N_don_1
				elif(y==N_y-1):	# Drain
					_res_phi[y,z] = phi[y,z] - Volt_thermal*np.log(N_don_3/n_intrin) - V_drain
					_res_ne[y,z] = ne_Si[y,z]-N_don_3
				else: 
					_res_phi[y,z] = np.sum(_Jaco_phi_phi[i]*phi_Si) - coeff*(-N_don_yvec[y]+ne_Si[y,z])
					#_res_ne[y,z] = (ne_Si[y,z]+ne_Si[y+1,z])*(phi_Si[y+1,z]-phi_Si[y,z])/(2*del_width) - (ne_Si[y,z]+ne_Si[y-1,z])*(phi_Si[y,z]-phi_Si[y-1,z])/(2*del_width) + Volt_thermal*(2*ne_Si[y,z]-ne_Si[y-1,z]-ne_Si[y+1,z])/del_width
					#_res_ne[y,z] += (ne_Si[y,z]+ne_Si[y,z+1])*(phi_Si[y,z+1]-phi_Si[y,z])/(2*del_hight) - (ne_Si[y,z]+ne_Si[y,z-1])*(phi_Si[y,z]-phi_Si[y,z-1])/(2*del_hight) + Volt_thermal*(2*ne_Si[y,z]-ne_Si[y,z-1]-ne_Si[y,z+1])/del_hight
					_res_ne[y,z] = -CurrentDensity(ne_Si[y+1,z],ne_Si[y,z],(phi_Si[y+1,z]-phi_Si[y,z])/Volt_thermal,del_width) +CurrentDensity(ne_Si[y,z],ne_Si[y-1,z],(phi_Si[y,z]-phi_Si[y-1,z])/Volt_thermal,del_width)
					_res_ne[y,z] += -CurrentDensity(ne_Si[y,z+1],ne_Si[y,z],(phi_Si[y,z+1]-phi_Si[y,z])/Volt_thermal,del_hight) +CurrentDensity(ne_Si[y,z],ne_Si[y,z-1],(phi_Si[y,z]-phi_Si[y,z-1])/Volt_thermal,del_hight)
					_Jaco_ne_ne[i,y,z] += (phi_Si[y+1,z]-2*phi_Si[y,z]+phi_Si[y-1,z])/(2*del_width) + (phi_Si[y,z+1]-2*phi_Si[y,z]+phi_Si[y,z-1])/(2*del_hight)
					_Jaco_ne_ne[i,y-1,z] += (phi_Si[y-1,z]-phi_Si[y,z])/(2*del_width)
					_Jaco_ne_ne[i,y+1,z] += (phi_Si[y+1,z]-phi_Si[y,z])/(2*del_width)
					_Jaco_ne_ne[i,y,z-1] += (phi_Si[y,z-1]-phi_Si[y,z])/(2*del_hight)
					_Jaco_ne_ne[i,y,z+1] += (phi_Si[y,z-1]-phi_Si[y,z])/(2*del_hight)
					_Jaco_ne_phi[i,y,z] = -1*(2*ne_Si[y,z]+ne_Si[y-1,z]+ne_Si[y+1,z])/(2*del_width) -1*(2*ne_Si[y,z]+ne_Si[y,z-1]+ne_Si[y,z+1])/(2*del_hight)
					_Jaco_ne_phi[i,y-1,z] = (ne_Si[y,z]+ne_Si[y-1,z])/(2*del_width) 
					_Jaco_ne_phi[i,y+1,z] = (ne_Si[y,z]+ne_Si[y+1,z])/(2*del_width) 
					_Jaco_ne_phi[i,y,z-1] = (ne_Si[y,z]+ne_Si[y,z-1])/(2*del_hight) 
					_Jaco_ne_phi[i,y,z+1] = (ne_Si[y,z]+ne_Si[y,z+1])/(2*del_hight)
					_Jaco_ne_ne[i] = _Jaco_ne_ne[i]*rate
					_Jaco_ne_phi[i] = _Jaco_ne_phi[i]*rate
					_res_ne[y,z] = _res_ne[y,z]*rate
				i=i+1
		Jaco[0:N_y*(N_Si_z+1),0:N_y*(N_Si_z+1)] = np.reshape(np.swapaxes(_Jaco_phi_phi,1,2),(-1,N_y*(N_Si_z+1)))
		Jaco[0:N_y*(N_Si_z+1),N_y*(N_Si_z+1): ] = np.reshape(np.swapaxes(_Jaco_phi_ne,1,2),(-1,N_y*(N_Si_z+1)))
		Jaco[N_y*(N_Si_z+1): ,0:N_y*(N_Si_z+1)] = np.reshape(np.swapaxes(_Jaco_ne_phi,1,2),(-1,N_y*(N_Si_z+1)))
		Jaco[N_y*(N_Si_z+1): ,N_y*(N_Si_z+1): ] = np.reshape(np.swapaxes(_Jaco_ne_ne,1,2),(-1,N_y*(N_Si_z+1)))
		Res[0:N_y*(N_Si_z+1)] = _res_phi.T.flatten()
		Res[N_y*(N_Si_z+1): ] = _res_ne.T.flatten()
		update = np.linalg.solve(Jaco,-Res)
		update_phi = update[0:N_y*(N_Si_z+1)]*0.1
		update_ne = update[N_y*(N_Si_z+1):]*0.1
		phi[:,N_Oxi_z:N_z-N_Oxi_z] = phi[:,N_Oxi_z:N_z-N_Oxi_z] + update_phi.reshape((N_Si_z+1),N_y).T
		n_elec[:,N_Oxi_z:N_z-N_Oxi_z] = n_elec[:,N_Oxi_z:N_z-N_Oxi_z] + update_ne.reshape((N_Si_z+1),N_y).T
		norm_update_phi = np.linalg.norm(update_phi)
		norm_update_ne = np.linalg.norm(update_ne)
		print(norm_update_phi,norm_update_ne)
		if(norm_update_phi<1e-5):
			switch = False
		if((V_gate==0.2 or V_gate==0.5) and (V_drain==0.2 or V_drain==0.5)):
			np.savetxt('Step2_Phi_result_Vg'+str(V_gate)+'Vd'+str(V_drain)+'.out',phi.T,fmt="%1.2e")
	current=0.0
	for z in range(N_Si_z+1):
		current += -q_electron*mobility*CurrentDensity(n_elec[-1,N_Oxi_z+z],n_elec[-2,N_Oxi_z+z],(phi[-1,N_Oxi_z+z]-phi[-2,N_Oxi_z+z])/Volt_thermal,del_width)
	current_result.append((current / N_Si_z+1)*hight_Si*1e-9)

#np.savetxt('Step2_VI_result_Vg'+str(V_gate)+'.out',[[V_drain_grid.T],[np.array(current_result).T]])

## Plot_Phi
fig = plt.figure()
fontsize=12

ax = fig.gca(projection='3d')
Y_grid = np.linspace(1, N_y, N_y)
Z_grid = np.linspace(1, N_z, N_z)
Y_grid, Z_grid = np.meshgrid(Y_grid, Z_grid)

sulf=ax.plot_surface(Y_grid, Z_grid, phi.T, linewidth=2, antialiased=False, alpha=0.7)
sulf=ax.plot_surface(Y_grid, Z_grid, phi_ref.T, linewidth=2, antialiased=False, alpha=0.7)
ax.set_xlabel("Y",fontsize=fontsize)
ax.set_ylabel("Z",fontsize=fontsize)
ax.set_zlabel("$\mathrm{\Phi}$",fontsize=fontsize)

fig.savefig('Step1_Phi_result_'+str(V_gate)+'.png')
plt.show()

## Plot_n_elec
fig = plt.figure()
fontsize=12

ax = fig.gca(projection='3d')
Y_grid = np.linspace(1, N_y, N_y)
Z_grid = np.linspace(1, N_z, N_z)
Y_grid, Z_grid = np.meshgrid(Y_grid, Z_grid)

sulf=ax.plot_surface(Y_grid, Z_grid, n_elec.T, linewidth=2, antialiased=False, alpha=0.7)
ax.set_xlabel("Y",fontsize=fontsize)
ax.set_ylabel("Z",fontsize=fontsize)
ax.set_zlabel("$\mathrm{n_e}$",fontsize=fontsize)

fig.savefig('Step1_ne_result_'+str(V_gate)+'.png')