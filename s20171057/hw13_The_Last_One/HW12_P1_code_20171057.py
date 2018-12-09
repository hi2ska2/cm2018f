###### input values ######
structure = 'long'
material = 'Si'
spacing = 0.5			# nm
##########################

if (structure =='long'):
	width_1 = 100e-9       # unit of m
	width_2 = 400e-9       # unit of m
	width_3 = 100e-9       # unit of m
	N_don_1 = 5e23
	N_don_2 = 2E21			# m^-3
	N_don_3 = 5e23
elif (structure =='short'):
	width_1 = 40e-9       # unit of m
	width_2 = 40e-9       # unit of m
	width_3 = 40e-9       # unit of m
	N_don_1 = 5e25
	N_don_2 = 2E23			# m^-3
	N_don_3 = 5e25
else :
	print("Check inputs : structure")
	exit()

if (material == 'Si'):
	n_intrin = 1.075e16
	epsilon = 11.7
else :
	print("Check inputs : material")
	exit()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
del_x = spacing * 1e-9
q_electron = 1.602192E-19		#C
epsilon_0 = 8.854187817E-12		#F/m
k_B = 1.380662E-23				#J/K
Temp = 300						#K
Volt_thermal = k_B*Temp/q_electron
coeff = del_x**2*q_electron/epsilon_0
err = 1e-20

total_width = width_1 + width_2 + width_3
N_grid_1 = int(width_1/del_x + err)
N_grid_2 = int(width_2/del_x + err)
N_grid_3 = int(width_3/del_x + err)
total_N_grid = N_grid_1 + N_grid_2 + N_grid_3 + 1
print("Total width of the capacitor = ",total_width,"m")
print("Total grid in the width = ",total_N_grid)
print("The index of first interface = ",N_grid_1+1)					# Vector[N_grid_1]
print("The index of second interface = ",N_grid_1+N_grid_2+1)		# Vector[N_grid_1+N_grid_2]
x_grid = np.linspace(0, total_width, total_N_grid)
eye = np.eye(total_N_grid-2)

def MakeVector(N_grid_1,N_grid_2,N_grid_3,v_1,v_2,v_3):
	vector = np.block([
			v_1*np.ones(N_grid_1),
			np.array([(v_1+v_2)/2]),
			v_2*np.ones(N_grid_2-1),
			np.array([(v_2+v_3)/2]),
			v_3*np.ones(N_grid_3)
			])
	return vector

N_don_vector = MakeVector(N_grid_1,N_grid_2,N_grid_3,N_don_1,N_don_2,N_don_3)
phi = Volt_thermal * np.log(N_don_vector/n_intrin)

Jaco_Lap = np.zeros((total_N_grid,total_N_grid))
Jaco_Lap[0,0],Jaco_Lap[total_N_grid-1,total_N_grid-1] = 1,1
Jaco_Lap[1:-1,0:-2] += epsilon*eye
Jaco_Lap[1:-1,2:] += epsilon*eye
Jaco_Lap[1:-1,1:-1] += -2*epsilon*eye

res = np.zeros(total_N_grid)

## nonlinaer Poisson solver
for i in range(0,10):
	n_electron = n_intrin*np.exp(phi/Volt_thermal)
	res[1:-1] = epsilon * (phi[0:-2]-2*phi[1:-1]+phi[2:]) - coeff*(-N_don_vector[1:-1]+n_electron[1:-1])
	Jaco = Jaco_Lap.copy()
	Jaco[1:-1,1:-1] -= coeff*n_electron[1:-1]*eye/Volt_thermal
	update = np.linalg.solve(Jaco,-res)
	phi = phi + update

n_electron = n_intrin*np.exp(phi/Volt_thermal)

phi_nlPoisson = phi.copy()
n_electron_nlPoisson = n_electron.copy()

## for x = [[phi],[elec]]
Jaco_init = np.zeros((2*total_N_grid,2*total_N_grid))
Jaco_init[0:total_N_grid,0:total_N_grid] = Jaco_Lap.copy()
Jaco_init[1:total_N_grid-1,total_N_grid+1:-1] = -coeff*eye
const = Volt_thermal/del_x
Jaco_init[total_N_grid,total_N_grid],Jaco_init[-1,-1] = 1,1
Jaco_init[total_N_grid+1:-1,total_N_grid:-2] += -1*const*eye
Jaco_init[total_N_grid+1:-1,total_N_grid+2:] += -1*const*eye
Jaco_init[total_N_grid+1:-1,total_N_grid+1:-1] += 2*const*eye
res = np.zeros(2*total_N_grid)

Diffusion = 0.01
Bernoulli = lambda x : x/(np.exp(x+err)-1)
CurrentDensity = lambda n_plus, n_minus, phi_Vt : -q_electron*Diffusion*(n_plus*Bernoulli(phi_Vt)-n_minus*Bernoulli(-phi_Vt))
V_grid = np.linspace(0.0, 0.5, 11)
del_V = V_grid[1] - V_grid[0] 
result = []
for V_appl in V_grid:
	phi[-1] = Volt_thermal * np.log(N_don_3/n_intrin) + V_appl
## Poisson-continuity solver
	for i in range(0,10):
		sum_n_left = (n_electron[0:-2]+n_electron[1:-1])
		sum_n_right = (n_electron[1:-1]+n_electron[2:])
		del_phi_left = (phi[1:-1]-phi[0:-2])
		del_phi_right = (phi[2:]-phi[1:-1])
		res[total_N_grid],res[-1]=n_electron[0]-N_don_1, n_electron[-1]-N_don_3
		res[1:total_N_grid-1] = epsilon * (phi[0:-2]-2*phi[1:-1]+phi[2:]) - coeff*(-N_don_vector[1:-1]+n_electron[1:-1])
		res[total_N_grid+1:-1] = sum_n_right*del_phi_right/(2*del_x) - sum_n_left*del_phi_left/(2*del_x) + const*(2*n_electron[1:-1]-n_electron[0:-2]-n_electron[2:])
		Jaco = Jaco_init.copy()
		Jaco[total_N_grid+1:-1,total_N_grid:-2] += (-1*del_phi_left/(2*del_x))*eye
		Jaco[total_N_grid+1:-1,total_N_grid+2:] += (del_phi_right/(2*del_x))*eye
		Jaco[total_N_grid+1:-1,total_N_grid+1:-1] += ((del_phi_right-del_phi_left)/(2*del_x))*eye
		Jaco[total_N_grid+1:-1,0:total_N_grid-2] += (sum_n_left/(2*del_x))*eye
		Jaco[total_N_grid+1:-1,2:total_N_grid] += (sum_n_right/(2*del_x))*eye
		Jaco[total_N_grid+1:-1,1:total_N_grid-1] += -1*((sum_n_left+sum_n_right)/(2*del_x))*eye
		update = np.linalg.solve(Jaco,-res)
		phi = phi+update[0:total_N_grid]
		n_electron = n_electron+update[total_N_grid: ]
	#result.append((sum_n_right*del_phi_right-Volt_thermal*(n_electron[2:]-n_electron[1:-1]))[-1])
	result.append(CurrentDensity(n_electron[-1],n_electron[-2],(phi[-1]-phi[-2])/Volt_thermal))


fig = plt.figure(figsize=(9,3))
gs = GridSpec(6, 10)
subplot_n = fig.add_subplot(gs[1:5, 0:4])
subplot_err = fig.add_subplot(gs[1:5, 6:10])

subplot_n.semilogy(x_grid*1e+9,n_electron_nlPoisson,color='gray',label ='Poisson')
subplot_n.semilogy(x_grid*1e+9,n_electron,'ro',label ='Self-consist')
subplot_n.set_xlabel("position (nm)",fontsize=11)
subplot_n.set_ylabel('$\mathrm{n_{e}}$ ($\mathrm{m^{-3}}$)' ,fontsize=11)
subplot_n.grid(True)
subplot_n.set_xlim(0,total_width*1e+9)
subplot_n.tick_params(which='both',direction='in')
subplot_n.set_title("(a) electron density distribution",multialignment='left')

subplot_err.plot(V_grid,result,color='red')
subplot_err.set_xlabel("Applied Voltage (V)",fontsize=11)
subplot_err.set_ylabel('Current / Area (J/$\mathrm{m^{-2}}$)' ,fontsize=11)
subplot_err.grid(True)
#subplot_err.set_xlim(0,total_width*1e+9)
subplot_err.tick_params(which='both',direction='in')
subplot_err.set_title("(b) Current density versus applied voltage")

plt.suptitle("Poisson-Continuity :"+structure,fontsize=20)
plt.savefig("./result_P1_"+structure+".png")
plt.show()
