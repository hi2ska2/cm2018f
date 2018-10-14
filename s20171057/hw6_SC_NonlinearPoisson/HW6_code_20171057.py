###### input values ######
width_1 = 0.5e-9       # unit of m
epsilon_1 = 3.9
N_acc_1 = 0
n_intrin_1 = 0			# m^-3
width_2 = 5.0e-9       # unit of m
epsilon_2 = 11.7
N_acc_2 = 1.0E24	# m^-3
n_intrin_2 = 1.075e16			# m^-3
width_3 = 0.5e-9       # unit of m
epsilon_3 = 3.9
N_acc_3 = 0
n_intrin_3 = 0			# m^-3
del_x = 0.1e-9	       # unit of m
##########################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
q_electron = 1.602192E-19		#C
epsilon_0 = 8.854187817E-12		#F/m
k_B = 1.380662E-23				#J/K
Temp = 300						#K
Volt_thermal = k_B*Temp/q_electron
coeff = del_x**2*q_electron/epsilon_0

total_width = width_1 + width_2 + width_3
N_grid_1 = int(width_1/del_x)
N_grid_2 = int(width_2/del_x)
N_grid_3 = int(width_3/del_x)
total_N_grid = N_grid_1 + N_grid_2 + N_grid_3 + 1
print("Total width of the capacitor = ",total_width,"m")
print("Total grid in the width = ",total_N_grid)
x_grid = np.linspace(0, total_width, total_N_grid)

n_intrinsic_vector = np.block([
		n_intrin_1*np.ones(N_grid_1),
		np.array([(n_intrin_1+n_intrin_2)/2]),
		n_intrin_2*np.ones(N_grid_2-1),
		np.array([(n_intrin_2+n_intrin_3)/2]),
		n_intrin_3*np.ones(N_grid_3)
		])
N_acc_vector = np.block([
		N_acc_1*np.ones(N_grid_1),
		np.array([(N_acc_1+N_acc_2)/2]),
		N_acc_2*np.ones(N_grid_2-1),
		np.array([(N_acc_2+N_acc_3)/2]),
		N_acc_3*np.ones(N_grid_3)
		])


## from Laplacian equation
_jaco_1 = np.block([
        [np.array([1]),np.zeros(N_grid_1)],
        [np.transpose([np.block([np.array([epsilon_1]),np.zeros(N_grid_1-2)])]),-2*epsilon_1*np.eye(N_grid_1-1, dtype=int) + epsilon_1*np.eye(N_grid_1-1, k=1, dtype=int) + epsilon_1*np.eye(N_grid_1-1, k=-1, dtype=int),np.transpose([np.block([np.zeros(N_grid_1-2),np.array([epsilon_1])])])]
		])
_jaco_juction_12 = np.block([np.zeros(N_grid_1-1),np.array([epsilon_1,-epsilon_1-epsilon_2,epsilon_2]),np.zeros(N_grid_2 + N_grid_3-1)])
_jaco_2 = np.block([
        [np.transpose([np.block([np.array([epsilon_2]),np.zeros(N_grid_2-2)])]),-2*epsilon_2*np.eye(N_grid_2-1, dtype=int) + epsilon_2*np.eye(N_grid_2-1, k=1, dtype=int) + epsilon_2*np.eye(N_grid_2-1, k=-1, dtype=int),np.transpose([np.block([np.zeros(N_grid_2-2),np.array([epsilon_2])])])], 
		])
_jaco_juction_23 = np.block([np.zeros(N_grid_1+N_grid_2-1),np.array([epsilon_2,-epsilon_2-epsilon_3,epsilon_3]),np.zeros(N_grid_3-1)])
_jaco_3 = np.block([
        [np.transpose([np.block([np.array([epsilon_3]),np.zeros(N_grid_3-2)])]),-2*epsilon_3*np.eye(N_grid_3-1, dtype=int) + epsilon_3*np.eye(N_grid_3-1, k=1, dtype=int) + epsilon_3*np.eye(N_grid_3-1, k=-1, dtype=int),np.transpose([np.block([np.zeros(N_grid_3-2),np.array([epsilon_3])])])], 
		[np.zeros(N_grid_3),np.array([1])]
		])
jaco_Lap = np.block([
	    [_jaco_1, np.zeros((N_grid_1,N_grid_2+N_grid_3))],
	    [_jaco_juction_12],
		[np.zeros((N_grid_2-1,N_grid_1)), _jaco_2, np.zeros((N_grid_2-1,N_grid_3))],
	    [_jaco_juction_23],
		[np.zeros((N_grid_3,N_grid_1+N_grid_2)), _jaco_3],
		])

fig = plt.figure(figsize=(9,3))
gs = GridSpec(6, 10)
subplot_n = fig.add_subplot(gs[1:5, 0:4])
subplot_n_V = fig.add_subplot(gs[1:5, 6:10])

n_integrated = []
Vg_grid = np.linspace(-1,1,21)
for V_gate in Vg_grid:
	print(V_gate)
	phi_boundary= 0.33374-V_gate		# V
	phi = np.zeros(total_N_grid)	# initial potential	
	phi[0],phi[total_N_grid-1] = phi_boundary,phi_boundary	# initial potential	
	for i in range(0,1000):
		res = np.sum((jaco_Lap*phi), axis=1)
		res[0],res[total_N_grid-1] = 0,0
		jaco = jaco_Lap.copy()
		###########################################

		## from charge correction
		electron_2 = n_intrinsic_vector[N_grid_1:total_N_grid-N_grid_3]*np.exp(phi[N_grid_1:total_N_grid-N_grid_3]/Volt_thermal)	# electron density in second area
		
		res[N_grid_1:total_N_grid-N_grid_3] -= coeff*(N_acc_vector[N_grid_1:total_N_grid-N_grid_3]+electron_2)
		jaco[N_grid_1:total_N_grid-N_grid_3,N_grid_1:total_N_grid-N_grid_3] -= coeff*electron_2/Volt_thermal*np.eye(N_grid_2+1)
		###########################################
		update = np.matmul(np.linalg.inv(jaco),-1*res)
		phi += update
	
	subplot_n.semilogy(x_grid[N_grid_1:total_N_grid-N_grid_3], electron_2, label = "Vg="+str(V_gate))
	electron_2 = n_intrinsic_vector[N_grid_1:total_N_grid-N_grid_3]*np.exp(phi[N_grid_1:total_N_grid-N_grid_3]/Volt_thermal)	# final electron density
	n_integrated.append(np.trapz(electron_2, dx=del_x))

print(n_integrated)


subplot_n.set_xlabel("position (m)",fontsize=11)
subplot_n.set_ylabel('$\mathrm{n_{e}}$ ($\mathrm{m^{-3}}$)' ,fontsize=11)
subplot_n.grid(True)
subplot_n.set_xlim(0.5e-9,5.5e-9)
subplot_n.tick_params(which='both',direction='in')

subplot_n_V.semilogy(Vg_grid, n_integrated, 'ko')
subplot_n_V.set_xlabel("$\mathrm{V_{G}} (V)$",fontsize=11)
subplot_n_V.set_ylabel('$\mathrm{n_{integrated}}$ ($\mathrm{m^{-2}}$)' ,fontsize=11)
subplot_n_V.grid(True)
subplot_n_V.tick_params(which='both',direction='in')

plt.suptitle("Double Gate MOS : SC-result",fontsize=20)
plt.savefig("./DoubleGate_SC_result.png")
plt.show()
