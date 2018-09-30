###### input values ######
width_1 = 0.5       # unit of nm
epsilon_1 = 3.9
N_acc_1 = 0
n_intrin_1 = 0			# m^-3
width_2 = 5.0       # unit of nm
epsilon_2 = 11.7
N_acc_2 = 1.0E24	# m^-3
n_intrin_2 = 1.075e16			# m^-3
width_3 = 0.5       # unit of nm
epsilon_3 = 3.9
N_acc_3 = 0
n_intrin_3 = 0			# m^-3
del_x = 0.1	       # unit of nm
##########################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
q_electron = 1.602192E-19		#C
epsilon_0 = 8.854187817E-12		#F/m
k_B = 1.380662E-23				#J/K
Temp = 300						#K

total_width = width_1 + width_2 + width_3
N_grid_1 = int(width_1/del_x)
N_grid_2 = int(width_2/del_x)
N_grid_3 = int(width_3/del_x)
total_N_grid = N_grid_1 + N_grid_2 + N_grid_3 + 1
print("Total width of the capacitor = ",total_width,"nm")
print("Total grid in the width = ",total_N_grid)
x_grid = np.linspace(0, total_width, total_N_grid)

_matrix_1 = np.block([
        [np.array([1]),np.zeros(N_grid_1)],
        [np.transpose([np.block([np.array([epsilon_1]),np.zeros(N_grid_1-2)])]),-2*epsilon_1*np.eye(N_grid_1-1, dtype=int) + epsilon_1*np.eye(N_grid_1-1, k=1, dtype=int) + epsilon_1*np.eye(N_grid_1-1, k=-1, dtype=int),np.transpose([np.block([np.zeros(N_grid_1-2),np.array([epsilon_1])])])]
		])
_matrix_juction_12 = np.block([np.zeros(N_grid_1-1),np.array([epsilon_1,-epsilon_1-epsilon_2,epsilon_2]),np.zeros(N_grid_2 + N_grid_3-1)])
_matrix_2 = np.block([
        [np.transpose([np.block([np.array([epsilon_2]),np.zeros(N_grid_2-2)])]),-2*epsilon_2*np.eye(N_grid_2-1, dtype=int) + epsilon_2*np.eye(N_grid_2-1, k=1, dtype=int) + epsilon_2*np.eye(N_grid_2-1, k=-1, dtype=int),np.transpose([np.block([np.zeros(N_grid_2-2),np.array([epsilon_2])])])], 
		])
_matrix_juction_23 = np.block([np.zeros(N_grid_1+N_grid_2-1),np.array([epsilon_2,-epsilon_2-epsilon_3,epsilon_3]),np.zeros(N_grid_3-1)])
_matrix_3 = np.block([
        [np.transpose([np.block([np.array([epsilon_3]),np.zeros(N_grid_3-2)])]),-2*epsilon_3*np.eye(N_grid_3-1, dtype=int) + epsilon_3*np.eye(N_grid_3-1, k=1, dtype=int) + epsilon_3*np.eye(N_grid_3-1, k=-1, dtype=int),np.transpose([np.block([np.zeros(N_grid_3-2),np.array([epsilon_3])])])], 
		[np.zeros(N_grid_3),np.array([1])]
		])
Matrix = np.block([
	    [_matrix_1, np.zeros((N_grid_1,N_grid_2+N_grid_3))],
	    [_matrix_juction_12],
		[np.zeros((N_grid_2-1,N_grid_1)), _matrix_2, np.zeros((N_grid_2-1,N_grid_3))],
	    [_matrix_juction_23],
		[np.zeros((N_grid_3,N_grid_1+N_grid_2)), _matrix_3],
		])

fig = plt.figure(figsize=(9,5))
gs = GridSpec(6, 10)
subplot_phi = fig.add_subplot(gs[1:3, 0:4])
subplot_n = fig.add_subplot(gs[1:3, 6:10])
subplot_d_phi = fig.add_subplot(gs[4:6, 0:4])
subplot_d_phi_V = fig.add_subplot(gs[4:6, 6:10])
del_phi_max = []
for V_gate in np.linspace(1,0,11):
	phi_boundary= 0.33374-V_gate		#V
	rho_1 = (del_x**2)*(1e-18)*q_electron*N_acc_1/epsilon_0
	rho_2 = (del_x**2)*(1e-18)*q_electron*N_acc_2/epsilon_0
	rho_3 = (del_x**2)*(1e-18)*q_electron*N_acc_3/epsilon_0
	vector = np.block([
			np.array([phi_boundary]),
			rho_1*np.ones(N_grid_1-1),
			np.array([(rho_2+rho_1)/2]),
			rho_2*np.ones(N_grid_2-1),
			np.array([(rho_2+rho_3)/2]),
			rho_3*np.ones(N_grid_3-1),
			np.array([phi_boundary])])

	phi_init = np.matmul(np.linalg.inv(Matrix),vector)
	n_intrinsic_vector = np.block([
			n_intrin_1*np.ones(N_grid_1),
			np.array([(n_intrin_1+n_intrin_2)/2]),
			n_intrin_2*np.ones(N_grid_2-1),
			np.array([(n_intrin_2+n_intrin_3)/2]),
			n_intrin_3*np.ones(N_grid_3)
			])
	n_carrier = n_intrinsic_vector * np.exp(q_electron*phi_init/(k_B*Temp))

	# step for updated potential 
	del_rho = (del_x**2)*(1e-18)*q_electron*n_carrier/epsilon_0
	vector = vector + del_rho

	phi_updat = np.matmul(np.linalg.inv(Matrix),vector)
	del_phi_max = np.append(del_phi_max,np.amax(phi_init-phi_updat))

del_phi_0 = phi_init-phi_updat

subplot_phi.plot(x_grid, phi_init, 'b', label = "initial")
subplot_phi.plot(x_grid, phi_updat,  'r', label = "updated")
subplot_phi.set_xlabel("Width (nm)",fontsize=11)
subplot_phi.set_ylabel('$\mathrm{\phi}$ (V)' ,fontsize=11)
subplot_phi.grid(True)
subplot_phi.tick_params(which='both',direction='in')
subplot_phi.legend()

subplot_n.plot(x_grid, n_carrier, 'b')
subplot_n.set_xlabel("Width (nm)",fontsize=11)
subplot_n.set_ylabel(' n ($\mathrm{m^{-3}}$)' ,fontsize=11)
subplot_n.grid(True)
subplot_n.tick_params(which='both',direction='in')

subplot_d_phi.plot(x_grid, del_phi_0, 'b')
subplot_d_phi.set_xlabel("Width (nm)",fontsize=11)
subplot_d_phi.set_ylabel(' $\mathrm{\Delta\phi}$ (V)' ,fontsize=11)
subplot_d_phi.grid(True)
subplot_d_phi.tick_params(which='both',direction='in')


subplot_d_phi_V.plot(np.linspace(1,0,11), del_phi_max, 'b')
subplot_d_phi_V.set_xlabel("$\mathrm{V_{gate}}$ (V)",fontsize=11)
subplot_d_phi_V.set_ylabel(' max($\mathrm{\Delta\phi}$) (V)' ,fontsize=11)
subplot_d_phi_V.grid(True)
subplot_d_phi_V.tick_params(which='both',direction='in')

plt.suptitle("Double Gate MOS - fixed charge case",fontsize=20)
plt.savefig("./fixedCharge_result.png")
plt.show()
