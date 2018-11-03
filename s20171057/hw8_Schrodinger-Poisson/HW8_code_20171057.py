###### input values ######
L_x = 100e-9
L_y = 100e-9

L_z_1 = 0.5e-9
epsilon_1 = 3.9
N_acc_1 = 0
n_intrin_1 = 0

L_z_2 = 5e-9
epsilon_2 = 11.7
N_acc_2 = 1.0E24	# m^-3
n_intrin_2 = 1.075e16

L_z_3 = 0.5e-9
epsilon_3 = 3.9
N_acc_3 = 0
n_intrin_3 = 0
##########################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
 
def MakeVector(N_grid_1,N_grid_2,N_grid_3,n_1,n_2,n_3):
	n_vector = np.block([
			n_1*np.ones(N_grid_1),
			np.array([(n_1+n_2)/2]),
			n_2*np.ones(N_grid_2-1),
			np.array([(n_2+n_3)/2]),
			n_3*np.ones(N_grid_3)
			])
	return n_vector

fig = plt.figure(tight_layout=True,figsize=(9,8))
gs = GridSpec(10, 10)
subplot_n_C = fig.add_subplot(gs[1:5, 0:4])
subplot_n_Q = fig.add_subplot(gs[1:5, 5:9])
subplot_n_integ = fig.add_subplot(gs[5:9, 0:9])

h_bar = 6.626176e-34/(2*np.pi)
q_electron = 1.602192E-19		#C
m_0 = 9.109534e-31
k_B = 1.380662E-23				#J/K
Temp = 300						#K
epsilon_0 = 8.854187817e-12
Ec_Ei = 0.561004				# eV
Volt_thermal = k_B*Temp/q_electron

del_z = 0.1e-9
coeff_Poisson = del_z**2*q_electron/epsilon_0
total_width = L_z_1 + L_z_2 + L_z_3
N_grid_1 = int(L_z_1/del_z)
N_grid_2 = int(L_z_2/del_z)
N_grid_3 = int(L_z_3/del_z)
total_N_grid = N_grid_1 + N_grid_2 + N_grid_3 + 1
print("Total grid in the width = ",total_N_grid)
z_grid = np.linspace(0, total_width, total_N_grid)
m_matrix = np.array([[0.91, 0.19, 0.19],
					[0.19, 0.91, 0.19],
					[0.19, 0.19, 0.91]])

n_intrinsic_vector = MakeVector(N_grid_1,N_grid_2,N_grid_3,n_intrin_1,n_intrin_2,n_intrin_3)
N_acc_vector = MakeVector(N_grid_1,N_grid_2,N_grid_3,N_acc_1,N_acc_2,N_acc_3)

V_gate_init=0.0
phi = np.zeros(total_N_grid)
phi_boundary= 0.33374+V_gate_init		# V	
phi[0],phi[total_N_grid-1] = phi_boundary,phi_boundary	# initial potential	## nonlinear POISSON solver 

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

H_0 = -2*np.eye(N_grid_2-1) + np.eye(N_grid_2-1, k=1) + np.eye(N_grid_2-1, k=-1) 

V_gate_grid = np.linspace(0.0,1.0,21)
n_integrated_Classic = np.array([])
n_integrated_Quantum = np.array([])
E_tot_history = np.array([])
for V_gate in V_gate_grid:
	phi_boundary= 0.33374+V_gate		# V	
	phi[0],phi[total_N_grid-1] = phi_boundary,phi_boundary	# initial potential	## nonlinear POISSON solver 
	## Solving Poisson eq for initial potential
	for i in range(0,20):
		res = np.sum((jaco_Lap*phi), axis=1)
		res[0],res[total_N_grid-1] = 0.0, 0.0
		jaco = jaco_Lap.copy()
		###########################################
		## from charge correction
		electron_2 = n_intrinsic_vector[N_grid_1:total_N_grid-N_grid_3]*np.exp(phi[N_grid_1:total_N_grid-N_grid_3]/Volt_thermal)	# electron density in second area
		
		res[N_grid_1:total_N_grid-N_grid_3] -= coeff_Poisson*(N_acc_vector[N_grid_1:total_N_grid-N_grid_3]+electron_2)
		jaco[N_grid_1:total_N_grid-N_grid_3,N_grid_1:total_N_grid-N_grid_3] -= coeff_Poisson*electron_2/Volt_thermal*np.eye(N_grid_2+1)
		###########################################
		update = np.dot(np.linalg.inv(jaco),-1*res)
		phi += update

	electron_2 = n_intrinsic_vector[N_grid_1:total_N_grid-N_grid_3]*np.exp(phi[N_grid_1:total_N_grid-N_grid_3]/Volt_thermal)	# electron density in second area
	n_integrated_Classic = np.append(n_integrated_Classic,np.sum(electron_2)*del_z)
	subplot_n_C.semilogy(z_grid[N_grid_1:total_N_grid-N_grid_3]*1e9, electron_2*1e-6, color=(V_gate*0.9,0.2,0.8), label = "%.2f V" %V_gate)

	## Newton's step (Self-consistent loop)
	NBAND = 20
	n_elec=np.zeros(N_grid_2+1)
	E_older = 100
	del_E = 100
	while(del_E >1e-3):	## Energy cutoff : 1e-3 meV
		## Schrodinger eq solving >>
		n_total = 0
		E_total = 0
		n_z = 0
		potential = q_electron*Ec_Ei - q_electron*phi[N_grid_1+1:total_N_grid-N_grid_3-1]
		for valley in range(0,3):
			H_V = -2*m_matrix[valley,2]*m_0*(del_z/h_bar)**2*potential*np.eye(N_grid_2-1)
			Harmilton= H_0+H_V
			eigval, eigvec = np.linalg.eigh(Harmilton)
			E_z = np.flip(eigval/(-2*m_matrix[valley,2]*m_0*(del_z/h_bar)**2),0)[0:NBAND]
			band_square = np.flip(eigvec.T,0)[0:NBAND]**2/del_z

			coeff_Schro = 2*L_x*L_y/(2*np.pi)*np.sqrt(m_matrix[valley,0]*m_matrix[valley,1])*m_0/(h_bar**2)*k_B*Temp
			accupancy = coeff_Schro*np.log(1+np.exp(-E_z/(k_B*Temp)))

			E_total += np.sum(np.sum(accupancy*E_z))*6.242e+21
			density_band = 1/(L_x*L_y)*band_square
			n_z += np.matmul(accupancy , density_band)
		n_elec[1:N_grid_2] = n_z*2*2	# 2-fold degeneracy , Spin degeneracy
		print(E_total,"meV")
		## <<<<<<<<< end

		## Poisson eq solving >>
		res = np.sum((jaco_Lap*phi), axis=1)
		res[0],res[total_N_grid-1] = 0.0, 0.0
		jaco = jaco_Lap.copy()
		# from charge correction
		res[N_grid_1:total_N_grid-N_grid_3] -= coeff_Poisson*(N_acc_vector[N_grid_1:total_N_grid-N_grid_3]+n_elec)
		jaco[N_grid_1:total_N_grid-N_grid_3,N_grid_1:total_N_grid-N_grid_3] -= coeff_Poisson*n_elec/Volt_thermal*np.eye(N_grid_2+1)

		update = np.dot(np.linalg.inv(jaco),-1*res)
		phi += update
		del_E = np.abs(E_total - E_older)
		E_older = E_total
		## <<<<<<<<< end
	
	E_tot_history = np.append(E_tot_history,E_total)
	n_integrated_Quantum = np.append(n_integrated_Quantum,np.sum(n_elec)*del_z)
	subplot_n_Q.semilogy(z_grid[N_grid_1:total_N_grid-N_grid_3]*1e9, n_elec*1e-6, color=(V_gate*0.9,0.2,0.8), label = "%.2f V" %V_gate)
	## <<<<<<<<< Newton end
print(E_tot_history,"meV")
	

subplot_n_C.set_xlabel("position (nm)",fontsize=11)
subplot_n_C.set_ylabel('$\mathrm{n_{e}}$ ($\mathrm{cm^{-3}}$)' ,fontsize=11)
subplot_n_C.grid(True)
subplot_n_C.set_xlim(0.5,5.5)
subplot_n_C.set_ylim(1e13,1e21)
subplot_n_C.set_title("Electron density(Poisson)",fontsize=12)
subplot_n_C.legend(bbox_to_anchor=(1.0, 1.05), loc=2,fontsize=6.5)
subplot_n_C.tick_params(which='both',direction='in')

subplot_n_Q.set_xlabel("position (nm)",fontsize=11)
subplot_n_Q.set_ylabel('$\mathrm{n_{e}}$ ($\mathrm{cm^{-3}}$)' ,fontsize=11)
subplot_n_Q.grid(True)
subplot_n_Q.set_xlim(0.5,5.5)
subplot_n_Q.set_ylim(1e13,1e21)
subplot_n_Q.set_title("Electron density(Schro+Poiss)",fontsize=12)
subplot_n_Q.legend(bbox_to_anchor=(1.0, 1.05), loc=2,fontsize=6.5)
subplot_n_Q.tick_params(which='both',direction='in')

subplot_n_integ.semilogy(V_gate_grid, n_integrated_Classic*1.0e-4, 'k',label = "Poisson")
subplot_n_integ.semilogy(V_gate_grid, n_integrated_Quantum*1.0e-4, 'r',label = "Schro-Poiss")
subplot_n_integ.set_xlabel("$\mathrm{V_{Gate}} (V)$",fontsize=11)
subplot_n_integ.set_ylabel('$\mathrm{n_{integrated}}$ ($\mathrm{cm^{-2}}$)' ,fontsize=11)
subplot_n_integ.set_xlim(0,1)
subplot_n_integ.set_ylim(1e8,1e14)
subplot_n_integ.set_title("Integrated electron density",fontsize=12)
subplot_n_integ.legend(fontsize=10)
subplot_n_integ.grid(True)
subplot_n_integ.tick_params(which='both',direction='in')

plt.suptitle("Schrodinger-Poisson Solver",fontsize=20)
plt.savefig("./Schrodinger-Poisson_result.png")
plt.show()

