###### input values ######
Temp = 300						#K				
n_intrin = 1e10			# cm^-3	
order_N_plus_init = 10			# cm^-3
order_N_plus_fin = 18			# cm^-3
##########################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
q_electron = 1					# electron coulomb
k_B = 8.6173303E-5				# eV/K
Volt_Thermal = (k_B*Temp)/q_electron	#V

N_log = []
phi_log = []
for N_plus in np.logspace(order_N_plus_init,order_N_plus_fin,27):
	fix_order = np.log10(N_plus)+10
	N_plus_fix = N_plus*np.exp(-fix_order)

	function = lambda phi : N_plus_fix + n_intrin*(np.exp(-phi/Volt_Thermal-fix_order)-np.exp(phi/Volt_Thermal-fix_order))
	Jaco_function = lambda phi : (n_intrin/Volt_Thermal)*(-np.exp(-phi/Volt_Thermal-fix_order)-np.exp(phi/Volt_Thermal-fix_order))
	phi = 10.00
	update = -function(phi)/Jaco_function(phi)
	while(np.abs(update)>1.0e-8):
		phi = phi+update
		update = -function(phi)/Jaco_function(phi)
	
	N_log = np.append(N_log,N_plus)
	phi_log = np.append(phi_log,phi)

for N_plus in -1*np.logspace(order_N_plus_init,order_N_plus_fin,27):
	fix_order = np.log10(-N_plus)+10
	N_plus_fix = N_plus*np.exp(fix_order)

	function = lambda phi : N_plus_fix + n_intrin*(np.exp(-phi/Volt_Thermal+fix_order)-np.exp(phi/Volt_Thermal+fix_order))
	Jaco_function = lambda phi : (n_intrin/Volt_Thermal)*(-np.exp(-phi/Volt_Thermal+fix_order)-np.exp(phi/Volt_Thermal+fix_order))
	phi = -10.00
	update = -function(phi)/Jaco_function(phi)
	while(np.abs(update)>1.0e-8):
		phi = phi+update
		update = -function(phi)/Jaco_function(phi)
	
	N_log = np.append(N_log,N_plus)
	phi_log = np.append(phi_log,phi)


## Analytic solution
N_grid = np.logspace(order_N_plus_init,order_N_plus_fin,27)
Solution = lambda N : Volt_Thermal * np.arcsinh(N/(2*n_intrin))

## plot
fig = plt.figure(figsize=(9,3))
gs = GridSpec(6, 10)
subplot_phi = fig.add_subplot(gs[1:5, 0:4])
subplot_err = fig.add_subplot(gs[1:5, 6:10])

subplot_phi.semilogx(N_log, phi_log, 'bo', label = "Positive")
subplot_phi.semilogx(-N_log, phi_log, 'ro', label = "Negative")
subplot_phi.semilogx(N_grid, Solution(N_grid), 'k', label = "Analytic")
subplot_phi.semilogx(N_grid, Solution(-N_grid), 'k')
subplot_phi.set_xlabel("$\mathrm{N^+}$ ($\mathrm{cm^{-3}}$)",fontsize=11)
subplot_phi.set_ylabel('$\mathrm{\phi}$ (V)' ,fontsize=11)
subplot_phi.grid(True)
subplot_phi.set_ylim(-0.6,0.6)
subplot_phi.set_xlim(1e10,1e18)
subplot_phi.tick_params(which='both',direction='in')
subplot_phi.legend()

subplot_err.semilogx(N_log, (Solution(N_log)-phi_log)/Solution(N_log), 'ko')
subplot_err.set_xlabel("$\mathrm{N^+}$ ($\mathrm{cm^{-3}}$)",fontsize=11)
subplot_err.set_ylabel(' error (%)' ,fontsize=11)
subplot_err.set_xlim(1e10,1e18)
subplot_err.grid(True)
subplot_err.tick_params(which='both',direction='in')

plt.suptitle("Newton's method",fontsize=20)
plt.savefig("./Newton_result.png")
plt.show()
