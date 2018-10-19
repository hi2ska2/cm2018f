###### input values ######
L_x = 100e-9
m_xx = 0.19
L_y = 100e-9
m_yy = 0.19

L_z = 5e-9
m_zz = 0.91
##########################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(9,3))
gs = GridSpec(6, 10)
subplot_n_z = fig.add_subplot(gs[1:5, 0:4])
subplot_n_Ef = fig.add_subplot(gs[1:5, 6:10])

h_bar = 6.626176e-34/(2*np.pi)
q_electron = 1.602192E-19		#C
m_0 = 9.109534e-31
k_B = 1.380662E-23				#J/K
Temp = 300						#K
Volt_thermal = k_B*Temp/q_electron
coeff = 2*L_x*L_y/(2*np.pi)*np.sqrt(m_xx*m_yy)*m_0/(h_bar**2)*k_B*Temp

del_z = 0.1e-9
N_grid_z = int(L_z/del_z+1)
print("Total grid in the width = ",N_grid_z)
z_grid = np.linspace(0, L_z, N_grid_z)

quantum_n =  np.arange(1,21)
E_z = (h_bar**2)/(2*m_zz*m_0)*(np.pi*quantum_n/L_z)**2
E_fermi = np.linspace(0.1,-0.1,11)

n_integrated = np.array([])
for Ef in E_fermi:
	accupancy = coeff*np.log(1+np.exp(-(E_z-Ef*1.60217646e-19)/(k_B*Temp)))
	N_total = np.sum(accupancy)
	density_band = 2/(L_x*L_y*L_z)*(np.sin(np.outer(quantum_n,z_grid)*np.pi/L_z))**2
	density_z = np.matmul(accupancy , density_band)
	n_integrated = np.append(n_integrated,np.sum(density_z)*del_z)
	subplot_n_z.semilogy(z_grid*1e9, density_z*1.0e-6, color=(0.7,(Ef*3+0.4),(Ef*2+0.2)), label = "%.2f eV" %Ef)

subplot_n_z.set_xlabel("position (nm)",fontsize=11)
subplot_n_z.set_ylabel('$\mathrm{n_{e}}$ ($\mathrm{cm^{-3}}$)' ,fontsize=11)
subplot_n_z.grid(True)
subplot_n_z.set_xlim(0.0,5)
subplot_n_z.set_ylim(1e14,1e20)
subplot_n_z.legend(bbox_to_anchor=(1.0, 1), loc=2,fontsize=5)
subplot_n_z.tick_params(which='both',direction='in')

subplot_n_Ef.semilogy(E_fermi, n_integrated*1.0e-4, 'ko')
subplot_n_Ef.set_xlabel("$\mathrm{E_{fermi}} (eV)$",fontsize=11)
subplot_n_Ef.set_ylabel('$\mathrm{n_{integrated}}$ ($\mathrm{cm^{-2}}$)' ,fontsize=11)
subplot_n_Ef.grid(True)
subplot_n_Ef.tick_params(which='both',direction='in')

plt.suptitle("Electron density in a box",fontsize=20)
plt.savefig("./density_inBox_result.png")
plt.show()
