######################################
Reg = 2e6			# Ohm
Cap = 5e-12			# F
######################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

timestep = 101
f_grid = np.linspace(1.0,1000.0,101)
T_grid = np.linspace(0,1,timestep)
results_Ivin = np.zeros((np.size(f_grid),timestep))
results_IC = np.zeros((np.size(f_grid),timestep))
results_Ir = np.zeros((np.size(f_grid),timestep))
results_Vin = np.zeros((np.size(f_grid),timestep))
results_Vout = np.zeros((np.size(f_grid),timestep))
i=0
for freq in f_grid:
	result = (0,0,0,1,0)	# initial condition
	time_grid = T_grid/freq
	del_t = time_grid[1]-time_grid[0]
	Matrix = ((0,0,0,1,0),
			(0,1,0,-Cap/del_t,Cap/del_t),
			(0,0,1,0,-1/Reg),
			(1,1,0,0,0),
			(0,-1,1,0,0))
	Vec = np.zeros(5)
	j=0
	for time in time_grid:
		Vec[0] = np.cos(2*np.pi*freq*time)
		Vec[1] = -Cap/del_t*(result[3]-result[4])
		result = np.linalg.solve(Matrix,Vec)
		results_Ivin[i,j] = result[0]
		results_IC[i,j] = result[1]
		results_Ir[i,j] = result[2]
		results_Vin[i,j] = result[3]
		results_Vout[i,j] = result[4]
		j +=1
	i+=1

f_grid, T_grid = np.meshgrid(f_grid, T_grid)
## Analytic Solution 
omega = lambda f : 2*np.pi*f
I_exact = lambda f,T : omega(f)**2 *Reg*Cap**2/(1+(omega(f)*Reg*Cap)**2)*np.cos(2*np.pi*T) - omega(f)*Cap/(1+(omega(f)*Reg*Cap)**2)*np.sin(2*np.pi*T)

I_results = np.array(list(map(I_exact,f_grid,T_grid)))

## PLOT #######
fig = plt.figure()
fontsize = 12

ax = fig.gca(projection='3d')

sulf=ax.plot_wireframe(f_grid, T_grid , results_Ir.T, linewidth=2, rstride=10, cstride=10)
sulf=ax.plot_surface(f_grid, T_grid , I_results , cmap=cm.plasma, linewidth=2, antialiased=False, alpha=0.7)
ax.set_xlabel("frequency (Hz)",fontsize=fontsize)
ax.set_ylabel("period",fontsize=fontsize)
ax.set_zlabel("Current",fontsize=fontsize)

plt.show()
plt.savefig("./result_P2.png")
