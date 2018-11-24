from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, MultipleLocator, ScalarFormatter
import numpy as np
import sys

fig = plt.figure()



Z = np.loadtxt("f"+sys.argv[1]+".dat",dtype=np.float32,usecols=(2,),unpack=True)

X = np.linspace(0,30,301)
Y = np.linspace(0.1,1.0,91)

X, Y = np.meshgrid(X,Y)

Z_grid = np.reshape(Z,(91,301))

print(X.shape)
print(Y.shape)
print(Z_grid.shape)


fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, Z_grid, cmap=cm.plasma, linewidth=0, antialiased=False, alpha=0.8)
#ax.set_zlim(0,1.0)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))
#ax.set_zscale("log")
ax.zaxis.set_major_formatter(ScalarFormatter())
ax.zaxis.set_minor_formatter(ScalarFormatter())
plt.autoscale(enable=True,axis='z')
ax.yaxis.set_major_formatter(FormatStrFormatter('%.01f'))
#ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_xlabel("x",fontsize=12)
ax.set_ylabel("H",fontsize=12)
ax.set_zlabel(r"$f_1$",fontsize=12)

fig.colorbar(surf, shrink=0.7,aspect=5)

plt.show()

fig.savefig("f"+sys.argv[1]+".pdf")
