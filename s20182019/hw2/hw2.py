import numpy as np
import numpy.linalg as lin

N_given = [3, 48, 498]
a = 5.0 # 5 nm
m = 0.19*9.10938356e-31 # mass
hbar = 1.0545718e-34 # hbar
E = 0.0
dx = 0.0

E_exact = 0.5*hbar*hbar*(np.pi*np.pi/a/a)/m
E_exact *= 1.0e+18*6.242e+18
print("Exact E_g = {} eV".format(E_exact))

for N in N_given:
    fname = 'data_%d.dat'%(N,)
    f = open(fname,'w')
    A = np.zeros(shape=(N,N))
    dx = a/(N+1)
    A[0][0] = -2.0
    A[0][1] = 1.0
    A[N-1][N-1] = -2.0
    A[N-1][N-2] = 1.0
    for i in range(N-2):
        A[i+1][i] = 1.0
        A[i+1][i+1] = -2.0
        A[i+1][i+2] = 1.0
    D, V = lin.eig(A)   # D: eig. value, V: eig. vector
    Ds = np.sort(D)
    Vs = V[:, D.argsort()]
    Vs = Vs.transpose()
    E = -hbar*hbar*0.5*Ds[N-1]/dx/dx/m
    E *= 1.0e+18 # nm -> m
    E *= 6.242e+18
    print("N={}, E = {} eV".format(N+2,E))
    vmax = np.amin(Vs[N-1])
    for i in range(N):
        f.write("{}\t{}\n".format((i+1)*dx,Vs[N-1][i]/vmax))
    f.close()
    print("Error = {}%\n".format((E_exact-E)/E_exact*100))
