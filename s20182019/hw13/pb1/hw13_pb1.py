import numpy as np
from numpy.linalg import inv
import sys

q=1.602192e-19
eps0 = 8.854187817e-12
kB = 1.380662e-23
T = 300.0
VT = kB*T/q

N = int(sys.argv[1])+1
Dn = 0.01

Nd_high = 0
Nd_low = 0
L = 0
mid1 = 0
mid2 = 0

if int(sys.argv[2]) is 0: # Long
    Nd_high = 5.0e+23
    Nd_low = 2.0e+21
    L = 600 # 600 nm
    mid1 = 100
    mid2 = 500
elif int(sys.argv[2]) is 1: # Short
    Nd_high = 5.0e+25
    Nd_low = 2.0e+23
    L = 120 # 120 nm
    mid1 = 40
    mid2 = 80
else:
    print("Aborted: Choose one of them; 1 for long, 2 for short.")
    sys.exit()

st1 = []

if int(sys.argv[2]) is 0:
    st1 = "long"
elif int(sys.argv[2]) is 1:
    st1 = "short"


dx = (L/(N-1))

mid1 = int(mid1/dx)
mid2 = int(mid2/dx)

dx = dx * 1e-9

eps_si = 11.7
eps_ox = 3.9

Ni = 1.075e+16

Nd = np.ones(shape=(N))*Nd_low
Nd[0:mid1+1] = Nd_high
Nd[mid2:N] = Nd_high

coef = dx*dx*q/eps0

### Non-linear Poisson solver
phi = np.zeros(shape=(N),dtype=np.float64)
phi[:] = VT*np.log(Nd/Ni)

for iNewton in range(10):
    res = np.zeros(shape=(N),dtype=np.float64)
    J = np.zeros(shape=(N,N),dtype=np.float64)
    res[0] = phi[0] - VT*np.log(Nd[0]/Ni)
    J[0][0] = 1.0
    for i in range(1,N-1):
        res[i] = eps_si*(phi[i-1]-2.0*phi[i]+phi[i+1])
        J[i][i-1] = eps_si
        J[i][i] = -2.0*eps_si
        J[i][i+1] = eps_si
    res[N-1] = phi[N-1] - VT*np.log(Nd[N-1]/Ni)
    J[N-1][N-1] = 1.0
    for i in range(1,N-1):
        res[i] = res[i] - coef*(-Nd[i]+Ni*np.exp(phi[i]/VT))
        J[i][i] = J[i][i] - coef*Ni*np.exp(phi[i]/VT)/VT
    update = np.dot(inv(J),-res)
    phi = phi + update
    print(np.linalg.norm(update, np.inf))

#f = open("ED"+st1+"_"+str(N)+"_init.dat",'w')
#for i in range(N):
#    f.write("{}\t{}\n".format(i,phi[i]))
#f.close()

### Continuity equation of one step

elec = np.zeros(shape=(N),dtype=np.float64)
elec = Ni*np.exp(phi/VT)

res = np.zeros(shape=(N),dtype=np.float64)
J = np.zeros(shape=(N,N),dtype=np.float64)
for i in range(1,N-1):
    n_av = 0.5*(elec[i+1]+elec[i])
    dphidx = (phi[i+1]-phi[i])/dx
    delecdx = (elec[i+1]-elec[i])/dx
    Jn = n_av*dphidx - VT*delecdx
    res[i] = res[i] + Jn
    J[i][i] = J[i][i]+0.5*dphidx + VT/dx
    J[i][i+1] = J[i][i+1]+0.5*dphidx - VT/dx
    res[i+1] = res[i+1] - Jn
    J[i+1][i] = J[i+1][i] - 0.5*dphidx - VT/dx
    J[i+1][i+1] = J[i+1][i+1] - 0.5*dphidx - VT/dx

res[0] = elec[0] - Nd[0]
J[0][0] = 1.0
res[N-1] = elec[N-1] - Nd[N-1]
J[N-1][N-1] = 1.0

update = np.dot(inv(J),-res)
elec = elec + update


### Coupled equation

V_grid = np.linspace(0,1,101)
E_cur = np.zeros(shape=(101),dtype=np.float64)
ind1 = 0
print(V_grid)
for V in V_grid:
    J_c = np.zeros(shape=(N),dtype=np.float64)
    for iNewton in range(10):
        res = np.zeros(shape=(2*N),dtype=np.float64)
        J = np.zeros(shape=(2*N,2*N),dtype=np.float64)
        res[0] = phi[0] - VT*np.log(Nd[0]/Ni) - V
        J[0][0] = 1.0
        for i in range(1,N-1):
            res[2*i] = eps_si*(phi[i-1]-2.0*phi[i]+phi[i+1]) + coef*(Nd[i]-elec[i])
            J[2*i][2*i-2] = eps_si
            J[2*i][2*i] = -2.0*eps_si
            J[2*i][2*i+2] = eps_si
            J[2*i][2*i+1] = -coef
        res[2*N-2] = phi[N-1] - VT*np.log(Nd[N-1]/Ni)
        J[2*N-2][2*N-2] = 1.0

        for i in range(N-1):
            n_av = 0.5*(elec[i+1]+elec[i])
            dphidx = (phi[i+1]-phi[i])/dx
            delecdx = (elec[i+1]-elec[i])/dx
            Jn = n_av*dphidx - VT*delecdx
            res[2*i+1] = res[2*i+1] + Jn
            J[2*i+1][2*i+3] = J[2*i+1][2*i+3] + 0.5*dphidx - VT/dx
            J[2*i+1][2*i+1] = J[2*i+1][2*i+1] + 0.5*dphidx + VT/dx
            J[2*i+1][2*i+2] = J[2*i+1][2*i+2] + n_av/dx
            J[2*i+1][2*i] = J[2*i+1][2*i] - n_av/dx
            res[2*i+3] = res[2*i+3] - Jn
            J[2*i+3][2*i+3] = J[2*i+3][2*i+3] - 0.5*dphidx + VT/dx
            J[2*i+3][2*i+1] = J[2*i+3][2*i+1] - 0.5*dphidx - VT/dx
            J[2*i+3][2*i+2] = J[2*i+3][2*i+2] - n_av/dx
            J[2*i+3][2*i] = J[2*i+3][2*i] + n_av/dx
        res[1] = elec[0] - Nd[0]
        J[1][1] = 1.0
        res[2*N-1] = elec[N-1] - Nd[N-1]
        J[2*N-1][2*N-1] = 1.0

        update = np.dot(inv(J),-res)
        phi = phi + update[::2]
        elec = elec + update[1::2]
        print(np.linalg.norm(update[::2], np.inf))

    for i in range(N-1):
        k = (phi[i+1] - phi[i])/VT
        J_c[i] = (q*Dn/dx)*(elec[i+1]*k/(np.exp(k)-1) + elec[i]*k/(np.exp(-k)-1))
    E_cur[ind1] = -J_c[N-2]
    ind1 = ind1 + 1

f = open("current_"+st1+".dat",'w')

for i in range(101):
    f.write("{}\t{}\n".format(i*0.01,-E_cur[i]))
    


