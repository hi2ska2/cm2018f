import numpy as np
from numpy.linalg import inv
import copy

R = 2.0e+6
C = 5.0e-12
fq = 5e+5    # Change the frequency!
dt = 1/fq/100 # 0.01 of a period, as in the example code
V0 = 1.0

A = np.zeros(shape=(5,5),dtype=np.float64)
A[0][:] = [0,0,0,1,0]
A[1][:] = [0,1,0,-C/dt,C/dt]
A[2][:] = [0,0,1,0,-1/R]
A[3][:] = [1,1,0,0,0]
A[4][:] = [0,-1,1,0,0]  # Variable order is {Ivin, Ic1, Ir1, Vin, Vout}.

b = np.zeros(shape=(5),dtype=np.float64)
sol = [0,0,0,1,0]   # Initial conditions.
N = 2000    # How many time points we record?

ind = int(fq)
f = open("sol_"+str(ind)+".dat",'w')
w = 2*np.pi*fq

for i in range(N):
    t = (i+1)*dt
    sol_old = copy.deepcopy(sol)
    b[0] = np.cos(w*t)
    b[1] = -C*(sol_old[3]-sol_old[4])/dt
    sol = np.dot(inv(A),b)
    f.write("{}\t{}\t{}\t{}\n".format(t,sol[0],sol[1],sol[2]))

f.close()

f1 = open("sol_"+str(ind)+"_exact.dat",'w')
for i in range(N):
    t = (i+1)*dt
    sole = (w*w*R*C*C)*V0*np.cos(w*t)/(1+(w*R*C)**2) -w*C*V0*np.sin(w*t)/(1+(w*R*C)**2)
    f1.write("{}\t{}\n".format(t,sole))

f1.close()
