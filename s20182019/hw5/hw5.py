import numpy as np
from numpy.linalg import inv

n_int = 1.0E+16 # Intrinsic carrier density (/m^3)
Vt = 1.0    # provided arbitrary.
N = 0.0

f = open("sol.dat",'w')
f.write("#N+\tphi_num\tphi_exact\tloop_number\n")

xi = 0.0

for i in range(18):
    if i < 9:
        N = 1.0E+16*(10**i)
        xi = 3.0+2.0*i
    else:
        N = -1.0E+16*(10**(i-9))
        xi = -(3.0+2.0*(i-9))
    xf = xi
    loop_index = 0
    dx = 0
    while(loop_index<50):
        loop_index += 1
        r = N - n_int*np.sinh(xi/Vt)
        dx = (-1.0)*r/(-1.0*n_int*np.cosh(xi/Vt))
#        print("x={}, dx={}, r={}".format(xi,dx,r))
        xf = xi + dx
        rf = N - n_int*np.sinh(xf/Vt)
        if (np.abs(r-rf)<1.0E-5) or (np.abs(dx)<1.0E-5):
            break
        else:
            xi = xf
    print("N+ = {}:\nAfter {} loops, x = {}, dx = {}".format(N,loop_index,xf,dx))
    x_real = Vt*np.arcsinh(N/n_int)
    print("x_exact={}".format(x_real))
    f.write("{}\t{}\t{}\t{}\n".format(N,xf,x_real,loop_index))

f.close()
