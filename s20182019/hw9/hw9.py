import numpy as np
from numpy.linalg import inv
import copy


ny = 9
nz = 5

b_init = np.zeros(shape=(ny*nz),dtype=np.float64)

A = np.zeros(shape=(ny*nz,ny*nz),dtype=np.float64)


### Bulk terms
for i in range(1,ny-1):
    for j in range(1,nz-1):
        ind1 = i + ny*j
        A[ind1][ind1-ny] = 1.0
        A[ind1][ind1-1] = 1.0
        A[ind1][ind1] = -4.0
        A[ind1][ind1+1] = 1.0
        A[ind1][ind1+ny] = 1.0

### Bottom boundary
#for i in range(1,ny-1):
#    A[i][i-1] = 0.5
#    A[i][i] = -2.0
#    A[i][i+1] = 0.5
#    A[i][i+ny] = 1.0

### Left boundary
for i in range(1,nz-1):
    ind1 = i*ny
    A[ind1][ind1-ny] = 0.5
    A[ind1][ind1] = -2.0
    A[ind1][ind1+1] = 1.0
    A[ind1][ind1+ny] = 0.5

### Right boundary
for i in range(1,nz-1):
    ind1 = (i+1)*ny -1
    A[ind1][ind1-ny] = 0.5
    A[ind1][ind1-1] = 1.0
    A[ind1][ind1] = -2.0
    A[ind1][ind1+ny] = 0.5

### Top boundary
for i in range(2,ny-2):
    ind1 = (nz-1)*ny+i
    A[ind1][ind1-ny] = 1.0
    A[ind1][ind1-1] = 0.5
    A[ind1][ind1] = -2.0
    A[ind1][ind1+1] = 0.5

A[(nz-1)*ny][(nz-1)*ny] = 1.0
A[(nz-1)*ny+1][(nz-1)*ny+1] = 1.0
A[nz*ny-1][nz*ny-1] = 1.0
A[nz*ny-2][nz*ny-2] = 1.0

for i in range(ny):
    A[i][i] = 1.0

#A[0][0] = -1.0
#A[0][1] = 0.5
#A[0][ny] = 0.5

#A[ny-1][ny-1] = -1.0
#A[ny-1][ny-2] = 0.5
#A[ny-1][2*ny-1] = 0.5

#A[(nz-1)*ny][(nz-1)*ny] = -1.0
#A[(nz-1)*ny][(nz-1)*ny+1] = 0.5
#A[(nz-1)*ny][(nz-2)*ny] = 0.5

#A[nz*ny-1][nz*ny-1] = -1.0
#A[nz*ny-1][nz*ny-2] = 0.5
#A[nz*ny-1][nz*ny-1-ny] = 0.5



d = 1 # Delta, distance between one node and the nearest node

##### First : Original position #####

b = copy.deepcopy(b_init)

b[ny*nz-1] = 1
b[ny*nz-2] = 1
A_h = A[:]


psi1 = np.dot(inv(A_h),b)

f = open("psi1.dat",'w')
for i in range(nz):
    for j in range(ny):
        f.write("{}\t{}\t{}\n".format(j,i,psi1[j+ny*i]))

f.close()

##### Second : psi = 1 at left top

b = copy.deepcopy(b_init)

b[ny*(nz-1)] = 1
b[ny*(nz-1)+1] = 1
print(b)

psi2 = np.dot(inv(A),b)

f = open("psi2.dat",'w')
for i in range(nz):
    for j in range(ny):
        f.write("{}\t{}\t{}\n".format(j,i,psi2[j+ny*i]))

f.close()

##### Third : psi = 1 at bottom

b = copy.deepcopy(b_init)

for i in range(ny):
    b[i] = 1

psi3 = np.dot(inv(A),b)

f = open("psi3.dat",'w')
for i in range(nz):
    for j in range(ny):
        f.write("{}\t{}\t{}\n".format(j,i,psi3[j+ny*i]))

f.close()

##### 4th : psi = 1 at all boundaries above

b = copy.deepcopy(b_init)

b[ny*(nz-1)] = 1
b[ny*(nz-1)+1] = 1

b[ny*(nz)-1] = 1
b[ny*(nz)-2] = 1

for i in range(ny):
    b[i] = 1

psi4 = np.dot(inv(A),b)

f = open("psi4.dat",'w')
for i in range(nz):
    for j in range(ny):
        f.write("{}\t{}\t{}\n".format(j,i,psi4[j+ny*i]))

f.close()
