import numpy as np


for i in range(11):
    phi = np.loadtxt("ed"+str(i)+".dat",usecols=(2,))
    f = open("ed"+str(i)+"_x.dat",'w')
    for i in range(121):
        f.write("{}\t{}\n".format(i,phi[121*30+i]))

    f.close()
