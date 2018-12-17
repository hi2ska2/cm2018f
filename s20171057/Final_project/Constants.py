import numpy as np

hight_Oxide = 0.5			# nm
hight_Si = 5					# nm
width_1 = 40				# nm
width_2 = 40				# nm
width_3 = 40				# nm
del_width = 1.0			# nm
del_hight = 0.1			# nm
epsilon_ox = 3.9
epsilon_si = 11.7
N_don_1 = 5e25
N_don_2 = 2E23				# m^-3
N_don_3 = 5e25

err = 1e-20
total_width = width_1 + width_2 + width_3
total_hight = 2*hight_Oxide + hight_Si

N_Oxi_z = int(hight_Oxide/del_hight + err)
N_Si_z = int(hight_Si/del_hight + err)
N_Si_1_y = int(width_1/del_width + err)
N_Si_2_y = int(width_2/del_width + err)
N_Si_3_y = int(width_3/del_width + err)

N_y = N_Si_1_y + N_Si_2_y + N_Si_3_y +1
N_z = N_Oxi_z + N_Si_z + N_Oxi_z +1

print("Total size of the device = ",total_width,'x',total_hight,"nm")
print("Total grids in the size = ",N_y,'x',N_z)

n_intrin = 1.075e16
q_electron = 1.602192E-19		#C
epsilon_0 = 8.854187817E-12		#F/m
k_B = 1.380662E-23				#J/K
Temp = 300						#K
Volt_thermal = k_B*Temp/q_electron
coeff_rate = (del_hight/del_width)**2		# y축 방향 epsilon의 rate**2로 나눠줘야함
del_hight = del_hight*1e-9
del_width = del_width*1e-9
coeff = (del_hight)**2*q_electron/epsilon_0

def MakeVector(N_grid_1,N_grid_2,N_grid_3,v_1,v_2,v_3):
	vector = np.block([
			v_1*np.ones(N_grid_1),
			np.array([(v_1+v_2)/2]),
			v_2*np.ones(N_grid_2-1),
			np.array([(v_2+v_3)/2]),
			v_3*np.ones(N_grid_3)
			])
	return vector
N_don_yvec = MakeVector(N_Si_1_y,N_Si_2_y,N_Si_3_y,N_don_1,N_don_2,N_don_3)
epsilon_zvec = MakeVector(N_Oxi_z,N_Si_z,N_Oxi_z,epsilon_ox,epsilon_si,epsilon_ox)

_jaco = np.zeros((N_y*N_z,N_y,N_z))

i = 0
for z in range(N_z):
	for y in range(N_y):
		if(y==0 or y==N_y-1 or z==0 or z==N_z-1):	#Boundary condition
			if((z==0 or z==N_z-1) and (y>=N_Si_1_y and y<N_y-N_Si_3_y)):	# Gate
				_jaco[i,y,z]= 1
			elif(y==0 and (z>=N_Oxi_z and z<N_z-N_Oxi_z)):	# Source
				_jaco[i,y,z]= 1
			elif(y==N_y-1 and (z>=N_Oxi_z and z<N_z-N_Oxi_z)):	# Drain
				_jaco[i,y,z]= 1
			elif(z==0):	# Nuemann BC in Oxide_1
				if(y==0):
					_jaco[i,y,z] = -0.5*epsilon_ox -0.5*epsilon_ox*coeff_rate
					_jaco[i,y+1,z] = 0.5*epsilon_ox*coeff_rate
					_jaco[i,y,z+1] = 0.5*epsilon_ox
				elif(y==N_y-1):
					_jaco[i,y,z] = -0.5*epsilon_ox -0.5*epsilon_ox*coeff_rate
					_jaco[i,y-1,z] = 0.5*epsilon_ox*coeff_rate
					_jaco[i,y,z+1] = 0.5*epsilon_ox
				else:
					_jaco[i,y,z+1],_jaco[i,y-1,z],_jaco[i,y+1,z]= epsilon_ox, 0.5*epsilon_ox*coeff_rate, 0.5*epsilon_ox*coeff_rate
					_jaco[i,y,z]= -1*epsilon_ox -1*epsilon_ox*coeff_rate
			elif(z==N_z-1):	# Nuemann BC in Oxide_2
				if(y==0):
					_jaco[i,y,z] = -0.5*epsilon_ox -0.5*epsilon_ox*coeff_rate
					_jaco[i,y+1,z] = 0.5*epsilon_ox*coeff_rate
					_jaco[i,y,z-1] = 0.5*epsilon_ox
				elif(y==N_y-1):
					_jaco[i,y,z] = -0.5*epsilon_ox -0.5*epsilon_ox*coeff_rate
					_jaco[i,y-1,z] = 0.5*epsilon_ox*coeff_rate
					_jaco[i,y,z-1] = 0.5*epsilon_ox
				else:
					_jaco[i,y,z-1],_jaco[i,y-1,z],_jaco[i,y+1,z]= epsilon_ox, 0.5*epsilon_ox*coeff_rate, 0.5*epsilon_ox*coeff_rate
					_jaco[i,y,z]=-1*epsilon_ox -1*epsilon_ox*coeff_rate
			elif(y==0):
				_jaco[i,y,z+1],_jaco[i,y,z-1],_jaco[i,y+1,z]= 0.5*epsilon_ox, 0.5*epsilon_ox, epsilon_ox*coeff_rate
				_jaco[i,y,z]=-1*epsilon_ox -1*epsilon_ox*coeff_rate
			elif(y==N_y-1):
				_jaco[i,y,z+1],_jaco[i,y,z-1],_jaco[i,y-1,z]= 0.5*epsilon_ox, 0.5*epsilon_ox, epsilon_ox*coeff_rate
				_jaco[i,y,z]=-1*epsilon_ox -1*epsilon_ox*coeff_rate
			else:
				print("error")
				exit()
		elif(z==N_Oxi_z-1):
			_jaco[i,y,z-1],_jaco[i,y,z+1],_jaco[i,y-1,z],_jaco[i,y+1,z]= epsilon_ox, epsilon_ox, epsilon_ox*coeff_rate, epsilon_ox*coeff_rate
			_jaco[i,y,z]= -2*epsilon_ox -2*epsilon_ox*coeff_rate
		elif(z==N_Oxi_z+1):
			_jaco[i,y,z-1],_jaco[i,y,z+1],_jaco[i,y-1,z],_jaco[i,y+1,z]= epsilon_si, epsilon_si, epsilon_si*coeff_rate, epsilon_si*coeff_rate
			_jaco[i,y,z]= -2*epsilon_si -2*epsilon_si*coeff_rate
		elif(z==N_z-N_Oxi_z-2):
			_jaco[i,y,z-1],_jaco[i,y,z+1],_jaco[i,y-1,z],_jaco[i,y+1,z]= epsilon_si, epsilon_si, epsilon_si*coeff_rate, epsilon_si*coeff_rate
			_jaco[i,y,z]= -2*epsilon_si -2*epsilon_si*coeff_rate
		elif(z==N_z-N_Oxi_z):
			_jaco[i,y,z-1],_jaco[i,y,z+1],_jaco[i,y-1,z],_jaco[i,y+1,z]= epsilon_ox, epsilon_ox, epsilon_ox*coeff_rate, epsilon_ox*coeff_rate
			_jaco[i,y,z]= -2*epsilon_ox -2*epsilon_ox*coeff_rate
		else:
			_jaco[i,y,z-1],_jaco[i,y,z+1],_jaco[i,y-1,z],_jaco[i,y+1,z]= epsilon_zvec[z-1], epsilon_zvec[z+1], epsilon_zvec[z]*coeff_rate, epsilon_zvec[z]*coeff_rate
			_jaco[i,y,z]= -2*epsilon_zvec[z] -2*epsilon_zvec[z]*coeff_rate
		i=i+1

np.save('./jaco.npy', _jaco)