import numpy as np
from numpy import sqrt
import math
import time

# Inputs
N_Particles = input("please, Enter number of the particles=")
N_time_step = input("please, Enter number of time steps=")

start_time = time.time()

# Constants
epsilon = float(120*1.38*10**(-23))
sigma = float(3.4*10**(-10))
Mass_Argon = float(6.63352*10**(-26))
K_B = float(1.38*10**(-23))


# Box and Grids
Box_Length = float(10.229*sigma)
N_Grids = math.ceil(N_Particles**(.333334))
Grid = (Box_Length/N_Grids) 

# simulation assumptions
time_step = float(0.2*10**(-14))
timemax = N_time_step*time_step
targetT = float(90)

#initial conditions

x_Array=[]
y_Array=[]
z_Array=[]

vx_Array = np.array([])
vy_Array = np.array([])
vz_Array = np.array([])

forceX01 = 0
forceY01 = 0
forceZ01 = 0
	
KE = 0
PE = 0


#saving constants and assumptions
initials = open ('initials.txt','w')

initials.writelines ("N_Particles")
initials.write("\n")
initials.writelines (str('%.5e'%N_Particles))
initials.write("\n")

initials.writelines ("N_time_step")
initials.write("\n")
initials.writelines (str('%.5e'%N_time_step))
initials.write("\n")

initials.writelines ("Box_Length")
initials.write("\n")
initials.writelines (str('%.5e'%Box_Length))
initials.write("\n")

initials.writelines ("epsilon")
initials.write("\n")
initials.writelines (str('%.5e'%epsilon))
initials.write("\n")

initials.writelines ("sigma")
initials.write("\n")
initials.writelines (str('%.5e'%sigma))
initials.write("\n")

initials.writelines ("Mass_Argon")
initials.write("\n")
initials.writelines (str('%.5e'%Mass_Argon))
initials.write("\n")

initials.writelines ("K_B")
initials.write("\n")
initials.writelines (str('%.5e'%K_B))
initials.write("\n")

initials.writelines ("time_step")
initials.write("\n")
initials.writelines (str('%.5e'%time_step))
initials.write("\n")

initials.writelines ("targetT")
initials.write("\n")
initials.writelines (str('%.5e'%targetT))
initials.write("\n")
initials.close()

#setting initial values

v = sqrt(K_B*targetT/Mass_Argon)

for i in range (0,N_Particles):
	vx_Array = np.append(vx_Array,v*pow(-1,i))
	vy_Array = np.append(vy_Array,v*pow(-1,i))
	vz_Array = np.append(vz_Array,v*pow(-1,i))
	if Grid * i + (Grid/2) > Box_Length:
		i = i - (N_Grids)*math.floor(i/N_Grids)
	x_Array.append(Grid * i + (Grid/2))

	
for i in range (0,N_Particles):
	if i == N_Grids*math.floor(i/N_Grids)-1:
		i = i - (N_Grids)*math.floor(i/N_Grids)
	i = math.floor(i/N_Grids)
	if Grid * i + (Grid/2) > Box_Length:
		i = i - (N_Grids)*math.floor(i/N_Grids)
	y_Array.append(Grid * i + (Grid/2))

	
for i in range (0,N_Particles):
	if i == N_Grids*N_Grids*math.floor(i/N_Grids)-1:
		i = i - (N_Grids*N_Grids)*math.floor(i/(N_Grids*N_Grids))
	i = math.floor(i/(N_Grids*N_Grids))
	if Grid * i + (Grid/2) > Box_Length:
		i = i - (N_Grids*N_Grids)*math.floor(i/(N_Grids*N_Grids))
	z_Array.append(Grid * i + (Grid/2))	



	
		
result = open ('output01.txt', 'w')
output = open ('output02.txt', 'w')


# starting Main Loop

for t in range (0,N_time_step):

	# Calculation of force 
	for j in range ( 0 , N_Particles):

		if t != 0 :
			for i in range (0 , N_Particles):
				if i!=j:	
					delx = float(x_Array[i]-x_Array[j])
					dely = float(y_Array[i]-y_Array[j])
					delz = float(z_Array[i]-z_Array[j])
					r = sqrt((delx)**2+(dely)**2+(delz)**2)
					
					forceX01  = forceX01  + ( 4. * epsilon * ( -12 * pow(sigma,12.)/pow (r,13.) + 6. * pow(sigma,6.)/pow(r,7.))*(delx/r))
					forceY01  = forceY01  + ( 4. * epsilon * ( -12 * pow(sigma,12.)/pow (r,13.) + 6. * pow(sigma,6.)/pow(r,7.))*(dely/r))
					forceZ01  = forceZ01  + ( 4. * epsilon * ( -12 * pow(sigma,12.)/pow (r,13.) + 6. * pow(sigma,6.)/pow(r,7.))*(delz/r))		


					

	# calculation of Velocities
			
			vx_Array [j] = vx_Array [j] + 0.5 * time_step * (forceX01/Mass_Argon) 
			vy_Array [j] = vy_Array [j] + 0.5 * time_step * (forceY01/Mass_Argon) 
			vz_Array [j] = vz_Array [j] + 0.5 * time_step * (forceY01/Mass_Argon) 		

		

		
	# Calculation of position
		if t != 0:
			x_Array [j] = x_Array [j] + time_step * vx_Array [j] + 0.5  * time_step**(2) * (forceX01 /Mass_Argon) 
			y_Array [j] = y_Array [j] + time_step * vy_Array [j] + 0.5  * time_step**(2) * (forceY01 /Mass_Argon)
			z_Array [j] = z_Array [j] + time_step * vz_Array [j] + 0.5  * time_step**(2) * (forceY01 /Mass_Argon)

			
			
	# Boundary condition	
		while x_Array [j] > Box_Length:
			x_Array [j] = x_Array[j]-Box_Length
		while x_Array [j] < 0:
			x_Array [j] = x_Array[j]+Box_Length

		while y_Array [j] > Box_Length:
			y_Array [j] = y_Array[j]-Box_Length
		while y_Array [j] < 0:
			y_Array [j] = y_Array[j]+Box_Length

		while z_Array [j] > Box_Length:
			z_Array [j] = z_Array[j]-Box_Length
		while z_Array [j] < 0:
			z_Array [j] = z_Array[j]+Box_Length
		
		result.writelines (str('%.5e'%t))		
		result.write("	")
		result.writelines (str('%.5e'%x_Array[j]))
		result.write("	")
		result.writelines (str('%.5e'%y_Array[j]))
		result.write("	")
		result.writelines (str('%.5e'%z_Array[j]))
		result.write ("\n")
				
	
	
	
	# Calculation of PE and KE
	for j in range (0 , N_Particles):
		for i in range (j , N_Particles):
			if i !=j :	
				delx = float(x_Array[i]-x_Array[j])
				dely = float(y_Array[i]-y_Array[j])
				delz = float(z_Array[i]-z_Array[j])
				r = float(sqrt((delx)**2+(dely)**2+(delz)**2))
				PE  = PE + (4. * epsilon * (pow(sigma/r,12)-pow(sigma/r,6)))
				
				

	
	vx_Array = vx_Array-(np.sum(vx_Array)/len(vx_Array))
	vy_Array = vy_Array-(np.sum(vy_Array)/len(vy_Array))	
	vz_Array = vz_Array-(np.sum(vz_Array)/len(vz_Array))	

	Vsquared =  (vx_Array**2 + vy_Array**2 + vz_Array**2)
	KE = .5*Mass_Argon*np.sum(Vsquared)
		
		
	# Calculation of system Temprature
	sysT = 2*KE/(3*N_Particles*K_B)
	
	# Checking the target T and updating velocities and KE
	if t != 0:
		if sysT > targetT +0.1 or sysT < targetT -0.1:
			vx_Array = vx_Array*sqrt(targetT/sysT)
			vy_Array = vy_Array*sqrt(targetT/sysT)
			vz_Array = vz_Array*sqrt(targetT/sysT)

			
			Vsquared =  vx_Array**2 + vy_Array**2 + vz_Array**2
			KE = .5 * Mass_Argon * np.sum(Vsquared)
			sysT = 2*KE/(3*N_Particles*K_B)
		
	print "time step =",t+1,"\r",

		

	
	output.writelines (str('%.5e'%t))
	output.write("	")
	
	KE = (KE /1000)/(N_Particles/(6.02*10**23))
	output.writelines (str('%.5e'%KE))
	output.write("	")
	
	PE = (PE /1000)/(N_Particles/(6.02*10**23))	
	output.writelines (str('%.5e'%PE))
	output.write("	")
	
	E = PE + KE
	output.writelines (str('%.5e'%E))
	output.write("	")
	
	output.writelines (str('%.5e'%sysT))
	output.write("	")
	output.write ("\n")
	forceX01 = 0
	forceY01 = 0
	PE = 0
	KE = 0

result.close()
output.close()
print "\n","finished"
# your code
elapsed_time = time.time() - start_time
print "elapsed time:",math.ceil(elapsed_time) , "seconds"
print "elapsed time:",math.ceil(elapsed_time)/60 , "minutes"
print "elapsed time:",math.ceil(elapsed_time)/3600 , "hours"

