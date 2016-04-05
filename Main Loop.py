import scipy.stats as stats
import numpy as np
from numpy import sqrt
from decimal import Decimal as D
import math
import time
start_time = time.time()



N_Particles = input("please, Enter number of the particles=")
N_time_step = input("please, Enter number of time steps=")


# Constants
epsilon = float(120*1.38*10**(-23))
sigma = float(3.4*10**(-10))
Mass_Argon = float(6.63352*10**(-26))
K_B = float(1.38*10**(-23))


# Box and Grids
Box_Length = float(float(10.229)*sigma)
N_Grids = math.ceil(N_Particles**(.333334))
Grid = (Box_Length/N_Grids) 

# simulation assumptions
time_step = float(0.2*10**(-14))
timemax = float(N_time_step*time_step)
targetT = (90)

#initial condition

x_Array=[]
y_Array=[]
z_Array=[]

vx_Array = []
vy_Array = []
vz_Array = []

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

for i in range (0,N_Particles):
	vx_Array.append(float(sqrt(K_B*targetT/Mass_Argon)))
	vy_Array.append(float(sqrt(K_B*targetT/Mass_Argon)))
	vz_Array.append(float(sqrt(K_B*targetT/Mass_Argon)))
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




for t in range (0,N_time_step):


	for j in range ( 0 , N_Particles):

		if t != 0 :
			for i in range (0 , N_Particles):
				if i!=j:	
					delx = float(x_Array[i]-x_Array[j])
					dely = float(y_Array[i]-y_Array[j])
					delz = float(z_Array[i]-z_Array[j])
					r = sqrt((delx)**2+(dely)**2+(delz)**2)
					
					forceX01  = forceX01  + ( float(4) * epsilon * (( float(-12) * (sigma**float(12)/(r)**float(13))) + float(6) *(sigma**float(6)/(r)**float(7))))*(delx/(r))
					forceY01  = forceY01  + ( float(4) * epsilon * (( float(-12) * (sigma**float(12)/(r)**float(13))) + float(6) *(sigma**float(6)/(r)**float(7))))*(dely/(r))
					forceZ01  = forceZ01  + ( float(4) * epsilon * (( float(-12) * (sigma**float(12)/(r)**float(13))) + float(6) *(sigma**float(6)/(r)**float(7))))*(delz/(r))		


					

			# calculation of Velocities
			
			vx_Array [j] = vx_Array [j] + float(0.5) * time_step * (forceX01/Mass_Argon) 
			vy_Array [j] = vy_Array [j] + float(0.5) * time_step * (forceY01/Mass_Argon) 
			vz_Array [j] = vz_Array [j] + float(0.5) * time_step * (forceY01/Mass_Argon) 		

		#Calculation of KE


		if t != 0:

			if sysT > targetT +2 or sysT < targetT -2:
				vx_Array [j] = float(((targetT/sysT)**float(0.5)*vx_Array[j])-(sum(vx_Array)/len(vx_Array)))
				vy_Array [j] = float(((targetT/sysT)**float(0.5)*vy_Array[j])-(sum(vy_Array)/len(vy_Array)))
				vz_Array [j] = float(((targetT/sysT)**float(0.5)*vz_Array[j])-(sum(vz_Array)/len(vz_Array)))				


	# Calculation of position
		if t != 0:
			x_Array [j] = x_Array [j] + time_step * vx_Array [j] + float(0.5)  * time_step**float(2) * (forceX01 /Mass_Argon) 
			y_Array [j] = y_Array [j] + time_step * vy_Array [j] + float(0.5)  * time_step**float(2) * (forceY01 /Mass_Argon)
			z_Array [j] = z_Array [j] + time_step * vz_Array [j] + float(0.5)  * time_step**float(2) * (forceY01 /Mass_Argon)

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

				
	
		
		
		print "time step =",t+1, "PARTICLE # ",j+1," out of", N_Particles,"\r",
	
	
	#Calculation of PE
	for j in range (0 , N_Particles):
		for i in range (j , N_Particles):
			if i !=j :	
				delx = float(x_Array[i]-x_Array[j])
				dely = float(y_Array[i]-y_Array[j])
				delz = float(z_Array[i]-z_Array[j])
				r = float(sqrt((delx)**2+(dely)**2+(delz)**2))
				PE  = PE + (4 * epsilon * (((sigma/(r))**float(12))-((sigma/(r))**float(6))))
		KE = KE	+ 0.5 * Mass_Argon * (vx_Array[j]**2 + vy_Array[j]**2+vz_Array[j]**2)
		sysT = 2*KE/(3*N_Particles*K_B)
		print "time step =",t+1, "PARTICLE # ",j+1," out of", N_Particles,"\r",
		

		result.writelines (str('%.5e'%t))
		result.write("	")
		result.writelines (str('%.5e'%x_Array[j]))
		result.write("	")
		result.writelines (str('%.5e'%y_Array[j]))
		result.write("	")
		result.writelines (str('%.5e'%z_Array[j]))
		result.write("	")
		result.writelines (str('%.5e'%vx_Array[j]))
		result.write("	")
		result.writelines (str('%.5e'%vy_Array[j]))
		result.write("	")
		result.writelines (str('%.5e'%vz_Array[j]))
		
		result.write ("\n")
	
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

