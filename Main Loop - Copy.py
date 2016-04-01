import scipy.stats as stats
import numpy as np
from numpy import sqrt
from decimal import Decimal as D
import random


Particles = input("please, Enter number of the particles=")
timestep = input("please, Enter number of time steps=")


# Constants
epsilon = float(120*1.38*10**(-23))
sigma = float(3.4*10**(-10))
Mass_Argon = float(6.63352*10**(-26))
K_B = float(1.38*10**(-23))


# simulation assumptions
Box_Length = float(float(10.229)*sigma)
timeint = float(0.2*10**(-14))
timemax = float(timestep*timeint)
N_Particles = Particles
targetT = (90)

#initial condition
Grid = (Box_Length/sqrt(N_Particles))
x_Array=[]
y_Array=[]
z_Array=[]

vx_Array=[]
vy_Array=[]
vz_Array=[]

forceX01 = 0
forceY01 = 0
forceZ01 = 0
	
KE = 0
PE = 0

for i in range(0,(N_Particles)):
	for j in range(0,(N_Particles)):
		if len(x_Array) < N_Particles:
			if Grid*j < Box_Length and Grid*i <Box_Length:
				x_Array.append(Grid * j)
				y_Array.append(Grid * i)
				z_Array.append(Grid * i)


for i in range(0,N_Particles):
	vx_Array.append (float(sqrt(K_B*targetT/Mass_Argon)))
	vy_Array.append (float(sqrt(K_B*targetT/Mass_Argon)))
	vz_Array.append (float(sqrt(K_B*targetT/Mass_Argon)))
	vx_Array [i] = float(vx_Array[i])
	vy_Array [i] = float(vy_Array[i])
	vz_Array [i] = float(vy_Array[i])

result = open ('output01.txt', 'w')
output = open ('output02.txt', 'w')

for t in np.arange (0,timemax,timeint):
	for j in range ( 0 , N_Particles):
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
		
		vx_Array [j] = vx_Array [j] + float(0.5) * timeint * (forceX01/Mass_Argon) 
		vy_Array [j] = vy_Array [j] + float(0.5) * timeint * (forceY01/Mass_Argon) 
		vz_Array [j] = vz_Array [j] + float(0.5) * timeint * (forceY01/Mass_Argon) 		

		#Calculation of KE
		
		
		KE = KE + (float(0.5) *Mass_Argon* vx_Array [j] **float(2)) + (float(0.5) *Mass_Argon* vy_Array [j] **float(2)) + (float(0.5) *Mass_Argon* vz_Array [j] **float(2))
		
		sysT = 2*KE/(3*N_Particles*K_B)
		
		if sysT > targetT +2 or sysT < targetT -2:
			vx_Array [j] = float(random.normalvariate(float(((targetT/sysT)**float(0.5)*vx_Array[j])-(sum(vx_Array)/len(vx_Array))),float(0.5)))
			vy_Array [j] = float(random.normalvariate(float(((targetT/sysT)**float(0.5)*vy_Array[j])-(sum(vy_Array)/len(vy_Array))),float (0.5)))
			vz_Array [j] = float(random.normalvariate(float(((targetT/sysT)**float(0.5)*vz_Array[j])-(sum(vz_Array)/len(vz_Array))),float (0.5)))				


	# Calculation of position

		x_Array [j] = x_Array [j] + timeint * vx_Array [j] + float(0.5)  * timeint**float(2) * (forceX01 /Mass_Argon) 
		y_Array [j] = y_Array [j] + timeint * vy_Array [j] + float(0.5)  * timeint**float(2) * (forceY01 /Mass_Argon)
		z_Array [j] = z_Array [j] + timeint * vz_Array [j] + float(0.5)  * timeint**float(2) * (forceY01 /Mass_Argon)

	# Boundary condition	
		if x_Array [j] > Box_Length or x_Array [j] < 0:
			x_Array [j] = float(random.random())*Box_Length
		
		if y_Array [j] > Box_Length or y_Array [j] < 0:
			y_Array [j] = float(random.random())*Box_Length

		if z_Array [j] > Box_Length or z_Array [j] < 0:
			z_Array [j] = float(random.random())*Box_Length
	
		
		#Calculation of PE
		
		if j != (N_Particles-1):
			for i in range (j+1 , N_Particles):
					
					delx = float(x_Array[i]-x_Array[j])
					dely = float(y_Array[i]-y_Array[j])
					delz = float(z_Array[i]-z_Array[j])
					r = float(sqrt((delx)**2+(dely)**2+(delz)**2))
					PE  = PE + (4 * epsilon * (((sigma/(r))**float(12))-((sigma/(r))**float(6))))
					
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
	
	E = PE + KE
	output.writelines (str('%.5e'%t))
	output.write("	")
	output.writelines (str('%.5e'%KE))
	output.write("	")
	output.writelines (str('%.5e'%PE))
	output.write("	")
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