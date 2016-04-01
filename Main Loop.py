import scipy.stats as stats
import numpy as np
from numpy import sqrt
from decimal import Decimal as D
import random


Particles = input("please, Enter number of the particles=")
timestep = input("please, Enter number of time steps=")


# Constants
epsilon = D(120*1.38*10**(-23))
sigma = D(3.4*10**(-10))
Mass_Argon = D(6.63352*10**(-26))
K_B = D(1.38*10**(-23))


# simulation assumptions
Box_Length = D(D(10.229)*sigma)
timeint = D(0.2*10**(-14))
timemax = D(timestep*timeint)
N_Particles = D(Particles)
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
	vx_Array [i] = D(vx_Array[i])
	vy_Array [i] = D(vy_Array[i])
	vz_Array [i] = D(vy_Array[i])

result = open ('output01.txt', 'w')
output = open ('output02.txt', 'w')

for t in np.arange (0,timemax,timeint):
	for j in range ( 0 , N_Particles):
		for i in range (0 , N_Particles):
			if i!=j:	
				delx = D(x_Array[i]-x_Array[j])
				dely = D(y_Array[i]-y_Array[j])
				delz = D(z_Array[i]-z_Array[j])
				r = sqrt((delx)**2+(dely)**2+(delz)**2)
				
				forceX01  = forceX01  + ( D(4) * epsilon * (( D(-12) * (sigma**D(12)/(r)**D(13))) + D(6) *(sigma**D(6)/(r)**D(7))))*(delx/(r))
				forceY01  = forceY01  + ( D(4) * epsilon * (( D(-12) * (sigma**D(12)/(r)**D(13))) + D(6) *(sigma**D(6)/(r)**D(7))))*(dely/(r))
				forceZ01  = forceZ01  + ( D(4) * epsilon * (( D(-12) * (sigma**D(12)/(r)**D(13))) + D(6) *(sigma**D(6)/(r)**D(7))))*(delz/(r))		


		# calculation of Velocities
		
		vx_Array [j] = vx_Array [j] + D(0.5) * timeint * (forceX01/Mass_Argon) 
		vy_Array [j] = vy_Array [j] + D(0.5) * timeint * (forceY01/Mass_Argon) 
		vz_Array [j] = vz_Array [j] + D(0.5) * timeint * (forceY01/Mass_Argon) 		

		#Calculation of KE
		
		
		KE = KE + (D(0.5) *Mass_Argon* vx_Array [j] **D(2)) + (D(0.5) *Mass_Argon* vy_Array [j] **D(2)) + (D(0.5) *Mass_Argon* vz_Array [j] **D(2))
		
		sysT = 2*KE/(3*N_Particles*K_B)
		
		if sysT > targetT +2 or sysT < targetT -2:
			vx_Array [j] = D(random.normalvariate(float(((targetT/sysT)**D(0.5)*vx_Array[j])-(sum(vx_Array)/len(vx_Array))),float(0.5)))
			vy_Array [j] = D(random.normalvariate(float(((targetT/sysT)**D(0.5)*vy_Array[j])-(sum(vy_Array)/len(vy_Array))),float (0.5)))
			vz_Array [j] = D(random.normalvariate(float(((targetT/sysT)**D(0.5)*vz_Array[j])-(sum(vz_Array)/len(vz_Array))),float (0.5)))				


	# Calculation of position

		x_Array [j] = x_Array [j] + timeint * vx_Array [j] + D(0.5)  * timeint**D(2) * (forceX01 /Mass_Argon) 
		y_Array [j] = y_Array [j] + timeint * vy_Array [j] + D(0.5)  * timeint**D(2) * (forceY01 /Mass_Argon)
		z_Array [j] = z_Array [j] + timeint * vz_Array [j] + D(0.5)  * timeint**D(2) * (forceY01 /Mass_Argon)

	# Boundary condition	
		if x_Array [j] > Box_Length or x_Array [j] < 0:
			x_Array [j] = D(random.random())*Box_Length
		
		if y_Array [j] > Box_Length or y_Array [j] < 0:
			y_Array [j] = D(random.random())*Box_Length

		if z_Array [j] > Box_Length or z_Array [j] < 0:
			z_Array [j] = D(random.random())*Box_Length
	
		
		#Calculation of PE
		
		if j != (N_Particles-1):
			for i in range (j+1 , N_Particles):
					
					delx = D(x_Array[i]-x_Array[j])
					dely = D(y_Array[i]-y_Array[j])
					delz = D(z_Array[i]-z_Array[j])
					r = D(sqrt((delx)**2+(dely)**2+(delz)**2))
					PE  = PE + (4 * epsilon * (((sigma/(r))**D(12))-((sigma/(r))**D(6))))
					
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