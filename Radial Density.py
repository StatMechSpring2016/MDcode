import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from decimal import Decimal as D


# Constants
epsilon = D(120*1.38*10**(-23))
sigma = D(3.4*10**(-10))
Mass_Argon = D(6.63352*10**(-26))
K_B = D(1.38*10**(-23))


# simulation assumptions
Box_Length = D(D(10.229)*sigma)
timeint = D(0.2*10**(-14))
targetT = (90)

file = open("output01.txt", 'r')
lines = file.readlines()
file.close()

file2 = open("output02.txt", 'r')
lines2 = file2.readlines()
file2.close() 



time =[]
KE = []
PE =[]
E = []
sysT =[]

xpos=[]
ypos=[]
zpos =[]
r = []

for line in lines2:
    parts = line.split()
    if len(parts) > 1:
		time.append (D(parts [0]))
		KE.append (D(parts [1]))
		PE.append (D(parts [2]))
		E.append (D(parts [3]))
		sysT.append (D(parts [4]))
		
for line in lines:
    parts = line.split() 
    if len(parts) > 1:
		xpos.append (D(parts [1]))
		ypos.append (D(parts [2]))
		zpos.append (D(parts [3]))

N_particles = len(xpos)/len(time)	

	
#Radial Densitysity


Density = []
N_Devisions = 10

for i in range (0,len(xpos)):
	r.append (sqrt(xpos[i]**2+ypos[i]**2+zpos[i]**2))		

for n in range (0,N_Devisions):
	Density.append(0)
	

for i in range (0,N_Devisions):
	for n in range(0,N_particles):
		if r[n] <= (i+1) * Box_Length/N_Devisions and r[n]> i*Box_Length/N_Devisions:
			Density[i] = Density[i] + 1

for n in range (0,N_Devisions):
	if n == 0:
		Density[n]=Mass_Argon*Density[n]/(D(3.1415)*((n+1)*Box_Length/N_Devisions)**D(2))
	else:
		Density[n]=Mass_Argon*Density[n]/((D(3.1415)*((n+1)*Box_Length/N_Devisions)**D(2))-(D(3.1415)*n*(Box_Length/N_Devisions)**D(2)))

		
		
plt.subplot (111)
plt.plot (Density,'r--')
plt.title ('Radial Density')

plt.show()