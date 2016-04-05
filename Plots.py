import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from decimal import Decimal as D


file0 = open("initials.txt", 'r')
lines0 = file0.readlines()
file0.close()

file = open("output01.txt", 'r')
lines = file.readlines()
file.close()

file2 = open("output02.txt", 'r')
lines2 = file2.readlines()
file2.close() 

N_particles = float(lines0 [1])
N_time_step = float(lines0 [3])
Box_Length = float(lines0 [5])
epsilon = float(lines0 [7])
sigma = float(lines0 [9])
Mass_Argon = float(lines0 [11])
K_B = float(lines0 [13])
time_step = float(lines0 [15])
targetT = float(lines0 [17])

time =[]
KE = []
PE =[]
E = []
sysT =[]

xpos=[]
ypos=[]
zpos =[]
r = []

N_particles = int (N_particles)

for line in lines2:
    parts = line.split()
    if len(parts) > 1:
		time.append (float(parts [0]))
		KE.append (float(parts [1]))
		PE.append (float(parts [2]))
		E.append (float(parts [3]))
		sysT.append (float(parts [4]))
		
for line in lines:
    parts = line.split() 
    if len(parts) > 1:
		xpos.append (float(parts [1]))
		ypos.append (float(parts [2]))
		zpos.append (float(parts [3]))


	
####### Radial Densitysity ##########


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
		Density[n]=Mass_Argon*Density[n]/((3.1415)*((n+1)*Box_Length/N_Devisions)**(2))
	else:
		Density[n]=Mass_Argon*Density[n]/(((3.1415)*((n+1)*Box_Length/N_Devisions)**(2))-((3.1415)*n*(Box_Length/N_Devisions)**(2)))

####### end of Radial Densitysity ##########


		
plt.subplot (221)
plt.plot (sysT,'r--')
plt.title ('system temperature')

plt.subplot (222)
plt.plot (KE,'r--')
plt.title ('KE')

plt.subplot (223)
plt.plot (PE,'r--')
plt.title ('PE')

plt.subplot (224)
plt.plot (E,'r--')
plt.title ('E')

plt.show()