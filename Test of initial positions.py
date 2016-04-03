import scipy.stats as stats
import numpy as np
from numpy import sqrt
from decimal import Decimal as D
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

N_Particles = input("please, Enter number of the particles=")

Box_Length = 10
N_Grids = math.ceil(N_Particles**(.333334))

#initial condition
Grid = (Box_Length/N_Grids) 
x_Array=[]
y_Array=[]
z_Array=[]

print N_Grids
print Grid
#print math.floor (10.95) # round down
#print math.ceil (10.95)	#round up

for i in range (0,N_Particles):
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


	
#for i in range(0,(N_Particles)):
#	print i
#	for j in range(0,(N_Particles)):
#		print j
#		for k in range (0,(N_Particles)):
#			print k
#			if len(x_Array) < N_Particles:
#				if Grid*j < Box_Length and Grid*i <Box_Length and Grid*k <Box_Length:
#					x_Array.append(Grid * k + (Grid/3))
#					y_Array.append(Grid * j + (Grid/3))
#					z_Array.append(Grid * i + (Grid/3))
					

print len(x_Array), len(y_Array), len(z_Array)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


ax.scatter(x_Array, y_Array, z_Array, c='r', marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
