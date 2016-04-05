import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from decimal import Decimal as D
from mpl_toolkits.mplot3d import Axes3D


file0 = open("initials.txt", 'r')
lines0 = file0.readlines()
file0.close()

file = open("output01.txt", 'r')
lines = file.readlines()
file.close()



N_particles = float(lines0 [1])
N_time_step = float(lines0 [3])
Box_Length = float(lines0 [5])
epsilon = float(lines0 [7])
sigma = float(lines0 [9])
Mass_Argon = float(lines0 [11])
K_B = float(lines0 [13])
time_step = float(lines0 [15])
targetT = float(lines0 [17])



xpos=[]
ypos=[]
zpos =[]

N_particles = int (N_particles)

for line in lines:
    parts = line.split() 
    if len(parts) > 1:
		xpos.append (float(parts [1]))
		ypos.append (float(parts [2]))
		zpos.append (float(parts [3]))


N_time_step=int(N_time_step)
for n in range (1,N_time_step):
	for i in range (0, N_particles):
		xpos.pop(0)
		ypos.pop(0)
		zpos.pop(0)
		print n*100/N_time_step , "%" , "\r",

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


ax.scatter(xpos, ypos, zpos, c='r', marker='o')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

