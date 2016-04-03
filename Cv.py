import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from decimal import Decimal as D
import sys


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

N_particles = len(xpos)/len(time)	

#Calculation of Cv

remove = input ("How many time steps do you want to remove from the beginning?")
if remove > len(time):
	sys.exit("Error ")
	

for i in range (0,remove):
	E.pop(0)
	sysT.pop(0)
	
E2=[]

avgE = np.sum(E)/len(E)
for i in range (0,len(E)):
	E2.append(E[i]**2)
avgE2 = np.sum(E2)/len(E2)
avgT = np.sum(sysT)/len(sysT)

Cv = (avgE2 - avgE**2)/(K_B * avgT**2)
print "Cv = ", Cv, "J/kg/K","\n" , "Average Temperature = ", avgT , "K"

