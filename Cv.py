import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt
from decimal import Decimal as D
import sys


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

