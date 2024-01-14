#this program should calculate the temperature distribution in a nuclear fuel rod 

import numpy as np
import matplotlib.pyplot as plt

fuelrad=8.19e-3/2                          #m
gapthickness=0.18e-3                       #m
cladthickness=0.57e-3                      #m
centeraltemp=2573                          #k
boundary1temp=2573                         #k             FUEL-GAP boundary
boundary2temp=300+273                      #k             temp at the outter radius of the clad boundary
Q = 164700*5/(np.pi*(8.19e-3/2)**2*3.66)   #Watt/m**3     volumetric heat generation


#============ intialize for r's and T's ================================

r1 = np.linspace(1e-100,fuelrad, 1000)# nodalisation 
T1 = np.zeros(np.size(r1)) # number of elements 
r2 = np.linspace(fuelrad, fuelrad+gapthickness, 1000)
T2 = np.zeros(np.size(r2))
r3 = np.linspace(fuelrad+gapthickness, fuelrad+gapthickness+cladthickness, 1000)
T3 = np.zeros(np.size(r3))

#=========================FUEL==========================================

h1 = r1[1] - r1[0]                         #step inside the fuel

#initial conditions 
T1[0] = T1[1] = centeraltemp
 
i=0 

while i < np.size(r1)-2:
    
    T1[i+2] = (T1[i+1]*(h1/r1[i]+2) - T1[i] - h1**2*Q/10.2)/(1+h1/r1[i])
    #print(T1[i+2])    
    i += 1

#==========================GAP===========================================

i=0 

T2[0]=T1[-1]
T2[-1]=400+273
slope=(T2[0]-T2[-1])/(r2[0]-r2[-1])
h2=r2[1]-r2[0]
T2[1]=slope*h2+T2[0]

while i < np.size(r2)-2:
    
    T2[i+2] = (T2[i+1]*(h2/r2[i]+2) - T2[i])/(1+h2/r2[i])
    i += 1

#==========================CLAD===========================================


i=0 

T3[0]=T2[-1]
T3[-1]=300+273
slope=(T3[0]-T3[-1])/(r3[0]-r3[-1])
h3=r3[1]-r3[0]
T3[1]=slope*h3+T3[0]

while i < np.size(r2)-2:
    
    T3[i+2] = (T3[i+1]*(h3/r3[i]+2) - T3[i])/(1+h3/r3[i])
    i += 1

#========================KELVEN-TO-DegC==================================

i=0

while i != np.size(r2):
    
    T1[i] = T1[i]-273
    T2[i] = T2[i]-273
    T3[i] = T3[i]-273
    i += 1

#=======================PLOTTING THE RESULTS=================================

plt.plot(r1,T1)
plt.plot(r2,T2)
plt.plot(r3,T3)
plt.xlabel('radius "m"')
plt.ylabel('Temperature "degC"')
plt.grid()
plt.show()

