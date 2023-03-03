# Data manipulation
import pandas as pd
import numpy as np
import scipy as sp
from scipy import integrate

# Data viz
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
from matplotlib import cm

# Apply some styling
plt.style.use("ggplot")
rcParams['figure.figsize'] = (12, 6)


###############################################################################################################################

par1=[]

#Opens the file containing parameters
f = open("parametersUPu.txt", "r") 
for line in f:
        index = line.split()[-1]  #Goes into new line after getting the value of the index
        param1 = float(index)  #Set the value of the parameter a float
        par1 = np.append(par1,param1) #For every iteration write the new variable
lambda_t, lambda_f, v, tau, nu  = par1  #List of paramerters name
f.close()

eta=v*(nu-1)/lambda_f
mu=lambda_t*v/3

#In order to compute the critical radius we need to find the zeros of the function:
def equation(R):
    eqr=-1+(R*np.sqrt(eta/mu)*(np.tan(R*np.sqrt(eta/mu)))**(-1))+((3/2)*(R/lambda_t))   
    return eqr

print("The critical radius is ", sp.optimize.fsolve(equation,0.1) , "m")
print(" ")
r1=(sp.optimize.fsolve(equation,0.1)+0.001) #Above the critical radius

def kap(k): #Substitute r1 in the previous equation and solve for k
    eqk=-1+k*r1*(1/np.tan(k*r1))+(3/2)*(r1/lambda_t)
    return eqk
k=sp.optimize.fsolve(kap,10)

alpha=mu*k**2-eta 
A=r1/np.sin(k*r1)


def n(r,t): #define the neutron density
    n = A*np.exp(-alpha*t)*np.sin(k*r)/r
    return n


#We plot the results

R_arr=np.linspace(-r1, r1, 100) #Range on r
T_arr=np.linspace(0, 3E-6, 100) #Range on t


R, T = np.meshgrid(R_arr, T_arr)

fig = plt.figure()

ax = plt.axes(projection='3d')

ax.plot_surface(R, T, n(R, T), cmap =cm.Spectral_r,\
                edgecolor ='black', color = "white")

ax.set_xlabel('r')
ax.set_ylabel('t')
ax.set_zlabel('n(r,z)')

#fig.savefig('3dSphericalNeumann.jpeg', facecolor=fig.get_facecolor(), transparent=True)

plt.show()