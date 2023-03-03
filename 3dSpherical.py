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

par1=[] #Inizialize empty list

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

print("The critical radius is ", np.pi*np.sqrt(mu/eta) , "m")
print(" ")
r1=np.pi*np.sqrt(mu/eta)+0.001   #Above the critical radius


integrand = lambda R: (2/r1)*R*(1-((R/r1)**2))*np.sin(p*np.pi*R/r1)  #Function f(x) we want to integrate

a_pValue=[]  #Only the values of the coefficients
a_pError=[]  #Only the errors of the coefficients

threshold = 1E-6 
#Need to compute a_p coefficients above the threshold


coeff = open("coeffSphericalD", "w+") #Write the coefficients in a file

    
for p in range(1, 45, 2):
    value, error = integrate.quad(integrand, 0, r1) #Compute the integral if it is bigger than the threshold
    a_pValue.append(value)
    a_pError.append(error)
    print("a_{}".format(p), "=", value,"+-", error)
    coeff.write("a_{} = {} +- {} \n".format(p, value, error)) 

coeff.close()
    
def n(r,t): #Define the neutron density
    n = 0
    for p in range(len(a_pValue)):
        n += a_pValue[p]/r * np.exp(eta*t-mu*(p**2*np.pi**2*t/r1**2))*np.sin(p*np.pi*r/r1)
    return abs(n)  


#We plot the results

R_arr=np.linspace(-r1, r1, 100)  #Range on r
T_arr=np.linspace(0, 2E-6, 100)  #Range on t


R, T = np.meshgrid(R_arr, T_arr)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(R, T, n(R, T), cmap =cm.Spectral_r,\
                edgecolor ='black', color = "white")

ax.set_xlabel('r')
ax.set_ylabel('t')
ax.set_zlabel('n(r,t)')

#plt.savefig('3dSpherical.jpeg', facecolor=fig.get_facecolor(), transparent=True)  
plt.show()