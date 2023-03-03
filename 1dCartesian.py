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

#calculating eta and nu
eta=v*(nu-1)/lambda_f
mu=lambda_t*v/3

A = 1  #Parameters A and Lambda using gaussian initial condition
Lambda = 100  
    
print("The critical lenght is ", np.pi*np.sqrt(mu/eta), "m" )
print(" ")
L = np.pi*np.sqrt(mu/eta)+0.001  #Above the critical lenght  L_supercritical

def f(X):  # Gaussian initial conditions for f(x)
    f_ic = A * np.exp((-4 * Lambda * (X - 0.5 * L) ** 2 / (L ** 2)))
    return f_ic


integrand = lambda X: (2 / L) * f(X) * np.sin(p * np.pi * X / L) 
#Function f(x) we want to integrate

threshold = 1E-6 
#Need to compute a_p coefficients above the threshold

coeff = open("coefficients1d", "w+") 
#Write the coefficients in a file named:

print("The coefficients of a_p are:")
a_pValue = [] #Initialized empty 
a_pError = [] 

#Coefficients a_p where p is an int number p starts from 1 with step 2 (p=p+2) because integral is 0 when p is even
for p in range(1, 50, 2):
    value, error = integrate.quad(integrand, 0, L) #Compute the integral if it is bigger than the threshold
    a_pValue.append(value)
    a_pError.append(error)
    print("a_{}".format(p), "=", value,"+-", error)
    coeff.write("a_{} = {} +- {}  \n".format(p, value, error)) 
    
    if np.abs(value) < threshold:
        break
        
coeff.close()
             
def n(x, t): #Define the neutron density as
    n = 0
    for i in range(len(a_pValue)):
        P = 2 * i + 1
        n += a_pValue[i] * np.exp(eta*t - mu*(P**2)*(np.pi**2)*t/(L**2))* np.sin(P * np.pi * x/L)
    return n
#Plotting the results

#Range on x             
X_arr=np.linspace(0, L, 100) 
#Range on t
T_arr=np.linspace(0, 0.00001, 100)

X, T = np.meshgrid(X_arr, T_arr)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, T, n(X, T), cmap =cm.Spectral_r,\
                edgecolor ='black', color = "white")

ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('n(x,t)')

#plt.savefig('1dCartesian.jpeg', facecolor=fig.get_facecolor(), transparent=True)  
plt.show()

