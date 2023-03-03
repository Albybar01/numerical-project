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

print("The critical lenght is ", np.pi*np.sqrt(2*mu/eta), "m" )
print(" ")
L = np.pi*np.sqrt(2*mu/eta)+0.001  #Above the critical lenght  L_supercritical

def f_2D(X,Y):  #Initial conditions for 2D
    f_ic = (16*X*Y/(L**2))*(1-(X/L))*(1-(Y/L))
    return f_ic

integrand = lambda X,Y:(4/(L**2))*f_2D(X,Y)*np.sin(p*np.pi*X/L)*np.sin(q*np.pi*Y/L) 
#Function f(x) we want to integrate

a_pqValue=[] #Only the values of the coefficients
a_pqError=[] #Only the errors of the coefficients

a_qValue=[] #Only the values of the coefficients at fixed q
a_qError=[] #Only the errors of the coefficients at fixed q

threshold = 1E-6

#Need to compute a_pq coefficients above the threshold


coeff = open("coefficients2d", "w+") 
#Write the coefficients in a file

print("The coefficients of a_pq are:")
#Coefficients a_pq where p and q are int number pq starts from 1 with step 2 (p=p+2) because integral is 0 when pq is even

for p in range(1, 45, 2):
    
    for q in range(1, 45, 2):
        value, error = integrate.dblquad(integrand, 0, L, 0, L) #Compute the integral if it is bigger than the threshold
        a_qValue.append(value)
        a_qError.append(error)
        print("a_{},{}".format(p,q), "=", value,"+-", error)
        coeff.write("a_{},{} = {} +- {} \
                               \n".format(p, q, value, error)) 
    
        a_pqValue.append(a_qValue) 
        a_pqError.append(a_qError) 
    
        if np.abs(value) < threshold:
            break
        
            if len(a_qValue) == 0: break
        
        else: 
            a_qValue = [] 
            a_qError = []
        
coeff.close()

t=1E-7  #Time is fixed at
def n(x, y): #Define the neutron density
    n = 0
    for i in range(len(a_pqValue)):
        P = 2 * i + 1
        for j in range(len(a_pqValue[i])):
            Q = 2 * j + 1
            n += a_pqValue[i][j] * np.exp(eta*t -mu*np.pi**2*((P**2)/(L**2)+ (Q**2)/(L**2))*t) * np.sin(P*np.pi*x/L)* np.sin(Q*np.pi*y/L)
    return n
#Plotting the results

X_arr=np.linspace(0, L, 1000)
#Range on x 
Y_arr=np.linspace(0, L, 1000)
#Range on y
X, Y = np.meshgrid(X_arr, Y_arr)
fig = plt.figure()

ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, n(X, Y), cmap =cm.Spectral_r,\
                edgecolor ='black', color = "white")

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('n(x,y)')

#plt.savefig('2dCartesian.jpeg', facecolor=fig.get_facecolor(), transparent=True)  
plt.show()