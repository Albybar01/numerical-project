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

L = np.pi*np.sqrt(3*mu/eta)+0.001  #Above the critical lenght  L_supercritical

#Above the critical radius
r1=np.sqrt((eta*L**2-np.pi**2*mu)*mu)*L*sp.special.jn_zeros(0, 1)/(eta*L**2 - np.pi**2*mu)+0.001 

print("The critical length is ", np.pi*np.sqrt(3*mu/eta), "m")
print(" ")

print("The critical radius is ", np.sqrt((eta*L**2-np.pi**2*mu)*mu)*L*sp.special.jn_zeros(0, 1)/(eta*L**2 - np.pi**2*mu) ,"m")
print(" ")

alpha=sp.special.jn_zeros(0, 100) #List with the zeros of first type bessel function

#Function f(x) we want to integrate
integrand = lambda Z,R: (4/(L*(r1**2)*(sp.special.jv(1, alpha[q-1])**2)))*sp.special.jv(0, alpha[q-1]*R/r1)*R*\
                        (1-((R**2)/(r1**2)))*(np.sin(np.pi*Z/L))**2     
                         
    
a_1qValue=[] #Only the values of the coefficients
a_1qError=[] #Only the errors of the coefficients


threshold = 0.0001 
#Need to compute a_p coefficients above the threshold


coeff = open("coeffCilindrical", "w+") #Write the coefficients in a file


print("The coefficients of a_1q are:")
#Coefficients a_1q where q is an int number q starts from 1 with step 2 (p=p+2) because integral is 0 when q is even


p = 1
    
for q in range(1, 45, 2):
    value, error = integrate.dblquad(integrand, 0, r1, 0, L) #Compute the integral if it is bigger than the threshold
    a_1qValue.append(value)
    a_1qError.append(error)
    print("a_{},{}".format(1,q), "=", value,"+-", error)
    coeff.write("a_{},{} = {} +- {} \n".format(1, q, value, error)) 
        
coeff.close()
    
t=1E-5 #we fix a particular time 
def n(r, z): #define the neutron density
    n = 0
    for q in range(len(a_1qValue)):
        n += a_1qValue[q] * sp.special.jv(0, alpha[q-1]*r/r1)*np.sin(np.pi*z/L)*np.exp(((eta*(r1**2)*(L**2)-mu*((alpha[q-1]**2)*\
        (L**2)+(np.pi**2)*(r1**2))) / ((r1**2)*(L**2)))*t )
    return abs(n)
    
    
#We plot the results

R_arr=np.linspace(-r1, r1, 100) #Range on r
Z_arr=np.linspace(0, L, 100) #Range on z


R, Z = np.meshgrid(R_arr, Z_arr)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(R, Z, n(R, Z), cmap =cm.Spectral_r,\
                edgecolor ='black', color = "white")

ax.set_xlabel('r')
ax.set_ylabel('z')
ax.set_zlabel('n(r,z)')

#plt.savefig('3dCylindrical.jpeg', facecolor=fig.get_facecolor(), transparent=True)  
plt.show()