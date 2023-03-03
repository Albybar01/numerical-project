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
print("The critical lenght is ",np.pi*np.sqrt(3*mu/eta) , "m" )
print(" ")

def f_3D(X,Y,Z):  #Initial conditions of f for 3D
    f_ic = (8*X*Y*Z/(L**3))*(1-(X/L))*(1-(Y/L))*(1-(Z/L))
    return f_ic

integrand = lambda X,Y,Z: (8/(L**3))*f_3D(X,Y,Z)*np.sin(p*np.pi*X/L)*np.sin(q*np.pi*Y/L)*np.sin(r*np.pi*Z/L)
#Function f(x) we want to integrate

a_pqrValue=[] #Only the values of the coefficients
a_pqrError=[] #Only the errors of the coefficients

a_qrValue=[] #Only the values of the coefficients at fixed p
a_qrError=[] #Only the errors of the coefficients at fixed p

a_rValue=[] #Only the values of the coefficients at fixed p and q
a_rError=[] #Only the errors of the coefficients at fixed p and q

threshold = 1E-6

#Need to compute a_p coefficients above the threshold


print("The coefficients of a_pqr are:")
#Coefficients a_pqr where p is an int number p starts from 1 with step 2 (p=p+2) because integral is 0 when p is even


coeff = open("coefficients3d", "w+") 
#Write the coefficients in a file


for p in range(1, 13, 2):
    
    for q in range(1, 13, 2):
        
        for r in range(1, 13, 2):
            value, error = integrate.tplquad(integrand, 0, L, 0, L, 0, L) #Compute the integral if it is bigger than the threshold
            if np.abs(value) < threshold:
                break
            a_rValue.append(value)
            #Add the value of the integral to the list
            a_rError.append(error)
            #Add the error of the integral to the list 
            print("a_{},{},{}".format(p,q,r), "=", value,"+-", error)  #Show the results
            coeff.write("a_{},{},{} = {} +- {} \n".format(p, q, r, value, error))

        
    if len(a_rValue) == 0: break
        
    a_qrValue.append(a_rValue)
    a_qrError.append(a_rError)
            
    a_rValue = []
    a_rError = []
    #Reset the two lists  
    
    if len(a_qrValue) == 0: break

    a_pqrValue.append(a_qrValue)
    a_pqrError.append(a_qrError) 
    a_qrValue = []
    a_qrError = []
    #Reset the two lists       
        
coeff.close()


t=2E-7 #Time is fixed at

def n(x, y, z): #Define the neutron density
    n = 0
    for i in range(len(a_pqrValue)):
        P = 2 * i + 1
        for j in range(len(a_pqrValue[i])):
            Q = 2 * j + 1
            for k in range(len(a_pqrValue[i][j])):
                R = 2 * k +1            
                n += a_pqrValue[i][j][k] *np.exp(eta*t -mu*np.pi**2*((P**2)/(L**2) + (Q**2)/(L**2) + (R**2)/(L**2))*t) * \
                    np.sin(P*np.pi*x/L)* np.sin(Q*np.pi*y/L)* np.sin(R*np.pi*z/L)
    return n

Z=L/2 #We fix a particular value of Z
#Plotting the results

X_arr=np.linspace(0, L, 1000)
#Range on x 
Y_arr=np.linspace(0, L, 1000)
#Range on y

X, Y = np.meshgrid(X_arr, Y_arr)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, n(X, Y, Z), cmap =cm.Spectral_r,\
                edgecolor ='black', color = "white")

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('n(x,y)')

#plt.savefig('3dCartesian.jpeg', facecolor=fig.get_facecolor(), transparent=True)  
plt.show()