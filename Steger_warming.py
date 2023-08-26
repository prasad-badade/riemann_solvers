import numpy as np
from numpy.linalg import eig,inv
import matplotlib.pyplot as plt
import pandas as pd

gamma=1.4
n = 400
cfl = 0.8
dx=(1.0/400)

rho_l,P_l,u_l=1.0,1.0,0.0
rho_r,P_r,u_r=0.125,0.1,0.0

e_l=u_l**2/2.0+P_l/(rho_l*(gamma-1.0))
e_r=u_r**2/2.0+P_r/(rho_r*(gamma-1.0))

#  U vector =[rho,rho*u,rho*e]
# giving initial value
U=np.array([[rho_l,rho_l*u_l,rho_l*e_l]]*(200)+[[rho_r,rho_r*u_r,rho_r*e_r]]*(200)) 

U_new=np.copy(U)

t=0
rho=[0]*n
u=[0]*n
p=[0]*n
e=[0]*n
H=[0]*n 
a=[0]*n 
itr=0
gamma=1.4  
lmdplus=[0]*3
lmdminus=[0]*3
f_plus=[[0]*3 for i in range (n)]
f_minus=[[0]*3 for i in range (n)]
f_p=[[0]*3 for i in range (n)]
f_m=[[0]*3 for i in range (n)]
lamda=np.array([0.0]*3)

while(t<0.15):
    l_max=0.0
    for k in range(n):                      #rho,u and p values
        rho[k]=U[k][0]
        u[k]=U[k][1]/U[k][0]
        p[k]=(U[k][2]/U[k][0]-u[k]**2 /2)*(gamma-1)*rho[k]
        
        e[k]=((u[k]**2)/2)+ (p[k]/((gamma-1)*rho[k]))
        H[k]=(rho[k]*e[k]+p[k])/rho[k]
        
        a[k]=(gamma*p[k]/rho[k])**(0.5)
        #lamda is notation for eigen values
        lamda[0]=u[k]-a[k]
        lamda[1]=u[k]
        lamda[2]=u[k]+a[k]
        l_max=max(l_max,max(abs(lamda)))
        
        for i in range(3):             #finding lamdaplus and lamdaminus
            lmdplus[i]=max(lamda[i],0)
            lmdminus[i]=min(lamda[i],0)      
        
        f_plus[k][0]=rho[k]/(2*gamma)*(lmdplus[0]+2*(gamma-1)*lmdplus[1]+lmdplus[2])
        f_plus[k][1]=rho[k]/(2*gamma)*(lmdplus[0]*(u[k]-a[k])+2*(gamma-1)*u[k]*lmdplus[1]+lmdplus[2]*(u[k]+a[k]))
        f_plus[k][2]=rho[k]/(2*gamma)*(lmdplus[0]*(H[k]-u[k]*a[k])+(gamma-1)*u[k]**2*lmdplus[1]+lmdplus[2]*(H[k]+u[k]*a[k]))
        f_minus[k][0]=rho[k]/(2*gamma)*(lmdminus[0]+2*(gamma-1)*lmdminus[1]+lmdminus[2])
        f_minus[k][1]=rho[k]/(2*gamma)*(lmdminus[0]*(u[k]-a[k])+2*(gamma-1)*u[k]*lmdminus[1]+lmdminus[2]*(u[k]+a[k]))
        f_minus[k][2]=rho[k]/(2*gamma)*(lmdminus[0]*(H[k]-u[k]*a[k])+(gamma-1)*u[k]**2*lmdminus[1]+lmdminus[2]*(H[k]+u[k]*a[k]))
    
    dt=cfl*dx/l_max    
    for k in range(1,n-1):
        for i in range(3):
            f_p[k][i]=f_plus[k][i]+f_minus[k+1][i]
            f_m[k][i]=f_plus[k-1][i]+f_minus[k][i]
        for i in range(3):
            U_new[k][i]=U[k][i]-(dt/dx)*(f_p[k][i]-f_m[k][i])
   
    U=U_new.copy()
    itr+=1
    t+=dt

x_plot=np.linspace(0,1,n)

std_data=pd.read_csv("StandardSod_0_15.txt",sep=' ')
x_std=np.array(std_data['x'])
rho_std=np.array(std_data['rho'])
u_std=np.array(std_data['u'])
P_std=np.array(std_data['p'])
internal_energy_std=np.array(std_data['ie'])


den=[[0]*1 for i in range(n)]
for k in range(400):
    den[k][0]=U[k][0]
plt.scatter(x_plot,den,marker='.',label='Steger-Warming')
plt.plot(x_std,rho_std,label='Analytic sol')
plt.legend()
plt.xlabel("x")
plt.ylabel("rho")
plt.title("Density")
plt.show()


x_plot=np.linspace(0,1,n)
uu=[[0]*1 for i in range(n)]
for k in range(400):
    uu[k]=U[k][1]/U[k][0]
plt.scatter(x_plot,uu,marker='.',label='Steger-Warming')
plt.plot(x_std,u_std,label='Analytic sol')
plt.legend()
plt.xlabel("x")
plt.ylabel("velocity")
plt.title("Velocity")
plt.show()

x_plot=np.linspace(0,1,n)
pp=[[0]*1 for i in range(n)]
for k in range(400):
    pp[k]=(U[k][2]/U[k][0]-u[k]**2 /2)*(gamma-1)*rho[k]
plt.scatter(x_plot,pp,marker='.',label='Steger-Warming')
plt.plot(x_std,P_std,label='Analytic sol')
plt.legend()
plt.xlabel("x")
plt.ylabel("pressure")
plt.title("Pressure")
plt.show()

x_plot=np.linspace(0,1,n)
ee=[[0]*1 for i in range(n)]
for k in range(400):
    ee[k]=(U[k][2]/U[k][0]-u[k]**2 /2)
plt.scatter(x_plot,ee,marker='.',label='Steger-Warming')
plt.plot(x_std,internal_energy_std,label='Analytic sol')
plt.legend()
plt.xlabel("x")
plt.ylabel("internal Energy")
plt.title("Internal Energy")
plt.show()

    