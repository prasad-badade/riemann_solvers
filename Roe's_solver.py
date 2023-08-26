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
gamma=1.4
d=[0]*3
itr=0
rho=[0]*n
u=[0]*n
p=[0]*n
f=[[0,0,0]]*n
F=[[0,0,0]]*(n-1)
eps=1e-6
while(t<0.15):
    l_max=0
    for k in range(n):      #rho,u and p values
        rho[k]=U[k][0]       
        u[k]=U[k][1]/U[k][0]
        p[k]=(U[k][2]/U[k][0]-u[k]**2 /2)*(gamma-1)*rho[k]
    for k in range(1,400):
        gamma=1.4
        e_r=((u[k]**2)/2)+ (p[k]/((gamma-1)*rho[k]))
        e_l=((u[k-1]**2)/2)+(p[k-1]/((gamma-1)*rho[k-1]))
        H_r=(rho[k]*e_r+p[k])/rho[k]
        H_l=(rho[k-1]*e_l+p[k-1])/rho[k-1]
        
        #finding u_delta,H_delta,a_delta
        u_delta=((((rho[k])**0.5)*u[k])+(((rho[k-1])**0.5)*u[k-1]))/(((rho[k])**0.5)+((rho[k-1])**0.5))
        H_delta=((((rho[k])**0.5)*H_r)+(((rho[k-1])**0.5)*H_l))/(((rho[k])**0.5)+((rho[k-1])**0.5))
        a_delta=(gamma-1)*(H_delta-((u_delta**2)/2))
        # f is a flux vector
        f[k-1]=[rho[k-1]*u[k-1],((rho[k-1]*(u[k-1]**2))+p[k-1]),(u[k-1]*(rho[k-1]*e_l+p[k-1]))]

        # lamda is notation for eigen value
        lamda=[0]*3
        lamda[0]=u_delta-a_delta
        lamda[1]=u_delta
        lamda[2]=u_delta+a_delta
            
        l_max=max(l_max,max(max(lamda),-min(lamda)))
        
        # z is notation for eigen vectors
        z=[[0]*3 for i in range (3)]
        z[0]=[1,u_delta-a_delta,(H_delta-u_delta*a_delta)]
        z[1]=[1,u_delta,((u_delta**2)/2)]
        z[2]=[1,u_delta+a_delta,(H_delta+u_delta*a_delta)]
        z=np.transpose(z)
        
        # findind delta_u1,delta_u2,delta_u3 as du
        du_0=rho[k]-rho[k-1]
        du_1=rho[k]*u[k]-rho[k-1]*u[k-1]
        du_2=rho[k]*e_r-rho[k-1]*e_l
        dell=[0]*3
        
        #finding delta values as dell
        dell[1]=(gamma-1)*((du_0*(H_delta-(u_delta**2))+u_delta*du_1-du_2))/(a_delta**2)
        dell[0]=((du_0*(u_delta+a_delta)-du_1-a_delta*dell[1]))/(2*a_delta)
        dell[2]=du_0-(dell[0]+dell[1])
        F[k-1]=f[k-1].copy()
       
        for j in range(3):
            s = 0
            for i in range(3):
                if(lamda[i]>0):
                    break
                s += lamda[i]*dell[i]*z[j][i]
            F[k-1][j] +=s      #finding flux at face
            
        
    dt=cfl*dx/l_max
    for i in range(1,n-1):
        for j in range(3):
            U_new[i][j]=U[i][j]-(dt/dx)*(F[i][j]-F[i-1][j])
    U=U_new.copy()
    itr+=1
    t+=dt
    
# x_plot=np.linspace(0,1,n)
# den=[Ui[0] for Ui in U]
# plt.scatter(x_plot,den,marker='.')
    
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
plt.scatter(x_plot,den,marker='.',color='#1f77b4',label='roe solver')
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
plt.scatter(x_plot,uu,marker='.',label='roe solver')
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
plt.scatter(x_plot,pp,marker='.',label='roe solver')
plt.plot(x_std,P_std,label='Analytic sol')
plt.legend()
plt.xlabel("x")
plt.ylabel("pressure")
plt.title("pressure")
plt.show()

x_plot=np.linspace(0,1,n)
ee=[[0]*1 for i in range(n)]
for k in range(400):
    ee[k]=(U[k][2]/U[k][0]-u[k]**2 /2)
plt.scatter(x_plot,ee,marker='.',label='roe solver')
plt.plot(x_std,internal_energy_std,label='Analytic sol')
plt.legend()
plt.xlabel("x")
plt.ylabel("internal Energy")
plt.title("Internal Energy")
plt.show()    

