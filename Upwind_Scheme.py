import numpy as np
import matplotlib.pyplot as plt
import math as math

dx=0.02
dt=0.001
u=[0]*1000
u_new=[0]*1000
def fx(x):
    return ((x)**2)/2
uplus,uminus=0,0
y=-10
z=0
x=[0]*1000
pi=3.142
for i in range (-500,0,1):
    x[i]=y
    # print("x",x[i])
    y=y+dx
for i in range (0,500,1):
    x[i]=z
    z=y+dx
    
    
for i in range (-500,0):
    u[i]=((math.exp(-((x[i]-2)**2)/2))/(2*pi)**(0.5))-((math.exp(-((x[i]+2.0)**2)/2))/(2*pi)**(0.5))
    # print(u[i])
   #u[i]=1
for i in range (0,500):
   # u[i]=0 
   u[i]=((math.exp(-((x[i]-2)**2)/2))/(2*pi)**(0.5))-((math.exp(-((x[i]+2.0)**2)/2))/(2*pi)**(0.5))
t=0
while t<20:
    for i in range (-499,500,1):
        s=(u[i]+u[i+1])/2
        if u[i]>u[i+1]:
            if s>0:
                uplus=u[i]
            else:
                uplus=u[i+1]
        else:
            if u[i]<=0:
                uplus=u[i];
            elif (u[i]<0 and 0<u[i+1]):
                uplus=0
            elif u[i+1]>=0:
                uplus=u[i+1]
        if u[i-1]>u[i]:
            s=(u[i-1]+u[i])/2
            if s>0:
                uminus=u[i-1]
            else:
                uminus=u[i]
        else:
            if u[i]<=0:
                uminus=u[i-1];
            elif (u[i]<0 and 0<u[i+1]):
                uminus=0
            elif u[i]>=0:
                uminus=u[i+1]
        u_new[i]=u[i]-((dt/dx)*(fx(uplus)-fx(uminus)))
        t+=dt
    for i in range (-500,500,1):
        u[i]=u_new[i]
x=np.linspace(-10,10,1000)            
plt.plot(x,u)    
            
            
        
            
