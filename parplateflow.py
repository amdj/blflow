#!/usr/bin/python

# Boundary layer flow
from numpy import *
import time 
import matplotlib
matplotlib.use('TkAgg')
# from matplotlib.pylab import *
import pylab as p



# import matplotlib.animation as animation
def K(t):                                 #Forcing function
    return (1-exp(-0.1*t))*cos(t)


s=10
#Define domain
n=50                                      #Number of gridpoints
y=linspace(-1,1,n)

dy=y[1]-y[0]
dt=0.0005

l=(dt/(s**2*dy**2))

hnu=cosh(sqrt(1j)*s*y)/cosh(sqrt(1j)*s)
fnu=tanh(sqrt(1j)*s)/(sqrt(1j)*s)


def u_ex(tn):
    return (((1-hnu)/(1-fnu))*exp(1j*(tn))/1j).real
    

def u_np1(un,tn,dt):
    Kn=K(tn)
    unp1=un
    unp1[0]=0                                  #Velocity zero ver here
    unp1[1:-1]=un[1:-1]+dt*Kn+l*(un[0:-2]-2*un[1:-1]+un[2:])
    unp1[-1]=0                     #Boundary other side
    return unp1

un0=zeros(n,float)
t=0
un=un0
# un.append(un0)

# Make the plot
p.ion()
linefd, = p.plot(un0,y)
linee, = p.plot(un0,y)
p.legend(('Finite difference','Periodic exact'))
p.ylim(-1,1)
p.xlim(-1.5,1.5)
p.ylabel('y')
p.xlabel('u')
p.grid('on')
i=0
uold=un
while(True):
    t+=dt
    uold=un
    un=u_np1(uold,t,dt)
    if(i%100==0):
        linefd.set_xdata(un)
        linee.set_xdata(u_ex(t))
        p.draw()
        # print("Time:",t)
    i+=1




