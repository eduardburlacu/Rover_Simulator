import numpy as np
import matplotlib.pyplot as plt

def Euler_integrator(dt=0.01,corrected=False):
    global k
    global m
    global v0
    t=np.arange(0,1000,dt)
    x=np.zeros(len(t))
    v=np.zeros(len(t))
    v[0]=v0
    for i in range(len(t)-1):
        a=-k*x[i]/m
        v[i+1]=v[i]+a*dt
        if corrected: x[i+1]=x[i]+v[i]* dt
        else: x[i+1]=x[i]+v[i]* dt+a*dt**2/2
    return t,x,v
def Verlet_integrator(dt=0.01):
    global k
    global m
    global v0
    t=np.arange(0,1000,dt)
    x=np.zeros(len(t))
    v=np.zeros(len(t))
    v[0]=v0
    x[1]=x[0]+v0*dt
    for i in range(2,len(t)):
        x[i]=(2-k/m* dt**2)*x[i-1]-x[i-2]
        v[i-1]=(x[i]-x[i-1])/dt
    v[len(t)-1]=(x[len(t)-1]-x[len(t)-2])/dt
    return x,v
dt=2 #This is the critical value after which x starts diverging. A goog usage value is dt=0.12.
k=1
m=1
v0=1
t,x_e,v_e= Euler_integrator(dt=dt)
x_v,v_v=Verlet_integrator(dt=dt)

x_true=v0 * np.sqrt(m/k) * np.sin(np.sqrt(k/m)*t)
#plt.plot(t,x_e,label='Euler Integrator')
plt.plot(t,x_true,label='Analytical Solution')
plt.plot(t,x_v,label='Verlet Integrator')
#plt.ylim(-1.5,1.5)
plt.grid(visible=True)
plt.xlabel('Time')
plt.ylabel('Position')
plt.legend()
plt.show()
