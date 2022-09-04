import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm as norm

def create_excentricity(scenario):
    if scenario==1 or scenario==2:
        return None
    else:
        e=float(input('Enter excentricity(any number >=0) -->'))
        return e

def create_v0(x0,scenario,e):
    ''''
    Select scenario and desired excentricity of the orbit. Outputs the desired initial velocity as 3x1 numpy array.
    '''
    if scenario==1: return np.array([0,0,0])
    v_circular=np.sqrt(k/norm(x0))
    if scenario==2: return np.array([0, 0, v_circular])
    if e is not None: return np.array([0,0,  np.sqrt(k / norm(x0) * (1 + e))])
    else: raise ValueError

def create_period(v0,x0):
    a=1/( 2/norm(x0)- norm(v0)**2/k)
    return 2* np.pi * np.sqrt( a**3/k)

def sample_time(v0,x0):
    if e<1:
        return create_period(v0,x0)
    else: return int(3* np.pi * np.sqrt(k/norm(x0)))

def Euler_integrator_3d(v0,x0,dt=10):

    t=np.arange(0,sample_time(v0,x0),dt)
    x=np.zeros((len(t),3))
    v=np.zeros((len(t),3))
    v[0]=v0
    x[0]=x0
    for i in range(len(t)-1):
        a = -k / ((norm(x[i])) ** 3) * x[i]
        v[i+1]=v[i]+a*dt
        x[i+1]=x[i]+v[i]* dt
    return t,x,v

def Verlet_integrator_3d(v0,x0,dt=10):

    t=np.arange(0,sample_time(v0,x0),dt)
    x=np.zeros((len(t),3))
    v=np.zeros((len(t),3))
    v[0]=v0
    x[0]=x0
    x[1]=x0+v0*dt
    for i in range(2,len(t)):
        a = -k / ((norm(x[i-1])) ** 3) * x[i-1]
        x[i]=2*x[i-1]-x[i-2]+a*dt**2
        v[i-1]=(x[i]-x[i-1])/dt
    v[len(t)-1]=(x[len(t)-1]-x[len(t)-2])/dt
    return x,v

if __name__=='__main__':

    G=6.67E-11
    M=6.42E23
    k=G*M
    scenario=int(input('Enter Scenario-->'))
    e=create_excentricity(scenario)
    x0=np.array([1E8/np.sqrt(2),1E8/np.sqrt(2),0])
    v0=create_v0(x0=x0,scenario=scenario,e=e)
    if scenario!=1:
        spatial = int(input('Enter number of dimensions of the desired plot (2 or 3)--> '))
        if spatial not in {2, 3}: raise ValueError
    else:spatial=None
    t,x_e,v_e=Euler_integrator_3d(v0,x0)
    x_v,v_v=Verlet_integrator_3d(v0,x0)

    if scenario==1:
        plt.plot(t, x_e[:, 0], label='Euler Integrator 3D')
        plt.plot(t, x_v[:, 0], label='Verlet Integrator 3D')
        plt.xlim(0,0.25E6)
        plt.ylim(-1.75E8,1.1E8)
        plt.xlabel('t')
        plt.ylabel('x')

    else:
        plt.scatter(0, 0, label='Mars', c='red')
        if spatial==2:
            plt.plot(x_e[:, 0], x_e[:, 2], color='blue',label='Euler Integrator 3D')
            plt.plot(x_v[:, 0], x_v[:, 2], color='green',label='Verlet Integrator 3D')
            plt.xlabel('x')
            plt.ylabel('z')
            plt.legend()
            plt.show()

        else:
            import pandas as pd
            import plotly.express as px
            df=pd.DataFrame(x_v,columns=['x','y','z'])
            df.head(n=12)
            fig=px.line_3d(df,'x','y','z')
            fig.show()
