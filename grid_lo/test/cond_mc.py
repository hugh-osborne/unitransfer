# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 16:20:10 2016

@author: scsyml
"""

#import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ
import time
from scipy.stats import poisson
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('T_max', nargs='?', default=0.2, type=float)
parser.add_argument('T_report', nargs='?', default=1e-3, type=float)
args = parser.parse_args()


#import writemesh
#import uuid
#import mesh

start=time.time() #start time

# Constants taken from Apfaltrer et al. (2006). The constants come from Appendix B
param_dict={

'tau_m':   20e-3,
'E_r': -65e-3,
'E_e':  0e-3,
'tau_s':  5e-3,
'g_max': 0.8,
'V_min':-66.0e-3,
'V_max':  -55.0e-3,
'V_th': -55.0e-3, #'V_max', # sometimes used in other scripts
'N_V': 200,
'w_min':0.0,
'w_max':  10.0,
'N_w': 20,

}

def Cond(V,g, t, param_dict):
    '''Derivatives of the Cond model'''


    p = (-(V-param_dict['E_r'])  - g*(V-param_dict['E_e']))/param_dict['tau_m']
    q = -g/param_dict['tau_s']
    return p, q


def simulation(rate, t_end = 0.2, g_max = 10.0, h = 0.05):
    dt=5e-6 #time step
    nu=rate #input rate




    N_pop=10000 #number of neurons

    T_max = t_end


    out_step = args.T_report
    factor=int(np.around(out_step/dt))



    N_steps=int(np.around(T_max/dt)) #number of time-steps)

    V_save=np.zeros((N_pop,N_steps/factor))
    g_save=np.zeros((N_pop,N_steps/factor))



    V=param_dict['E_r']*np.ones(N_pop)
    g=param_dict['w_min']*np.ones(N_pop)
#initial conditions

    t=0
    g_m = g_max*np.ones(N_pop)

    out=np.zeros(N_steps)
    for i in range(N_steps):

        [dV,dg]=Cond(V,g,t,param_dict)
        spikes=poisson.rvs(nu*dt,size=N_pop) #generating spikes

        V=V+dt*dV
        gs=g+dt*dg+h*spikes
        g = np.minimum(gs,g_m)
    #euler step
        if i % factor == 0:
            V_save[:,i/factor]=V
            g_save[:,i/factor]=g

        out[i]=np.sum(V>param_dict['V_th'])   #calculate number of firing neurons
        V[V>param_dict['V_th']]=param_dict['E_r']  #firing neurons reset


        t=t+dt


#binning firing neurons to get a smooth firing rate curve
#factor=50 #number of time-steps in each bin
    rate_smooth=np.zeros(N_steps/factor)
    for j in range(N_steps/factor):
        rate_smooth[j]=np.sum(out[factor*j:factor*(j+1)-1])/dt/factor/N_pop

    t_axis=dt*factor*np.arange(N_steps/factor)

    f= open('rates'+ '_' + str(nu),'w')
    for i, t in enumerate(t_axis):
        f.write(str(t_axis[i]) + '\t' + str(rate_smooth[i]) + '\n')

    end=time.time()
    print "Runtime for simulation is %(t)fs" %{'t':end-start}

    plt.plot(t_axis,rate_smooth)
    plt.show()
    title = str(nu) + '_' + str(g_max)
    np.savez(title,t_axis=t_axis,rate_smooth=rate_smooth,V_save=V_save,w_save=g_save,out_step=out_step)

    with open(title + '.rates','w') as rf:
        for t in t_axis:
            rf.write(str(t) + ' ')
        rf.write('\n')
        for f in rate_smooth:
            rf.write(str(f) + ' ')
        rf.write('\n')

simulation(5000)
