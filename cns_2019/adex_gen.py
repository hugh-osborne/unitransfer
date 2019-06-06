#!/usr/bin/env python3

import numpy as np
from scipy.integrate import odeint
import math
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from numpy.linalg import norm

def adEx(y,t):
    C = 281
    g_l = 30
    E_l = -70.6
    v_t = -50.4
    tau = 2.0
    alpha = 4.0
    tau_w = 144.0

    v = y[0];
    w = y[1];

    v_prime = (-g_l*(v - E_l) + g_l*tau*np.exp((v - v_t)/tau) - w) / C
    w_prime = (alpha*(v - E_l) - w) / tau_w

    return [v_prime, w_prime]

def adExBack(y,t):
    C = 281
    g_l = 30
    E_l = -70.6
    v_t = -50.4
    tau = 2.0
    alpha = 4.0
    tau_w = 144.0

    v = y[0];
    w = y[1];

    v_prime = (-g_l*(v - E_l) + g_l*tau*np.exp((v - v_t)/tau) - w) / C
    w_prime = (alpha*(v - E_l) - w) / tau_w

    return [-v_prime, -w_prime]

def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def isPolyConvex(ps):
    return np.abs(np.arccos(np.dot((ps[1]-ps[0])/norm(ps[1]-ps[0]),(ps[2]-ps[0])/norm(ps[2]-ps[0])))) < np.abs(np.arccos(np.dot((ps[1]-ps[0])/norm(ps[1]-ps[0]),(ps[3]-ps[0])/norm(ps[3]-ps[0])))) and np.abs(np.arccos(np.dot((ps[2]-ps[1])/norm(ps[2]-ps[1]),(ps[3]-ps[1])/norm(ps[3]-ps[1])))) < np.abs(np.arccos(np.dot((ps[2]-ps[1])/norm(ps[2]-ps[1]),(ps[0]-ps[1])/norm(ps[0]-ps[1]))))

def generateMesh():
    timestep = 0.001
    start_v = -80
    start_range_w = range(-720,300,20)

    start_points = [[start_v,w] for w in start_range_w]

    tspan = np.linspace(0, 80, 100)

    fig, ax = plt.subplots()
    lines = []
    for strip in range(len(start_points)-1):
        t_1 = odeint(adEx, start_points[strip], tspan)
        t_2 = odeint(adEx, start_points[strip+1], tspan)

        t_1 = np.array([t for t in t_1[:] if t[0] < -30.0 and t[1] < 300 and t[0] > -80.01 and t[1] > -800])
        t_2 = np.array([t for t in t_2[:] if t[0] < -30.0 and t[1] < 300 and t[0] > -80.01 and t[1] > -800])

        for cell in range(len(t_1[0:])-1):
            p1 = [t_1[cell][0],t_1[cell][1]]
            p2 = [t_1[cell+1][0],t_1[cell+1][1]]
            p3 = [t_2[cell][0],t_2[cell][1]]
            p4 = [t_2[cell+1][0],t_2[cell+1][1]]

            if PolyArea([p1[0],p2[0],p3[0],p4[0]], [p1[1],p2[1],p3[1],p4[1]]) < 0.001:
                break

            if not isPolyConvex(np.array([p1,p2,p4,p3])):
                break

            lines = lines + [[p1,p3],[p4,p2]]

        ax.plot(t_1[:,0], t_1[:,1], color='k')
        ax.plot(t_2[:,0], t_2[:,1], color='k')

    lc = mc.LineCollection(lines, linewidths=2)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)

    ############################

    start_v = -50
    start_range_w = range(-510,300,20)

    start_points = [[start_v,w] for w in start_range_w]

    tspan = np.linspace(0, 50, 100)

    lines = []
    for strip in range(len(start_points)-1):
        t_1 = odeint(adEx, start_points[strip], tspan)
        t_2 = odeint(adEx, start_points[strip+1], tspan)

        t_1 = np.array([t for t in t_1[:] if t[0] < -30.0 and t[1] < 300 and t[0] > -80.01 and t[1] > -800])
        t_2 = np.array([t for t in t_2[:] if t[0] < -30.0 and t[1] < 300 and t[0] > -80.01 and t[1] > -800])

        for cell in range(len(t_1[0:])-1):
            p1 = [t_1[cell][0],t_1[cell][1]]
            p2 = [t_1[cell+1][0],t_1[cell+1][1]]
            p3 = [t_2[cell][0],t_2[cell][1]]
            p4 = [t_2[cell+1][0],t_2[cell+1][1]]

            if PolyArea([p1[0],p2[0],p3[0],p4[0]], [p1[1],p2[1],p3[1],p4[1]]) < 0.001:
                break

            if not isPolyConvex(np.array([p1,p2,p4,p3])):
                break

            lines = lines + [[p1,p3],[p4,p2]]

        ax.plot(t_1[:,0], t_1[:,1], color='k')
        ax.plot(t_2[:,0], t_2[:,1], color='k')

    lc = mc.LineCollection(lines, linewidths=2)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)

    ####################################################

    start_v = -50
    start_range_w = range(-510,300,20)

    start_points = [[start_v,w] for w in start_range_w]

    tspan = np.linspace(0, 50, 100)

    lines = []
    for strip in range(len(start_points)-1):
        t_1 = odeint(adExBack, start_points[strip], tspan)
        t_2 = odeint(adExBack, start_points[strip+1], tspan)

        t_1 = np.array([t for t in t_1[:] if t[0] < -30.0 and t[1] < 300 and t[0] > -80.01 and t[1] > -800])
        t_2 = np.array([t for t in t_2[:] if t[0] < -30.0 and t[1] < 300 and t[0] > -80.01 and t[1] > -800])

        for cell in range(min(len(t_2[0:])-1,len(t_1[0:])-1)):
            p1 = [t_1[cell][0],t_1[cell][1]]
            p2 = [t_1[cell+1][0],t_1[cell+1][1]]
            p3 = [t_2[cell][0],t_2[cell][1]]
            p4 = [t_2[cell+1][0],t_2[cell+1][1]]

            if PolyArea([p1[0],p2[0],p3[0],p4[0]], [p1[1],p2[1],p3[1],p4[1]]) < 0.001:
                break

            if not isPolyConvex(np.array([p1,p2,p4,p3])):
                break

            lines = lines + [[p1,p3],[p4,p2]]

        ax.plot(t_1[:,0], t_1[:,1], color='k')
        ax.plot(t_2[:,0], t_2[:,1], color='k')

    lc = mc.LineCollection(lines, linewidths=2)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)

    ####################################################

    start_v = -40
    start_range_w = range(-700,300,20)

    start_points = [[start_v,w] for w in start_range_w]

    tspan = np.linspace(0, 50, 100)

    lines = []
    for strip in range(len(start_points)-1):
        t_1 = odeint(adExBack, start_points[strip], tspan)
        t_2 = odeint(adExBack, start_points[strip+1], tspan)

        t_1 = np.array([t for t in t_1[:] if t[0] < -30.0 and t[1] < 300 and t[0] > -80.01 and t[1] > -800])
        t_2 = np.array([t for t in t_2[:] if t[0] < -30.0 and t[1] < 300 and t[0] > -80.01 and t[1] > -800])

        for cell in range(min(len(t_2[0:])-1,len(t_1[0:])-1)):
            p1 = [t_1[cell][0],t_1[cell][1]]
            p2 = [t_1[cell+1][0],t_1[cell+1][1]]
            p3 = [t_2[cell][0],t_2[cell][1]]
            p4 = [t_2[cell+1][0],t_2[cell+1][1]]

            if PolyArea([p1[0],p2[0],p3[0],p4[0]], [p1[1],p2[1],p3[1],p4[1]]) < 0.001:
                break

            if not isPolyConvex(np.array([p1,p2,p4,p3])):
                break

            lines = lines + [[p1,p3],[p4,p2]]

        ax.plot(t_1[:,0], t_1[:,1], color='k')
        ax.plot(t_2[:,0], t_2[:,1], color='k')

    lc = mc.LineCollection(lines, linewidths=2)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)

    ################################################################

    plt.show()

generateMesh()
