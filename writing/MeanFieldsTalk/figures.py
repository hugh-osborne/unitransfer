import pylab
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import imageio
from scipy.stats import kde
from matplotlib import collections as mc

tau             = 0.01
V_threshold     = -35.0
epsilon         = 0.001
labda           = 0.0001
V_rest          = -57.5
V_min           = -70.0
V_max           = -34.99
N_grid          = 500
dt              = 0.0001
strip_w         = 0.005

#######################
###### LIF STRIP EXAMPLE
#
# ts = dt * np.arange(N_grid)
# pos_vs = V_rest + (V_threshold-V_rest)*np.exp(-ts/tau)
#
# f, (a0, a1) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1]})
# a0.set_ylim((V_rest,V_threshold))
# a0.plot(ts*400,pos_vs,'k')
# a0.set_xlabel('Time (dt)')
# a0.set_ylabel('Membrane Potential (mV)')
# a1.tick_params(axis='both',which='both',left=False,bottom=False,labelbottom=False,labelleft=False,labelsize=20)
# a1.set_yticks(pos_vs[0::20], minor=True)
# a1.yaxis.grid(b=True,which='minor')
# a1.set_ylim((V_rest,V_threshold))
# a1.set_xlim((0, 1))
# f.savefig('lif_strip_empty.svg', format='svg')
# f.savefig('lif_strip_empty.png',  dpi=250, format='png')
#
# images = []
# for anim in range(25):
#     t_ts_start = ts[anim*20]
#     t_ts = (dt * np.arange(21)) + t_ts_start
#     t_vs =  V_rest + (V_threshold-V_rest)*np.exp(-t_ts/tau)
#
#     f, (a0, a1) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1]})
#     a0.set_ylim((V_rest,V_threshold))
#     a0.plot(ts*400,pos_vs,'k')
#     a0.plot(t_ts*400,t_vs,'r',linewidth=4)
#     a0.set_xlabel('Time (dt)')
#     a0.set_ylabel('Membrane Potential (mV)')
#     a1.tick_params(axis='both',which='both',left=False,bottom=False,labelbottom=False,labelleft=False,labelsize=20)
#     a1.set_yticks(pos_vs[0::20], minor=True)
#     a1.yaxis.grid(b=True,which='minor')
#     a1.set_ylim((V_rest,V_threshold))
#     a1.set_xlim((0, 1))
#     a1.add_patch(
#     plt.Rectangle((0, t_vs[0]),1,t_vs[-1]-t_vs[0],color='r')
#     )
#     f.savefig('lif_gif/lif_strip_' + str(anim) + '.png', dpi=250, format='png')
#     images.append(imageio.imread('lif_gif/lif_strip_' + str(anim) + '.png'))
#
# imageio.mimsave('lif_strip.gif', images, duration=0.5)

####### LIF EXAMPLE
#######################

#######################
####### RANDOM NEURONS

# EMPTY STATE SPACE

# f = plt.figure()
# plt.xlim((V_rest,V_threshold))
# plt.ylim((0, 3))
# plt.xlabel('Time-dependent Variable 1')
# plt.ylabel('Time-dependent Variable 2')
# f.savefig('empty_state_space.svg', format='svg')
# f.savefig('empty_state_space.png',  dpi=250, format='png')
#
# f = plt.figure()
# points = np.random.randn(2,300)
# points[0] = (1.3 * points[0]) - 45
# points[1] = (0.3 * points[1]) + 1.8
# plt.xlim((V_rest,V_threshold))
# plt.ylim((0, 3))
# plt.yticks([])
# plt.xticks([])
# plt.xlabel('Time-dependent Variable 1')
# plt.ylabel('Time-dependent Variable 2')
# plt.scatter(points[0],points[1], marker='o', s=5)
# f.savefig('random_neurons.svg', format='svg')
# f.savefig('random_neurons.png',  dpi=250, format='png')
#
#
# f = plt.figure()
# nbins=300
# x = points[0]
# y = points[1]
# k = kde.gaussian_kde([x,y])
# xi, yi = np.mgrid[V_rest:V_threshold:nbins*1j, 0:3:nbins*1j]
# zi = k(np.vstack([xi.flatten(), yi.flatten()]))
#
# # Make the plot
# plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
# plt.xlabel('Time-dependent Variable 1')
# plt.ylabel('Time-dependent Variable 2')
# plt.yticks([])
# plt.xticks([])
# # plt.colorbar()
# f.savefig('random_neurons_prob.svg', format='svg')
# f.savefig('random_neurons_prob.png',  dpi=250, format='png')
#
# plt.show()
####### RANDOM NEURONS
#######################

#######################
####### IZH VECTOR FIELD

# X = np.arange(-86, -20, 5)
# Y = np.arange(-16, 21, 2)
# print(X.shape[0])
# V = np.zeros([Y.shape[0],0])
# W = np.zeros([Y.shape[0],0])
# for x in X.tolist():
#     v_row = []
#     w_row = []
#     for y in Y.tolist():
#         v_row = np.append(v_row, [(0.04*x**2 + 5*x + 140 - y + 10)], axis=0)
#         w_row = np.append(w_row, [(0.5 * (0.2*x - y))], axis=0)
#
#     v_row = np.reshape(v_row, (Y.shape[0],1))
#     w_row = np.reshape(w_row, (Y.shape[0],1))
#     V = np.append(V, v_row, axis=1)
#     W = np.append(W, w_row, axis=1)
#
# f, ax = plt.subplots()
# plt.xlim((-86,-20))
# plt.ylim((-16, 20.0001))
# q = ax.quiver(X, Y, V, W)
#
# f.savefig('izh_vector_no_nullclines.svg', format='svg')
# f.savefig('izh_vector_no_nullclines.png',  dpi=250, format='png')
#
# ## NULLCLINES
#
# v_null_ys = [0.04*x**2 + 5*x + 140 + 10 for x in X.tolist()]
# w_null_ys = [0.2*x for x in X.tolist()]
#
# ax.plot(X,v_null_ys,'r')
# ax.plot(X,w_null_ys,'b')
# plt.xlabel('Time-dependent Variable 1')
# plt.ylabel('Time-dependent Variable 2')
#
# f.savefig('izh_vector.svg', format='svg')
# f.savefig('izh_vector.png',  dpi=250, format='png')
#
# plt.show()

####### IZH VECTOR FIELD
#######################

#######################
####### IZH MESH

# images = []
# for anim in range(30):
#
#     starting_1_w = np.arange(-20, 30, 1)
#     starting_1_v = -86
#     dt = 0.1
#     t_series = np.arange(0, 10, dt)
#     quad_time = anim
#     quad_strip = 5
#
#     f, ax = plt.subplots()
#     plt.xlim((-86,-20))
#     plt.ylim((-16, 20))
#     plt.xlabel('Time-dependent Variable 1')
#     plt.ylabel('Time-dependent Variable 2')
#
#     strip_count = 0
#     for s in starting_1_w:
#         ps_w_1 = [s]
#         ps_v_1 = [starting_1_v]
#         ps_w_2 = [s+1]
#         ps_v_2 = [starting_1_v]
#         lines = [[(starting_1_v, s), (starting_1_v, s+1)]]
#         quad_points = []
#         time_count = 0
#         for t in t_series:
#             p1_v = ps_v_1[-1] + dt*(0.04*ps_v_1[-1]**2 + 5*ps_v_1[-1] + 140 - ps_w_1[-1] + 10)
#             p1_w = ps_w_1[-1] + dt*(0.2 * (0.2*ps_v_1[-1] - ps_w_1[-1]))
#             p2_v = ps_v_2[-1] + dt*(0.04*ps_v_2[-1]**2 + 5*ps_v_2[-1] + 140 - ps_w_2[-1] + 10)
#             p2_w = ps_w_2[-1] + dt*(0.2 * (0.2*ps_v_2[-1] - ps_w_2[-1]))
#             ps_v_1 = np.append(ps_v_1, p1_v)
#             ps_w_1 = np.append(ps_w_1, p1_w)
#             ps_v_2 = np.append(ps_v_2, p2_v)
#             ps_w_2 = np.append(ps_w_2, p2_w)
#             lines = lines + [[(p1_v, p1_w), (p2_v, p2_w)]]
#             if quad_strip == strip_count:
#                 if quad_time == time_count+1:
#                     quad_points = quad_points + [[p1_v, p1_w]] + [[p2_v, p2_w]]
#                 if quad_time == time_count:
#                     quad_points = quad_points + [[p2_v, p2_w]] + [[p1_v, p1_w]]
#             time_count = time_count + 1
#
#
#         lc = mc.LineCollection(lines)
#         ax.add_collection(lc)
#         if quad_strip == strip_count:
#             ax.add_patch(plt.Polygon(quad_points,color='r'))
#         ax.plot(ps_v_1,ps_w_1, 'k')
#         ax.plot(ps_v_2,ps_w_2, 'k')
#         strip_count = strip_count + 1
#
#     quad_time = 35 - anim
#     quad_strip = 14
#
#     starting_2_v = np.arange(-86, -30, 1)
#     v_null_ys = [0.04*x**2 + 5*x + 140 + 10 for x in starting_2_v]
#
#     strip_count = 0
#     for s in range(len(starting_2_v)-1):
#         ps_w_1 = [v_null_ys[s]]
#         ps_v_1 = [starting_2_v[s]]
#         ps_w_2 = [v_null_ys[s+1]]
#         ps_v_2 = [starting_2_v[s+1]]
#         lines = [[(starting_2_v[s], v_null_ys[s]), (starting_2_v[s+1], v_null_ys[s+1])]]
#         quad_points = []
#         time_count = 0
#         for t in t_series:
#             p1_v = ps_v_1[-1] - dt*(0.04*ps_v_1[-1]**2 + 5*ps_v_1[-1] + 140 - ps_w_1[-1] + 10)
#             p1_w = ps_w_1[-1] - dt*(0.2 * (0.2*ps_v_1[-1] - ps_w_1[-1]))
#             p2_v = ps_v_2[-1] - dt*(0.04*ps_v_2[-1]**2 + 5*ps_v_2[-1] + 140 - ps_w_2[-1] + 10)
#             p2_w = ps_w_2[-1] - dt*(0.2 * (0.2*ps_v_2[-1] - ps_w_2[-1]))
#             ps_v_1 = np.append(ps_v_1, p1_v)
#             ps_w_1 = np.append(ps_w_1, p1_w)
#             ps_v_2 = np.append(ps_v_2, p2_v)
#             ps_w_2 = np.append(ps_w_2, p2_w)
#             lines = lines + [[(p1_v, p1_w), (p2_v, p2_w)]]
#             if quad_strip == strip_count:
#                 if quad_time == time_count+1:
#                     quad_points = quad_points + [[p1_v, p1_w]] + [[p2_v, p2_w]]
#                 if quad_time == time_count:
#                     quad_points = quad_points + [[p2_v, p2_w]] + [[p1_v, p1_w]]
#             time_count = time_count + 1
#
#         lc = mc.LineCollection(lines)
#         ax.add_collection(lc)
#         if quad_strip == strip_count:
#             ax.add_patch(plt.Polygon(quad_points,color='r'))
#         ax.plot(ps_v_1,ps_w_1, 'k')
#         ax.plot(ps_v_2,ps_w_2, 'k')
#         strip_count = strip_count + 1
#
#     quad_time = anim
#     quad_strip = 40
#
#     starting_3_v = np.arange(-86, -30, 1)
#     v_null_ys = [0.04*x**2 + 5*x + 140 + 10 for x in starting_3_v]
#
#     strip_count = 0
#     for s in range(len(starting_3_v)-1):
#         ps_w_1 = [v_null_ys[s]]
#         ps_v_1 = [starting_3_v[s]]
#         ps_w_2 = [v_null_ys[s+1]]
#         ps_v_2 = [starting_3_v[s+1]]
#         lines = [[(starting_3_v[s], v_null_ys[s]), (starting_3_v[s+1], v_null_ys[s+1])]]
#         quad_points = []
#         time_count = 0
#         for t in t_series:
#             p1_v = ps_v_1[-1] + dt*(0.04*ps_v_1[-1]**2 + 5*ps_v_1[-1] + 140 - ps_w_1[-1] + 10)
#             p1_w = ps_w_1[-1] + dt*(0.2 * (0.2*ps_v_1[-1] - ps_w_1[-1]))
#             p2_v = ps_v_2[-1] + dt*(0.04*ps_v_2[-1]**2 + 5*ps_v_2[-1] + 140 - ps_w_2[-1] + 10)
#             p2_w = ps_w_2[-1] + dt*(0.2 * (0.2*ps_v_2[-1] - ps_w_2[-1]))
#             ps_v_1 = np.append(ps_v_1, p1_v)
#             ps_w_1 = np.append(ps_w_1, p1_w)
#             ps_v_2 = np.append(ps_v_2, p2_v)
#             ps_w_2 = np.append(ps_w_2, p2_w)
#             lines = lines + [[(p1_v, p1_w), (p2_v, p2_w)]]
#             if quad_strip == strip_count:
#                 if quad_time == time_count+1:
#                     quad_points = quad_points + [[p1_v, p1_w]] + [[p2_v, p2_w]]
#                 if quad_time == time_count:
#                     quad_points = quad_points + [[p2_v, p2_w]] + [[p1_v, p1_w]]
#             time_count = time_count + 1
#
#         lc = mc.LineCollection(lines)
#         ax.add_collection(lc)
#         if quad_strip == strip_count:
#             ax.add_patch(plt.Polygon(quad_points,color='r'))
#         ax.plot(ps_v_1,ps_w_1, 'k')
#         ax.plot(ps_v_2,ps_w_2, 'k')
#         strip_count = strip_count + 1
#
#     f.savefig('izh_gif/izh_mesh_' + str(anim) + '.png', dpi=300, format='png')
#     images.append(imageio.imread('izh_gif/izh_mesh_' + str(anim) + '.png'))
#
# imageio.mimsave('izh_mesh.gif', images, duration=0.5)

####### IZH MESH
#######################

#######################
####### IZH MESH QUAD TRANSITION

images=[]
for anim in range(30):
    starting_1_w = np.arange(-20, 30, 1)
    starting_1_v = -86
    dt = 0.1
    t_series = np.arange(0, 10, dt)
    quad_time = 0
    quad_strip = 5

    f, ax = plt.subplots()
    plt.xlim((-86,-20))
    plt.ylim((-16, 20))

    quad_time = 25
    quad_strip = 14

    starting_2_v = np.arange(-86, -30, 1)
    v_null_ys = [0.04*x**2 + 5*x + 140 + 10 for x in starting_2_v]

    strip_count = 0
    for s in range(len(starting_2_v)-1):
        ps_w_1 = [v_null_ys[s]]
        ps_v_1 = [starting_2_v[s]]
        ps_w_2 = [v_null_ys[s+1]]
        ps_v_2 = [starting_2_v[s+1]]
        lines = [[(starting_2_v[s], v_null_ys[s]), (starting_2_v[s+1], v_null_ys[s+1])]]
        quad_points = []
        quad_points_trans = []
        time_count = 0
        for t in t_series:
            p1_v = ps_v_1[-1] - dt*(0.04*ps_v_1[-1]**2 + 5*ps_v_1[-1] + 140 - ps_w_1[-1] + 10)
            p1_w = ps_w_1[-1] - dt*(0.2 * (0.2*ps_v_1[-1] - ps_w_1[-1]))
            p2_v = ps_v_2[-1] - dt*(0.04*ps_v_2[-1]**2 + 5*ps_v_2[-1] + 140 - ps_w_2[-1] + 10)
            p2_w = ps_w_2[-1] - dt*(0.2 * (0.2*ps_v_2[-1] - ps_w_2[-1]))
            ps_v_1 = np.append(ps_v_1, p1_v)
            ps_w_1 = np.append(ps_w_1, p1_w)
            ps_v_2 = np.append(ps_v_2, p2_v)
            ps_w_2 = np.append(ps_w_2, p2_w)
            lines = lines + [[(p1_v, p1_w), (p2_v, p2_w)]]
            if quad_strip == strip_count:
                if quad_time == time_count+1:
                    quad_points = quad_points + [[p1_v, p1_w]] + [[p2_v, p2_w]]
                    quad_points_trans = quad_points_trans + [[p1_v+1.5, p1_w]] + [[p2_v+1.5, p2_w]]
                if quad_time == time_count:
                    quad_points = quad_points + [[p2_v, p2_w]] + [[p1_v, p1_w]]
                    quad_points_trans = quad_points_trans + [[p2_v+1.5, p2_w]] + [[p1_v+1.5, p1_w]]
            time_count = time_count + 1

        lc = mc.LineCollection(lines)
        ax.add_collection(lc)
        if quad_strip == strip_count:
            aaas = [1.0,0.0,0.0,0.0,0.0,0.0,0.0]
            for a in range(1+(anim*2)):
                for l in range(6):
                    aaas[l] = aaas[l] - (aaas[l] * 0.05)
                    aaas[l+1] = aaas[l+1] + (aaas[l] * 0.05)
            for l in range(6):
                trans = [[v+(l*1.5),w] for [v,w] in quad_points]
                ax.add_patch(plt.Polygon(trans,facecolor='r', edgecolor='k', linewidth=2, alpha=aaas[l]))
        ax.plot(ps_v_1,ps_w_1, 'k')
        ax.plot(ps_v_2,ps_w_2, 'k')
        plt.xlim((-60.5,-51.5))
        plt.ylim((3, 5.1))
        plt.yticks([])
        plt.xticks([])
        strip_count = strip_count + 1

    quad_time = 0
    quad_strip = 40

    starting_3_v = np.arange(-86, -30, 1)
    v_null_ys = [0.04*x**2 + 5*x + 140 + 10 for x in starting_3_v]

    f.savefig('master_gif/master_' + str(anim) + '.png', dpi=300, format='png')
    images.append(imageio.imread('master_gif/master_' + str(anim) + '.png'))

imageio.mimsave('master_solve.gif', images, duration=0.5)

####### IZH MESH QUAD TRANSITION
#######################

#######################
####### GRID TRANSITION

# starting_1_v = np.arange(-86, -20, 1)
# starting_1_w = np.arange(-16, 20, 1)
#
# f, ax = plt.subplots()
# plt.xlim((-86,-20))
# plt.ylim((-16, 20))
# plt.xlabel('Time-dependent Variable 1')
# plt.ylabel('Time-dependent Variable 2')
#
# lines = []
# for s in starting_1_v:
#     lines = lines + [[(s, starting_1_w[0]), (s, starting_1_w[-1])]]
#
# for s in starting_1_w:
#     lines = lines + [[(starting_1_v[0], s), (starting_1_v[-1], s)]]
#
# quad_1 = [[-60,0],[-60,1],[-59,1],[-59,0]]
# quad_2 = [[-57.3,0],[-57.3,1],[-56.3,1],[-56.3,0]]
#
# lc = mc.LineCollection(lines)
# ax.add_collection(lc)
# ax.add_patch(plt.Polygon(quad_1,color='r'))
# ax.add_patch(plt.Polygon(quad_2,facecolor='r', edgecolor='k', linewidth=2, alpha=0.6))
#
# plt.show()

####### GRID TRANSITION
#######################
