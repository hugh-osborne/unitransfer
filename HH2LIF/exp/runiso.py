# Once the Python Shared Library has been built in MIIND,
# copy this file to the results directory (where the .cpp and .so files were
# created).

import pylab
import numpy
import matplotlib.pyplot as plt
import imp
from sklearn.decomposition import NMF,PCA
from sklearn.preprocessing import normalize

# Comment out MPI, comm and rank lines below if not using
# MPI
#######################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#######################

miind = imp.load_dynamic('libiso', './libisopy.so')
number_of_nodes = 2
simulation_length = 7 #s
miindmodel = miind.MiindModel(number_of_nodes, simulation_length, 100,50,50,50)

miindmodel.init([])

timestep = miindmodel.getTimeStep()
print('Timestep from XML : {}'.format(timestep))

# For MPI child processes, startSimulation runs the full simulation loop
# and so will not return until MPI process 0 has completed. At that point,
# we want to kill the child processes so that sim() is not called more than once.
if miindmodel.startSimulation() > 0 :
    quit()

constant_input = [1000,1000]
rg_e = []
rg_f = []
t = 0.0

for i in range(int(simulation_length/timestep)):
    t += timestep

    start_flexion = 1
    # INTERESTING! : with a flexion ramp of < ~0.1 (100ms), the double peak in the second synergy
    # disappears  
    flexion_ramp = 0.2
    end_flexion = 5
    if(t > start_flexion):
        if(t > end_flexion):
            if(t < end_flexion+flexion_ramp):
                constant_input[0] = (1.0-((t - end_flexion) / (flexion_ramp) ) ) * 1000
            else:
                constant_input[0] = 0
        else:
            if(t < start_flexion+flexion_ramp):
                constant_input[0] = ((t - start_flexion) / (flexion_ramp) ) * 1000
            else:
                constant_input[0] = 1000
    else:
        constant_input[0] = 0

    if(t > start_flexion):
        if(t > end_flexion):
            if(t < end_flexion+flexion_ramp):
                constant_input[1] = (1.0-((t - end_flexion) / (flexion_ramp) ) ) * 1000
            else:
                constant_input[1] = 0
        else:
            if(t < start_flexion+flexion_ramp):
                constant_input[1] = ((t - start_flexion) / (flexion_ramp) ) * 1000
            else:
                constant_input[1] = 1000
    else:
        constant_input[1] = 0

    res = miindmodel.evolveSingleStep(constant_input)
    rg_e.append(res[0])
    rg_f.append(res[1])

    np_rg_e = numpy.array(rg_e)
    np_rg_f = numpy.array(rg_f)


res_list = numpy.array(list(zip(np_rg_e.tolist(),np_rg_f.tolist())))

for i in range(100):
    noise = numpy.random.normal(0, 0, np_rg_e.shape)
    np_rg_e = np_rg_e + noise
    np_rg_e[np_rg_e < 0] = 0

    noise = numpy.random.normal(0, 0, np_rg_f.shape)
    np_rg_f = np_rg_f + noise
    np_rg_f[np_rg_f < 0] = 0

    double_np_rg_e = numpy.array([list(np_rg_e)]).T
    double_np_rg_f = numpy.array([list(np_rg_f)]).T
    res_list = numpy.concatenate((res_list, double_np_rg_e), axis=1)
    res_list = numpy.concatenate((res_list, double_np_rg_f), axis=1)

# rg_e = normalize(np_rg_e.reshape(-1,1), axis=0, norm='l2').ravel()
# rg_f = normalize(np_rg_f.reshape(-1,1), axis=0, norm='l2').ravel()
#
# print(max(rg_e.tolist()))
# print(max(rg_f.tolist()))



# np_rg_e = np_rg_e / max(np_rg_e.tolist())
# np_rg_f = np_rg_f / max(np_rg_f.tolist())

model = NMF(n_components=2, init='nndsvd', random_state=0)

W = model.fit_transform(numpy.array(res_list))
H = model.components_

V = numpy.matmul(W,H)
print(W)
plt.figure()
plt.plot(W[:,0])
plt.title("Synergy 1")
plt.show()

plt.figure()
plt.plot(W[:,1])
plt.title("Synergy 2")
plt.show()

print(H[0:2,0:2])

pca = PCA(n_components=2)
pca.fit(res_list[2000:7000])

print(pca.explained_variance_ratio_)
print(pca.components_)

plt.figure()
plt.plot(np_rg_e.tolist()[2000:7000],np_rg_f.tolist()[2000:7000], '.')
plt.title("time points")
plt.show()

plt.figure()
plt.subplot(211)
plt.plot(np_rg_e)
plt.title("Firing Rates")
plt.subplot(212)
plt.plot(np_rg_f)

plt.show()
