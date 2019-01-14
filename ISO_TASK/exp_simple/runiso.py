# Once the Python Shared Library has been built in MIIND,
# copy this file to the results directory (where the .cpp and .so files were
# created).

import pylab
import numpy
import matplotlib.pyplot as plt
import imp
from operator import add
import libisopy as miind
from sklearn.decomposition import NMF,PCA
from sklearn.preprocessing import normalize

# Comment out MPI, comm and rank lines below if not using
# MPI
#######################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#######################

from scipy.signal import butter, lfilter

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

number_of_nodes = 1
simulation_length = 7 #s
miindmodel = miind.MiindModel(number_of_nodes, simulation_length)

miindmodel.init([])

timestep = miindmodel.getTimeStep()
print('Timestep from XML : {}'.format(timestep))

# For MPI child processes, startSimulation runs the full simulation loop
# and so will not return until MPI process 0 has completed. At that point,
# we want to kill the child processes so that sim() is not called more than once.
if miindmodel.startSimulation() > 0 :
    quit()

# each MN and associated INT gets the same afferent Ia input
# afferent input represents the degree of stretch of each muscle
# change this to change the "angle" of the leg
# 90 degrees
#prop_input = [3000,3000,3000,3000,3000,3000,3000]
#0 degrees
prop_input = [3000,6000,3000,3000,3000,6000,6000]
bg_input = [3000,3000,3000,3000,3000,3000,3000]
# each MN and assicuated INT gets the same supraspinal input
supra_input = [1000,1000,1000,1000,1000,1000,1000]
outputs = []
t = 0.0

for z in range(int(simulation_length/timestep)):
    t += timestep

    start_flexion = [1,1,1,1,1,1,1]
    # INTERESTING! : with a flexion ramp of < ~0.1 (100ms), the double peak in the second synergy
    # disappears
    flexion_ramp = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
    end_flexion = [5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5]
    # Differences in the maximum rate of each supraspinal input indicates which
    # which muscles are the agonists + variation among muscle activation
    rate = [8000,5000,8000,8000,8000,8000,8000]
    #rate = list( map(add, rate, prop_input) )

    for i in range(7):
        if(t > start_flexion[i]):
            if(t > end_flexion[i]):
                if(t < end_flexion[i]+flexion_ramp[i]):
                    supra_input[i] = (1.0-((t - end_flexion[i]) / (flexion_ramp[i]) ) ) * rate[i]
                else:
                    supra_input[i] = 0
            else:
                if(t < start_flexion[i]+flexion_ramp[i]):
                    supra_input[i] = ((t - start_flexion[i]) / (flexion_ramp[i]) ) * rate[i]
                else:
                    supra_input[i] = rate[i]
        else:
            supra_input[i] = 0

    node_input = list( map(add, supra_input, prop_input) )
    node_input = list( map(add, node_input, bg_input) )
    node_input += numpy.abs(numpy.random.normal(0,1,len(node_input))*500)
    #node_input = supra_input
    # Miind XML set up to only return the output of MNs
    o = miindmodel.evolveSingleStep(node_input.tolist())
    #o += o * (numpy.random.normal(0,1,len(o)) * 5 / 100)
    if(z > 100):
        outputs.append(o)

miindmodel.endSimulation()

res_list = numpy.matrix(outputs)
res_list = res_list.transpose()

# normalise values per muscle
for i in range(5):
#	res_list[i] += numpy.random.normal(25,1,res_list[i].shape)*1
#	res_list[i] = butter_bandpass_filter(res_list[i], 20, 450, 1000, 2)
	xmax, xmin = res_list[i].max(), res_list[i].min()
	res_list[i] = res_list[i]/xmax

#res_list = numpy.append(res_list, numpy.matrix([numpy.zeros(50),numpy.zeros(50),numpy.zeros(50),numpy.zeros(50),numpy.zeros(50)]), axis=1)

plt.figure()
plt.subplot(511)
plt.plot((res_list[0].tolist())[0])
plt.title("Firing Rates")
plt.subplot(512)
plt.plot((res_list[1].tolist())[0])
plt.subplot(513)
plt.plot((res_list[2].tolist())[0])
plt.subplot(514)
plt.plot((res_list[3].tolist())[0])
plt.subplot(515)
plt.plot((res_list[4].tolist())[0])

plt.show()

res_list = numpy.absolute(res_list)

comps = 2

model = NMF(n_components=comps, init='random', solver='cd', max_iter=2000, random_state=1234)

W = model.fit_transform(numpy.array(res_list))
H = model.components_

V = numpy.matmul(W,H)

for i in range(comps):
	plt.figure()
	plt.bar([1,2,3,4,5],W[:,i])
	plt.title("Synergy " + str(i))
	plt.show()


plt.figure()
plt.subplot((comps*100) + 11)
plt.plot(H[0,:])
plt.title("Firing Rates")
for i in range(comps-1):
	plt.subplot((comps*100) + 12+i)
	plt.plot(H[i+1,:])
plt.show()

# pca = PCA(n_components=2)
# pca.fit(res_list[2000:7000])
#
# print(pca.explained_variance_ratio_)
# print(pca.components_)
#
# plt.figure()
# plt.plot(np_rg_e.tolist()[2000:7000],np_rg_f.tolist()[2000:7000], '.')
# plt.title("time points")
# plt.show()
