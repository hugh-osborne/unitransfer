# Once the Python Shared Library has been built in MIIND,
# copy this file to the results directory (where the .cpp and .so files were
# created).

import pylab
import numpy
import matplotlib.pyplot as plt
import imp

# Comment out MPI, comm and rank lines below if not using
# MPI
#######################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#######################

miind = imp.load_dynamic('libgrid', './libgrid.so')
number_of_nodes = 1 # get rid of this MPI nonsense.
simulation_length = 2 #s
miindmodel = miind.MiindModel(number_of_nodes, simulation_length, 0.1)

miindmodel.init([])

timestep = miindmodel.getTimeStep()
print('Timestep from XML : {}'.format(timestep))

# For MPI child processes, startSimulation runs the full simulation loop
# and so will not return until MPI process 0 has completed. At that point,
# we want to kill the child processes so that sim() is not called more than once.
if miindmodel.startSimulation() > 0 :
    quit()

constant_input = [15000,15000]
rg_e = []
t = 0.0

for i in range(int(simulation_length/timestep)):
    t += timestep

    constant_input[0] = 15000
    constant_input[1] = 0

    if (t > 0.3):
        constant_input[1] = 15000

    res = miindmodel.evolveSingleStep(constant_input)
    rg_e.append(res[0])

    np_rg_e = numpy.array(rg_e)

plt.figure()
plt.plot(np_rg_e)
plt.title("RG_E Firing Rate")
plt.show()
