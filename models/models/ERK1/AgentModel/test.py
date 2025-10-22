# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

times = np.loadtxt("/Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel/times.dat", comments="#", delimiter=",", unpack=False)
X = np.loadtxt("/Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel/X.dat", comments="#", delimiter=",", unpack=False)
Y = np.loadtxt("/Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel/Y.dat", comments="#", delimiter=",", unpack=False)
Z = np.loadtxt("/Users/pau/Desktop/PhD/EmbrionicModel/Nestor/code/AgentModel/Z.dat", comments="#", delimiter=",", unpack=False)

data = np.dstack([X,Y,Z])

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# Initialize scatters
scatters = [ ax.scatter(X[-1][i], Y[-1][i], Z[-1][i], s=1000) for i in range(data[0].shape[0]) ]

xmin = min(data[-1,:,0])
xmax = max(data[-1,:,0])
ymin = min(data[-1,:,1])
ymax = max(data[-1,:,1])
zmin = min(data[-1,:,2])
zmax = max(data[-1,:,2])
#Create 3d axes
#ax = plt.axes(xlim = (xmin,xmax),ylim=(ymin,ymax),zlim = (zmin,zmax),projection="3d")

# Setting the axes properties
ax.set_xlim3d([xmin, xmax])
ax.set_xlabel('X')

ax.set_ylim3d([ymin, ymax])
ax.set_ylabel('Y')

ax.set_zlim3d([zmin, zmax])
ax.set_zlabel('Z')

ax.grid(False)
