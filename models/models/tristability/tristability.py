import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


xmin = ymin = -0.3
xmax = ymax = 2.5
dx = 0.001
x = np.arange(xmin,xmax,dx)

n = 4

an = 1
bn = 1
kn1 = 0.5
kn2 = 0.5
deltan = 1

ag = 1
bg = 1
kg1 = 0.5
kg2 = 0.5
deltag =  1

def f(xy,t=0):
    x, y = xy
    xp = (an*kn1**n)/(kn1**n+y**n) + (bn*x**n)/(kn2**n+x**n) - deltan*x
    yp = (an*kg1**n)/(kg1**n+x**n) + (bg*y**n)/(kg2**n+y**n) - deltag*y
    return [xp, yp]

dX = 0.12
dY = 0.12
X,Y = np.meshgrid(np.arange(xmin,xmax,dX),np.arange(xmin,xmax,dY))
dxdt = (an*kn1**n)/(kn1**n+Y**n) + (bn*X**n)/(kn2**n+X**n) - deltan*X
dydt = (an*kg1**n)/(kg1**n+X**n) + (bg*Y**n)/(kg2**n+Y**n) - deltag*Y

Tmax = 1000
dt = 0.1
xy0 = [0.3, 0.4]

tvec = np.arange(0,Tmax,dt)
xy_out = odeint(f, xy0, tvec)

import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
})
mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
mpl.rc('font', size=14) 
mpl.rc('axes', labelsize=14) 
mpl.rc('xtick', labelsize=14) 
mpl.rc('ytick', labelsize=14) 
mpl.rc('legend', fontsize=16) 

dxdt = 2.0*(dxdt / np.sqrt(dxdt**2 + dydt**2))
dydt = 2.0*(dydt / np.sqrt(dxdt**2 + dydt**2))

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(14,6))
ax1.quiver(X,Y,dxdt,dydt)
ax1.plot(xy0[0],xy0[1],'o',markerfacecolor=[0.0, 0.0, 0.0],
         markeredgecolor='black',markersize=15,markeredgewidth=2, label=r"[$N_0$, $G_0$]", zorder=10)
# To draw the nulclines
dX = dx
dY = dx
X,Y = np.meshgrid(np.arange(xmin,xmax,dX),np.arange(xmin,xmax,dY))
dxdt = (an*kn1**n)/(kn1**n+Y**n) + (bn*X**n)/(kn2**n+X**n) - deltan*X
dydt = (an*kg1**n)/(kg1**n+X**n) + (bg*Y**n)/(kg2**n+Y**n) - deltag*Y

cmap = mpl.cm.get_cmap('tab10')

cs1 = ax1.contour(X, Y, dxdt, levels=[0], linewidths=4, colors=[cmap(0)], label="N-nullcline")
cs2 = ax1.contour(X, Y, dydt, levels=[0], linewidths=4, colors=[cmap(1)], label="G-nullcline")

labels = ["N-nullcline"]
for i in range(len(labels)):
    cs1.collections[i].set_label(labels[i])

labels = ["G-nullcline"]
for i in range(len(labels)):
    cs2.collections[i].set_label(labels[i])

ax1.plot(xy_out[:,0],xy_out[:,1],lw=3, label="trajectory", c=cmap(5))
ax1.set_ylim(ymin,ymax)
ax1.set_xlim(xmin,xmax)
ax1.set_xlabel(r'$N$',fontsize=20)
ax1.set_ylabel(r'$G$',fontsize=20,rotation=0,labelpad=20)
ax1.tick_params(labelsize=20)

ax2.plot(tvec,xy_out[:,0],label=r'N',lw=3,color=[0.0, 0.9, 0.0])
ax2.plot(tvec,xy_out[:,1],label=r'G',lw=3,color=[0.9, 0.0, 0.9])
ax2.set_xlabel(r'$t$',fontsize=20)
ax2.set_ylabel(r'$N$, $G$',fontsize=20)
ax2.set_xlim(0,10)
ax2.legend()
ax2.tick_params(labelsize=20)

Tmax = 1000
dt = 0.1
xy0 = [0.2, 2.0]

tvec = np.arange(0,Tmax,dt)
xy_out = odeint(f, xy0, tvec)
ax1.plot(xy_out[-1,0],xy_out[-1,1],'o',markerfacecolor=[0.9, 0.0, 0.9],
         markeredgecolor='black',markersize=15,markeredgewidth=2, label="PrE")

Tmax = 1000
dt = 0.1
xy0 = [2.0, 0.2]

tvec = np.arange(0,Tmax,dt)
xy_out = odeint(f, xy0, tvec)
ax1.plot(xy_out[-1,0],xy_out[-1,1],'o',markerfacecolor=[0.0, 0.8, 0.0],
         markeredgecolor='black',markersize=15,markeredgewidth=2, label="Epi")

Tmax = 1000
dt = 0.1
xy0 = [1.0, 1.0]

tvec = np.arange(0,Tmax,dt)
xy_out = odeint(f, xy0, tvec)
ax1.plot(xy_out[-1,0],xy_out[-1,1],'o',markerfacecolor=[0.7, 0.7, 0.7],
         markeredgecolor='black',markersize=15,markeredgewidth=2, label="ICM")

saddle1 = [1.5, 0.49393833]
ax1.plot(saddle1[0],saddle1[1],'o',markerfacecolor=[1.0, 1.0, 1.0],
         markeredgecolor='black',markersize=15,markeredgewidth=2, label="Saddle")

saddle2 = [0.48825832, 1.506]
ax1.plot(saddle2[0],saddle2[1],'o',markerfacecolor=[1.0, 1.0, 1.0],
         markeredgecolor='black',markersize=15,markeredgewidth=2)

ax1.legend(facecolor="white", framealpha=1.0)
ax1.set_aspect('equal')
ax1.spines[['right', 'top']].set_visible(False)
ax2.spines[['right', 'top']].set_visible(False)

plt.tight_layout()

import os

cwd = os.getcwd()
foldername = os.path.basename(cwd)
basepath = os.path.dirname(os.path.dirname(cwd))
save_dir = os.path.join(basepath, "results", foldername)
savepath = os.path.join(save_dir, "tristability.pdf")
os.makedirs(save_dir, exist_ok=True)
plt.savefig(savepath, dpi=300, bbox_inches="tight")
plt.show()