import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

xmin = ymin = -0.1
xmax = ymax = 2.2
dx = 0.01
n = 4
x = np.arange(xmin,xmax,dx)
xnulx = 1/(0.5+x**n)
ynuly = 1/(0.5+x**n)

def f(xy,t):
    x, y = xy
    xp = 1/(0.5+y**n) - x
    yp = 1/(0.5+x**n) - y
    return [xp, yp]

dX = 0.12
dY = 0.12
X,Y = np.meshgrid(np.arange(xmin,xmax,dX),np.arange(xmin,xmax,dY))
dxdt = 1/(0.5+Y**n) - X
dydt = 1/(0.5+X**n) - Y

Tmax = 1000
dt = 0.1
xy0 = [1.2, 1.6]

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
ax1.plot(xnulx,x,lw=4, label="N-nullcline", c="C0")
ax1.plot(x,ynuly,lw=4, label="G-nullcline", c="C1")
ax1.plot(xy_out[:,0],xy_out[:,1],lw=3, label="trajectory", c="C5")
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
xy0 = [1.2, 1.6]

tvec = np.arange(0,Tmax,dt)
xy_out = odeint(f, xy0, tvec)
ax1.plot(xy_out[-1,0],xy_out[-1,1],'o',markerfacecolor=[0.9, 0.0, 0.9],
         markeredgecolor='black',markersize=15,markeredgewidth=2, label="PrE")

Tmax = 1000
dt = 0.1
xy0 = [1.6, 1.2]

tvec = np.arange(0,Tmax,dt)
xy_out = odeint(f, xy0, tvec)
ax1.plot(xy_out[-1,0],xy_out[-1,1],'o',markerfacecolor=[0.0, 0.8, 0.0],
         markeredgecolor='black',markersize=15,markeredgewidth=2, label="Epi")


Tmax = 1000
dt = 0.1
xy0 = [1.0, 1.0]

tvec = np.arange(0,Tmax,dt)
xy_out = odeint(f, xy0, tvec)
ax1.plot(xy_out[-1,0],xy_out[-1,1],'o',markerfacecolor=[1.0, 1.0, 1.0],
         markeredgecolor='black',markersize=15,markeredgewidth=2, label="Saddle")
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
savepath = os.path.join(save_dir, "bistability.pdf")
os.makedirs(save_dir, exist_ok=True)
plt.savefig(savepath, dpi=300, bbox_inches="tight")
plt.show()