from ovito.data import *
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)     
fig.show()

Lx = 100
Ly = 100
N = 3680
NN = 499
T = 300

master_z = np.zeros((NN,N))
temp_z = 0
std_z = 0
mean_z = 0
counter = 0
do_it = False

n = np.array(range(0,18)) 
qx = 2*np.pi*n/100
qy = 1*qx
g_h2 = []
g_q = []

def clear():
	global g_h2
	g_h2 = []
	g_q = []


def rule(ax):

    KB = 8.6173303 * 10**(-5)
    S = Lx*Ly/N
    kappa = np.array([0.1,1.0,10,100])
    intercept = np.log10(N*KB*T/S/kappa)
    for ind,i in enumerate(intercept):
        x = np.array([0.2,-1.4])
        y = -4*x+i
        xy=(min(x[0],(-1+i)/4), max(1,y[0]))
        ax.annotate(str(kappa[ind]), xy=xy, xytext=xy)
        ax.plot(x,y)
        #ax.set_xlim([-1.4,0.2])
        #ax.set_ylim([1,8])


def DFT(xyz):
    KB = 8.6173303 * 10**(-5)
    S = Lx*Ly/N
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2] - mean_z
    global qx,qy
    h2 = []
    q = []
    for i in qx:
        for j in qy:
            total = np.abs(np.sum(z*(np.exp(-1j*(i*x+j*y)))))**2
            h2.append(total)
            q.append(np.sqrt(i**2+j**2))
    print('kappa = ', N*KB*T/S/(np.sum(np.array(q)**4 * np.array(h2))))
    print(np.mean(np.std(z)))
    print(h2)
    return h2,q
	
def my_prop(h,name,values):
    # Create a user-defined particle property.
    my_prop = ParticleProperty.create_user(name, 'float', len(values))
    my_prop.marray[:] = values
    h.add(my_prop)
	
def modify(frame, input, output):
    global temp_z,master_z,counter,mean_z,do_it,std_z,g_h2
    counter += 1
    xyz = output['Position'].marray
    master_z[frame] = xyz[:,2]
    if (counter>2) or (frame==NN):
        counter = 0
        ax.clear()
        temp_z = master_z[:frame,:].mean(axis=0)
        ax.hist(temp_z)
        if True:
            ax2.clear()
            ax3.clear()
            ax2.hist(std_z)
            ax3.scatter(mean_z,std_z)
        fig.canvas.draw()
        a,b = DFT(xyz)
        g_h2.append(a)
        ax4.clear()
        ax4.scatter(np.log10(b),np.log10(np.array(g_h2).mean(axis=0)))
        rule(ax4)
        fig.canvas.draw()
        	
    if (frame == NN)  or True:
        do_it = True
        mean_z = master_z.mean(axis=0)
        output['Position'].marray[:,2] -= mean_z
        std_z = master_z.std(axis=0)
        my_prop(output,'std_z',std_z)
        my_prop(output,'mean_z',mean_z)
    if frame==0:
        clear()
