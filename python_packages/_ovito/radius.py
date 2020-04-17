from ovito.data import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from matplotlib import cm

import math
pi = math.pi
sin = math.sin
cos = math.cos


def createSphere(r, N=10):
    lst = []
    thetas = [(2*pi*i)/N for i in range(N)]
    phis = [(pi*i)/N for i in range(N)]
    for theta in thetas:
        for phi in phis:
            x = r * sin(phi) * cos(theta)
            y = r * sin(phi) * sin(theta)
            z = r * cos(phi)
            lst.append([2,x, y, z])
    return lst


def fit_sphere(x, r2,x1,y1,z1):
	return (r2-(x[0]-x1)**2-(x[1]-z1)**2)-y1**2+2*y1*x[2]
		
		
def modify(frame, input, output):
	xyz = output['Position'].marray
	xlim = [-90,50]
	zlim = [-150,-10]
	step = 10
	Y = []
	X = []
	xx = []
	zz = []
	for i in range(0,int((xlim[1]-xlim[0])/step)):
		X.append(xlim[0]+(i+0.5)*step)
		Z = []
		for j in range(0,int((zlim[1]-zlim[0])/step)):
			Z.append(zlim[0]+(j+0.5)*step)
			xx.append(xlim[0]+(i+0.5)*step)
			zz.append(zlim[0]+(j+0.5)*step)
			
			mask =  (xyz[:,0]>(xlim[0]+i*step)) & (xyz[:,0]<(xlim[0]+(i+1)*step))
			mask1 =  (xyz[:,2]>(zlim[0]+j*step)) & (xyz[:,2]<(zlim[0]+(j+1)*step))
			Y.append(xyz[mask & mask1,1].mean())
	
	y = np.array(Y)
	x = np.array(X)
	z = np.array(Z)
	
	
	print(len(x),len(z),len(x)*len(z),len(y),)
	X, Z = np.meshgrid(x, z)
	Y = y.reshape(X.shape)

	data = np.zeros((len(x)*len(z),4))
	for i in range(0,len(x)*len(z)):
		data[i,:]=[1,xx[i],y[i],zz[i]]
	np.savetxt('Users/ravinder/Desktop/mem_data+'+str(frame)+'.data',data,comments='',header=str(len(data))+'\n')

	mask = ~np.isnan(y)
	xx = np.array(xx)[mask]
	zz = np.array(zz)[mask]
	
	popt, pcov = curve_fit(fit_sphere,[xx,zz,y[mask]],y[mask]**2,p0=[999999999999,1,1,1])
	
	print(popt,pcov)
	
	data = np.array(createSphere(np.sqrt(popt[0]), N=30))
	data[:,1:] += np.array(popt[1:])
	
#	Y = -np.sqrt(popt[0]-(X-popt[1])**2-(Z-popt[3])**2)+popt[2]	

#	a = -np.sqrt(popt[0]-(xx-popt[1])**2-(zz-popt[3])**2)
#	y = a + popt[2]	

#	data = np.zeros((2*sum(mask),4))

#	for i in range(0,sum(mask)):
#		data[2*i,:]=[2,xx[i],y[i],zz[i]]
#		data[2*i+1,:]=[2,xx[i],-a[i]+popt[2],zz[i]]

	print(len(data))	
	np.savetxt('Users/ravinder/Desktop/fit_data+'+str(frame)+'.data',data,comments='',header=str(len(data))+'\n')
	print(np.sqrt(popt[0]))
	
	plt.scatter(frame,1/np.sqrt(popt[0]))
	plt.draw()
	plt.show()	
	
	
		
	
	
	