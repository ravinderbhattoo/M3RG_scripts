from ovito.data import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import *

fig,ax = plt.subplots(1,1)
fig.show()
ho = np.zeros((179,1))
std = np.zeros((179,1))
kr = np.zeros((179,1))

def modify(frame, input, output):
	if 1:
		pos = output['Position'].marray
		ho[frame] = output['v_totho'].marray.mean()	
		std[frame] = output['v_totho'].marray.std()	
		kr[frame] = kurtosis(output['v_totho'].marray)	
		
		ax.clear()		
		ax.plot(range(179),ho,color='r')
		ax.plot(frame,ho[frame],'r*')

		ax.plot(range(179),std,color='b')
		ax.plot(frame,std[frame],'bs')

		ax.plot(range(179),kr,color='g')
		ax.plot(frame,kr[frame],'go')

		fig.canvas.draw()
	if frame==178:
		np.savetxt('/Users/ravinder/Main/Graphene/Cooling_rate/op_12_100.txt',np.hstack([ho.reshape(-1,1),std.reshape(-1,1),kr.reshape(-1,1),]))