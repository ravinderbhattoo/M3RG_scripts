from ovito.data import *
import numpy as np

def modify(frame, input, output):
	num = len(output['Position'].marray)
	radius = ParticleProperty.create(ParticleProperty.Type.Radius, num)
	radius.marray[:] *= 0
	radius.marray[:] += 0.00001	
	mask1 = (output['Damage'].marray <= 0.5)
	mask2 = (output['Particle Type'].marray==2)
	mask = mask1 & mask2
	radius.marray[mask] = np.ones((sum(mask),1))*0.001	
	output.add(radius)
