from ovito.data import *
import math	
import numpy as np
import ovito
#===========================================================================	
#Start Here
#===========================================================================	
	
def modify(frame, input, output):
	color_property = output.create_particle_property(ParticleProperty.Type.Color)
	for i in range(2,21):
		z_in = i*5.0
		z_fn = (i+1)*5.05
		mask = (output.position.array[:,2]<z_fn)*(output.position.array[:,2]>z_in )
		density = sum(mask)
		print(density)
		color_property.marray[mask] = (density, 0.0, 0.0)	
	
	print(output['Color'].array)