from ovito.data import *
import ovito
import math
import numpy as np
import matplotlib.pyplot as plt


def modify(frame, input, output):
	xyz = input['Position'].array
	min_x = np.min(xyz[:,0])
	min_y = np.min(xyz[:,1])
	max_x = np.max(xyz[:,0])
	max_y = np.max(xyz[:,1])
	size = 2
	data=np.zeros((int((max_x-min_x)/size)+1,int((max_y-min_y)/size)+1),dtype=np.complex_)
	
	for i in xyz:
		data[int( (i[0]-min_x)/size ), int( (i[1]-min_y)/size)] += i[2] + 1j
			
	for ind1,i in enumerate(data):
		for ind2,jj in enumerate(i):
			if jj.imag==0:
				s=0
				for xi in [-1,0,1]:
					for yi in [-1,0,1]:
						try:
							data[ind1,ind2] += data[ind1+xi,ind2+yi]
							s += 1
						except:
						    pass
				data[ind1,ind2] = data[ind1,ind2]/s			
	
	data2=np.zeros(data.shape)
	for ind1,i in enumerate(data):
		for ind2,jj in enumerate(i):
			try:
				data2[ind1,ind2]=jj.real/jj.imag
			except:
				print(ind1,ind2)
	np.nan_to_num(data2)		
	np.savetxt('/users/ravinder/Desktop/ovito',data2)


	len_data = data2.shape[0]*data2.shape[1]
	
	xyz2 = np.ndarray((len_data,3))
	for ind1,i in enumerate(data2):
		for ind2,jj in enumerate(i):
			xyz2[ind2+data2.shape[0]*ind1,:]=[ind1*size+size/2,ind2*size+size/2,jj]
	
	
	pos_prop = ParticleProperty.create(ParticleProperty.Type.Position, len_data)
	for ind,i in enumerate(xyz2):
		pos_prop.marray[ind] = (i[0], i[1], i[2])	
		
	
	data2 = DataCollection()
	data2.add(pos_prop)
#	cell = SimulationCell()
#	cell.matrix = [[200,0,0,0],
#                       [0,200,0,-400],
#                       [0,0,20,-10]]
#	cell.pbc = (True, True, True)
#	cell.display.line_width = 0.1
#	data2.add(cell)
	
	node = ovito.ObjectNode()
	node.source = data2
	node.add_to_scene()
 
	