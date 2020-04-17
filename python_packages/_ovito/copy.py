from ovito.data import *
import ovito
import math
import numpy as np
import matplotlib.pyplot as plt


def modify(frame, input, output):
	xyz = input['Position'].array
	data=np.zeros((100,100),dtype=np.complex_)
	for i in xyz:
		data[int(i[0]/2+50),int(i[1]/2+50)] += i[2] + 1j
		
	for k in range(0,3):
		for ind1,i in enumerate(data):
			for ind2,jj in enumerate(i):
				if jj.imag==0 and ind1>0 and ind1<99 and ind2>0 and ind2<99:
					try:
						data[ind1,ind2]=0.25*(data[ind1-1,ind2]+data[ind1+1,ind2]+
									data[ind1,ind2-1]+data[ind1,ind2+1])
					except:
						pass
	data2=np.ndarray((98,98))
	for ind1,i in enumerate(data[1:-2]):
		for ind2,jj in enumerate(i[1:-2]):
			data2[ind2,ind1]=jj.real/jj.imag
	np.nan_to_num(data2)		
	print(data2)
	np.savetxt('/users/ravinder/Desktop/ovito',data2)



	
	xyz2 = np.ndarray((98*98,3))
	for ind1,i in enumerate(data[1:-2]):
		for ind2,jj in enumerate(i[1:-2]):
			xyz2[ind2+98*ind1,:]=[ind1,ind2,jj.real/jj.imag]
	
	pos_prop = ParticleProperty.create(ParticleProperty.Type.Position, 98*98)
	for ind,i in enumerate(xyz2):
		pos_prop.marray[ind] = (2*i[0], 2*i[1]-400, i[2])	
		
	
	data2 = DataCollection()
	data2.add(pos_prop)
	cell = SimulationCell()
	cell.matrix = [[200,0,0,0],
               [0,200,0,-400],
               [0,0,20,-10]]
	cell.pbc = (True, True, True)
	cell.display.line_width = 0.1
	data2.add(cell)
	
	node = ovito.ObjectNode()
	node.source = data2
	node.add_to_scene()
 
	