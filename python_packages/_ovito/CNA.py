from ovito.data import *
import math	
import numpy as np
import ovito

#===========================================================================	
class add_filtered_bonds:
	def __init__(self, data):
		self.data = data
		self.bl = np.array([])
		self.array = data['Bonds'].array[:]
		self.pbc = data.bonds.pbc_vectors[:]
		xyz = data['Position'].marray
		for i in data.bonds.array:
			deltas = xyz[i[1]] - xyz[i[0]]
			self.bl = np.append(self.bl,math.sqrt(sum(deltas**2)))
		
	def filter(self, logic, *kwarg):
		bonds = Bonds()
		array = self.array[logic]
		pbc = self.pbc[logic]
		for ind,i in enumerate(array):
			bonds.add_full(i[0],i[1],pbc[ind])			
		bonds.display.width = 0.31
		bonds.display.color = (1,1,0)
		bonds.display.use_particle_colors = False
		self.bonds = bonds

#===========================================================================	
class add_filtered_atoms:
	def __init__(self, data, cutoff):
		self.data = data
		#self.input = input
		self.CN = np.array([])
		finder = CutoffNeighborFinder(cutoff, data)
		for index in range(data.number_of_particles):
			self.CN = np.append(self.CN,len([i for i in finder.find(index)]))
		
	def filter(self, logic, *kwarg):
		# Create the particle position property.
		self.pos_prop = ParticleProperty.create(ParticleProperty.Type.Position, sum(logic))
		start = 0
		xyz = self.data['Position'].marray
		for ind,i in enumerate(logic):
			if i:
				self.pos_prop.marray[start] = xyz[ind] 
				start += 1
		self.total_atoms = start 		
		# Create the particle type property and insert two atom types.
		self.type_prop = ParticleProperty.create(ParticleProperty.Type.ParticleType, sum(logic) )
		self.type_prop.type_list.append(ParticleType(id = 2, name = 'my_atom', color = (1.0,0.0,0.0)))
		start = 0
		for ind,i in enumerate(logic):
			if i:
				self.type_prop.marray[start] = 2
				start += 1
		
		
		
	
	
#===========================================================================	
#Start Here
#===========================================================================	
	
def modify(frame, input, output):

	my_bonds = add_filtered_bonds(output)
	my_bonds.filter(my_bonds.bl > 1.6)
	output.add(my_bonds.bonds)
	help(my_bonds.bonds)
	print('Total Bonds: ',(my_bonds.bonds.size)/2)
	my_atoms = add_filtered_atoms(output,1.92)
	my_atoms.filter(my_atoms.CN == 4)	
	output.add(my_atoms.type_prop)
	output.add(my_atoms.pos_prop)
	print('Total Atoms: ',my_atoms.total_atoms)
	
	
	