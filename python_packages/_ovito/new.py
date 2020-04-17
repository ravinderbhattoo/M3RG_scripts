from ovito.data import *

def modify(frame, input, output):
	vol = output['Atomic Volume'].array[:15240]
	stress = output['c_peratom1'].array[:15240]
	avg = (sum(stress/vol)/15240)
	print(sum(vol)/15240)	
	print(avg/15240/(6.023*10**23)*4184*10**(30) / 10**9, ' GPa')
	
	 
