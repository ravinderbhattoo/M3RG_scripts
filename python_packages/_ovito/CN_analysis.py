from ovito.data import *
import numpy as np


def modify(frame, input, output):
	type = output['Particle Type'].array
	bonds = output['Bonds'].array
	
	CN = np.zeros((len(type),1))
	
	for i in bonds:
		CN[i[0]] += 1
	
	si_CN = CN[type==2]
	o_CN = CN[type==3]

	print('Average :',np.average(si_CN))
	print('4 :',sum(si_CN==4))
	print('3 :',sum(si_CN==3))
	print('2 :',sum(si_CN==2))

	print('Average :',np.average(o_CN))
	print('2 :',sum(o_CN==2))
	print('1 :',sum(o_CN==1))
	print('0 :',sum(o_CN==0))
	
	