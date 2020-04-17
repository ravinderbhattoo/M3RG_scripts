import numpy as np
from peri.simulation import session as sess

# Create session object
s1 = sess.run(name='test3')
s1.show_graph=False

prop1 =  {'Lx': 360,'Ly': 76, 'Lz': 1.9, 'res': 2, 'id': 1, 'center': [0, 0, 0]} # id = 1 for 2DSheet (matrix)
type1 = 'cuboid'

prop2 = {'radius': 5, 'thickness': 0.9, 'res': 1, 'id': 2, 'center': [-152, -43, 0], 'rand': 0.03}
type2 = 'disk'

prop3 = {'radius': 5, 'thickness': 0.9, 'res': 1, 'id': 3, 'center': [152, -43, 0], 'rand': 0.03}
type3 = 'disk'

prop4 = {'radius': 10, 'thickness': 0.9, 'res': 1, 'id': 4, 'center': [0, 48, 0], 'rand': 0.03}
type4 = 'disk'


# Define sub objects (notch)
sub_objs_inner = []
h=-2
k=2
a=-38
b=-19
logi_inner = '(x> {0}) * (x< {1}) * (y> {2}) * (y< {3})'.format(h,k,a,b)



# store all sub objects (notch) as python list.
sub_objs_inner.append({'logical':logi_inner,'remove':True})



# Create object as defined above with sub_objs (inclusions)
s1.add_object(name='Beam',type_=type1, prop=prop1, sub_objs=sub_objs_inner, SO_node_list=False)

s1.add_object(name='Roller1',type_=type2, prop=prop2, SO_node_list=False)

s1.add_object(name='Roller2',type_=type3, prop=prop3, SO_node_list=False)

s1.add_object(name='Roller3',type_=type4, prop=prop4, SO_node_list=False)

# Create data matrix
s1.create_data()
s1.assemble_data()



# Define node list using logical expression
s1.add_node_list(name='right_end', objects=[],logical='x<-152')
# right most edge (width 10)
s1.add_node_list(name='left_end', objects=[],logical='x>152') # left most edge (width 10)


# Write node list to files
s1.write_node_lists()


# Write matrix data to file (input data file for peridigm remove first two line if required)
s1.write_data()


# Write matrix data to file (id based on node list for easy visualization )
s1.write_data_node_list()
