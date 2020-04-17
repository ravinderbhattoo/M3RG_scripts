from peri.simulation import session as sess
s1 = sess.run(name='test2')

prop =  {'Lx': 200, 'Ly': 200, 'Lz': 5, 'res': 1, 'id': 1, 'center': [0, 0, 0],"rand":0.1}
type_ = 'cuboid'

# Define sub objects (notch)
logi_inner1 = '(x> {0}) + (x< {1}) + (y> {2}) + (y< {3})'.format(50,-50,50,-50)

logi_inner2 = '(x> {0}) * (x< {1}) * (y> {2}) * (y< {3})'.format(-50,-40,-5,5)



# store all sub objects (notch) as python list.
sub_objs_inner = [{'logical':logi_inner1,'remove':True}, {'logical':logi_inner2,'remove':True}]

s1.add_object(name='bar',type_=type_,prop=prop,apply={'rot':dict(vector=[0,0,1],point=[0,0,0],angle=49)},sub_objs=sub_objs_inner)



s1.create_data()
s1.assemble_data()
s1.write_data()
s1.write_data_node_list()
