from peri.simulation import session as sess
s1 = sess.run('test1')



prop =  {'Lx': 100, 'Ly': 100, 'Lz': 100, 'res': 2, 'id': 1, 'center': [0, 0, 0]}
type_ = 'cuboid'
s1.add_object(name='cuboid_1',type_=type_,prop=prop,apply={'rot':dict(vector=[0,1,0],point=[0,0,0],angle=10)})

sub_objs = [{'logical':'z>15','id':3},{'logical':'z>30','id':4}]

prop = {'radius':7.6,'length':40,'res': 1 , 'id': 2,'center':[0,0,70],'rand':0.03}
type_ = 'bullet'
s1.add_object(type_=type_,prop=prop,name='bullet',sub_objs=sub_objs)

#s1.create_data()
s1.rotate_object(name='cuboid_1',vector=[0,1,0],point=[0,0,0],angle=10)
s1.rotate_object(name='bullet',vector=[0,1,0],point='center',angle=-10)

s1.slide_object(name='bullet',to=[-100,0,0])
s1.assemble_data()

s1.add_node_list(objects=['cuboid_1'],logical='(x>45) | (y>45) | (x<-45) | (y<-45)')
s1.add_node_list(objects=['cuboid_1'],logical='')
s1.remove_node_list(name='node_list_2')
s1.write_node_lists()

s1.write_data()
s1.write_data_node_list()
