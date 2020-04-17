r'''
Custom mesh generator for peridigm.


Example:

from peri.simulation import session as sess
s1 = sess.run()
prop =  {'Lx': 100, 'Ly': 100, 'Lz': 100, 'res': 2, 'id': 1, 'center': [0, 0, 0]}
type_ = 'cuboid'
s1.add_object(name='c1',type_=type_,prop=prop)
sub_objs = [{'logical':'Lz>15','id':3},{'logical':'Lz>30','id':4}]
prop = {'radius':7.6,'length':40,'res': 1 , 'id': 2,'center':[0,0,70],'rand':0.03}
type_ = 'bullet'
s1.add_object(type_=type_,prop=prop,name='b',sub_objs=sub_objs)

#s1.create_data()
s1.rotate_object(name='c1',vector=[0,1,0],point=[0,0,0],angle=10)
s1.rotate_object(name='b',vector=[0,1,0],point='center',angle=-10)
s1.slide_object(name='b',to=[0,0,0])
s1.assemble_data()
s1.add_node_list(objects=['c1'],logical='(x>45) | (y>45) | (x<-45) | (y<-45)')
s1.add_node_list(objects=['c1'],logical='')
s1.remove_node_list(name='node_list_2')
s1.write_node_lists()

s1.write_data()
s1.write_data_node_list()
s1.write_peri_input()

'''
from peri.logo import print_logo
import numpy as np
import random

print_logo()
print('\n\n\n\n')
def avail_samples():
    geom_samples = {}
    print('Avialable geometries: ')
    for name in geom_func:
        g_func = geom_func[name]
        print('\ttype_:',name)
        geom_samples[name] = g_func.sample()
        print('\tprop:',geom_samples[name])
    print('\n\n')
    return geom_samples


from copy import deepcopy
from peri.mesh.geom import *
geom_func = {}
for name in avail_geom.list_all():
    geom_func[name] = globals()[name]
geom_samples = avail_samples()


class run:
    r'''
    A class to run a session.

    '''
    def __init__(self,name = 'mesh'):
        r'''
        __init__ function for class.

        Input parameter:
        name	: Give a name to session.

        '''
        self.geom_samples = geom_samples
        self.data = None
        self.mask = None
        self.lcounter = 0
        self.assembled = False
        self.node_lists = {}
        self.obj_node_lists = {}
        self.sess_name  = name
        self.names = {}
        self.counter = 0
        self.valid_names = avail_geom.list_all()

    def graph(self):
        r'''
        Tree view of objects in session.

        Input parameter:

        '''
        if True:
            print('\n===============\nNode: ')
            count = 0
            print('\nObjects:')
            for obj in self.names:
                count += 1
                print('\t|---',str(count),': ',obj)
                for key in self.names[obj]:
                    if key=='data':
                        pass
                    else:
                        print('\t\t|---',key,' : ',self.names[obj][key])

            count = 0
            print('\nDefined Node Lists:')
            for lst in self.node_lists:
                count += 1
                print('\t|---',str(count),': ',lst)
                for key in self.node_lists[lst]:
                    if key=='data':
                        pass
                    else:
                        print('\t\t|---',key,' : ',self.node_lists[lst][key])

            count = 0
            print('\nObject Node Lists:')
            for lst in self.obj_node_lists:
                count += 1
                print('\t|---',str(count),': ',lst)
                for key in self.obj_node_lists[lst]:
                    if key=='data':
                        pass
                    else:
                        print('\t\t|---',key,' : ',self.obj_node_lists[lst][key])

            print('\n===============\n')

    def assemble_data(self):
        r'''
        Assemble data matrix.

        Input parameter:

        '''
        self.create_data()
        if self.assembled:
            print('Data already assembled.')
        else:
            print('Assembling data.....')
            self.data_list = [self.names[name]['data'] for name in self.names]
            self.data = np.vstack(self.data_list)
            self.data[:,:3] /= 1000
            self.data[:,4] /= 1.0e9
            self.assembled = True
            self.mask = self.data[:,0].copy()*0


    def add_node_list(self,name='none',logical='',objects=[]):
        r'''
        Add node lists to session.

        Input parameter:
        name      : Name of the node_list.

        '''
        if name == 'none':
            self.lcounter += 1
            name = 'node_list_' +str(self.lcounter)
        self.node_lists[name] = dict(logical=logical,objects=objects)
        self.graph

    def remove_node_list(self,name=''):
        r'''
        Remove node lists from session.

        Input parameter:
        name      : Name of the node_list.

        '''
        if name in self.node_lists:
            self.node_lists.pop(name)
            self.graph
        else:
            print(name,' does not exist.')
            self.graph

    def create_node_lists(self):
        if self.assembled:
            print('Creating node lists.....')
            x = self.data[:,0]*1000
            y = self.data[:,1]*1000
            z = self.data[:,2]*1000
            ids = self.data[:,3]
            #objects node lists
            for name in self.names:
                self.obj_node_lists[name] = {'info':'node_list_'+name}
                id = self.names[name]['prop']['id']
                mask = (ids == id)
                self.obj_node_lists[name]['data'] = np.array(range(1,len(self.data)+1))[mask]
                if self.names[name]['SO_node_list']:
                    for ind,sb_ob in enumerate(self.names[name]['Sub Objects']):
                        self.obj_node_lists[name+'_sb_obj_'+str(ind)] = {'info':'node_list_'+name+'_sb_obj_'+str(ind)}
                        id = sb_ob['id']
                        mask = (ids == id)
                        self.obj_node_lists[name+'_sb_obj_'+str(ind)]['data'] = np.array(range(1,len(self.data)+1))[mask]


            #defined node list
            count = 0
            for lst in self.node_lists:
                count += 1
                mask1 = True
                for obj in self.node_lists[lst]['objects']:
                    mask1 = (ids == self.names[obj]['prop']['id'])
                mask = eval(self.node_lists[lst]['logical'])
                mask = mask & mask1
                self.mask[mask] = count
                self.node_lists[lst]['data'] = np.array(range(1,len(self.data)+1))[mask]
            self.graph
            return True
        else:
            print('Node list can only be done after assembling of final data matrix!!!')
            return False

    def write_node_lists(self):
        if self.create_node_lists():
            print('Writing node list....')
            for lst in self.node_lists:
                np.savetxt('node_list_'+lst,self.node_lists[lst]['data'],fmt=['%d'])
            for lst in self.obj_node_lists:
                np.savetxt(self.obj_node_lists[lst]['info'],self.obj_node_lists[lst]['data'],fmt=['%d'])

    def write_data(self):
        r'''
        Write data files.

        Input parameter:

        '''
        self.assemble_data()
        print('Writting data.....')
        lines = len(self.data)
        filename = 'mesh_'+self.sess_name+'.data'
        np.savetxt(filename,self.data,fmt=['%e','%e','%e','%d','%e'],comments='',
        header="{}\n{}".format(lines,'Properties=pos:R:3:species:S:1:per_atom_volume:R:1'))
        print('Written to file: '+filename)

    def write_data_node_list(self):
        r'''
        Write data files.

        Input parameter:

        '''
        print('Re-writting data as node_lists.')
        self.assemble_data()
        print('Writting data.....')
        self.data[:,3] = self.mask*1
        lines = len(self.data)
        filename = 'mesh_' + self.sess_name+'_node_list'+'.data'
        np.savetxt(filename,self.data,fmt=['%e','%e','%e','%d','%e'],comments='',
        header="{}\n{}".format(lines,'Properties=pos:R:3:species:S:1:per_atom_volume:R:1'))
        print('Written to file: '+filename)

    def create_data(self):
        r'''
        Calculates data points for objects.

        Input parameter:

        '''
        for name in self.names:
            if self.names[name]['changed']:
                print('Creating data for: ', name+'.....')
                type_ = self.names[name]['type_']
                prop = self.names[name]['prop']
                apply = self.names[name]['apply']
                self.names[name]['data'] = np.array( geom_func[type_].create(prop=prop))
                if 'rot' in apply.keys():
                    self.rotate_object(name=name,**apply['rot'])
                if 'slide' in apply.keys():
                    self.slide_object(name=name,**apply['slide'])
                x = self.names[name]['data'][:,0] - prop['center'][0]
                y = self.names[name]['data'][:,1] - prop['center'][1]
                z = self.names[name]['data'][:,2] - prop['center'][2]
                ids = self.names[name]['data'][:,3]
                rmask = np.array([False for i in range(len(x))])
                for sb_ob1 in self.names[name]['Sub Objects']:
                    sb_ob = {'overlap':True,'remove':False}
                    sb_ob.update(sb_ob1)
                    mask1 = eval(sb_ob['logical'])
                    if sb_ob['remove']:
                        rmask = rmask + mask1
                        print(sum(rmask),' particles removed from ',name)
                    else:
                        mask2 = (ids==prop['id'])
                        if sb_ob['overlap']:
                            self.names[name]['data'][:,3][mask1] = sb_ob['id']
                        else:
                            if sum(mask1)==sum(mask1 & mask2):
                                self.names[name]['data'][:,3][mask1] = sb_ob['id']
                            else:
                                print('Overlap avoided. Logical exp: '+ sb_ob['logical'])
                self.names[name]['data'] = self.names[name]['data'][~rmask]
                self.names[name]['changed'] = False
                self.assembled = False

    def add_object(self,name = 'none', prop = {}, type_ = 'cuboid', sub_objs=[],SO_node_list=False,apply={}):
        r'''
        Adds object to session.

        Input parameter:
        name      : Give a name to object.
        prop      : Properties of the object.
        type_      : type_ of object.

        '''
        if type_ in self.valid_names:
            if name == 'none':
                self.counter += 1
                name = type_ + '_object_' +str(self.counter)
            self.names[name] = {'type_':type_,'prop':prop,'changed':True,'Sub Objects':sub_objs, 'SO_node_list':SO_node_list,'apply':apply}
            self.graph
        else:
            print(type_,' is not a valid geometry type_')
            print('Use one of these only:')
            print(self.valid_names)

    def make_copy(self,name):
        r'''

        '''
        for objs in list(self.names.keys()):
            if objs==name:
                self.names[name+'_copy'] = deepcopy(self.names[name])

    def remove_object(self,name = 'none'):
        r'''
        Remove object from session.

        Input parameter:
        name      : Name of the object.

        '''
        if name in self.names:
            self.names.pop(name)
            self.graph
        else:
            print(name,'(obj) does not exist.')
            self.graph


    def rotate_object(self,name='',vector=[],point=[],angle=0):
        if name in self.names:
            t = angle/180*np.pi
            obj = self.names[name]
            center = obj['prop']['center']
            if point == 'center':
                point = center
            try:
                data = obj['data'][:,:3]-point
            except:
                self.create_data()
                data = obj['data'][:,:3]-point
            vector = np.array(vector)
            unit_vector = vector/(np.sqrt(np.sum(vector**2)))
            l = unit_vector[0]
            m = unit_vector[1]
            n = unit_vector[2]
            c = np.cos(t)
            s = np.sin(t)
            rot = np.array([[l*l*(1-c)+c,m*l*(1-c)-n*s,n*l*(1-c)+m*s],
                            [l*m*(1-c)+n*s,m*m*(1-c)+c,n*m*(1-c)-l*s],
                            [l*n*(1-c)-m*s,m*n*(1-c)+l*s,n*n*(1-c)+c]])
            new_data = np.dot(rot,np.transpose(data))
            obj['data'][:,:3] = np.transpose(new_data) + point
            self.assembled = False
            print(name,'is rotated.')
        else:
            print('No such object in list.')
            print(self.names)


    def slide_object(self,name='none', by = True, to = True):
        if name in self.names:
            try:
                if to==True:
                    pass
                else:
                    self.names[name]['data'][:,:3] += (np.array(to) - np.array(self.names[name]['prop']['center']))
                if by==True:
                    pass
                else:
                    self.names[name]['data'][:,:3] += np.array(by)
                print(name,'is shifted.')

            except:
                self.create_data()
                if to==True:
                    pass
                else:
                    self.names[name]['data'][:,:3] += (np.array(to) - np.array(self.names[name]['prop']['center']))
                if by==True:
                    pass
                else:
                    self.names[name]['data'][:,:3] += np.array(by)
                print(name,'is shifted.')

        else:
            print('No such object in list.')
            print(self.names)




##############################################################################
