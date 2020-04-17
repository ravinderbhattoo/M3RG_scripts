r'''
Custom mesh generator for plate_bullet project
'''
import numpy as np
import random
import matplotlib.pyplot as plt

from peri.mesh.geom import bullet as mybullet
from peri.mesh.geom import cuboid

class plate:
    def __init__(self,parameters = dict(Lx = 0.10, Ly = 0.10, Lz = 0.02, res = 0.002, id = 1,center = (0,0,0))):
        #plate parameters
        self.parameters = {}
        self.exist = False
        self.parameters.update(parameters)

    def create(self):
        self.data = cuboid.create(self.parameters['Lx'],self.parameters['Ly'],self.parameters['Lz'],self.parameters['center'],self.parameters['res'],self.parameters['id'])
        self.exist = True

    def change(self,parameters):
        self.parameters.update(parameters)
        if self.exist:
            self.create()

    def get_data(self):
        if self.exist:
            pass
        else:
            self.create()
        return self.data

class bullet:
    def __init__(self,parameters = dict(dia = 0.01, length = 0.03, id = 2, res = 0.002, shape = 'none', center = (0,0,0))):
        #bullet parameters
        self.parameters = {}
        self.exist = False
        self.parameters.update(parameters)

    def create(self):
        self.data = mybullet.create(self.parameters['dia'],self.parameters['length'],self.parameters['res'],self.parameters['id'],self.parameters['center'],shape=self.parameters['shape'])
        self.exist = True

    def change(self,parameters):
        self.parameters.update(parameters)
        if self.exist:
            self.create()

    def get_data(self):
        if self.exist:
            pass
        else:
            self.create()
        return self.data

class mesh:
    def __init__(self):
        self.plates = {}
        self.bullets = {}
        self.node_lists = {}

    def add_plate(self,plate_id,parameters = {}):
        n = len(self.plates)
        if parameters == {}:
            self.plates[plate_id] = plate()
        else:
            self.plates[plate_id] = plate(parameters= parameters)
        self.show_plates()

    def remove_plate(self,plate_id):
        self.plates.pop(plate_id)
        self.show_plates()

    def change_plate(self,plate_id,parameters = {}):
        self.plates[plate_id].change(parameters)
        self.show_plates()

    def show_plates(self):
        if self.plates == {}:
            print('No plates')
        else:
            for p in self.plates:
                print('\n',p,':')
                for key in self.plates[p].parameters:
                    print('\t',key,': ',self.plates[p].parameters[key])




    def add_node_list(self,name,parameters = {}):
        n = len(self.plates)
        if parameters == {}:
            self.plates[plate_id] = plate()
        else:
            self.plates[plate_id] = plate(parameters= parameters)
        self.show_plates()

    def remove_node_list(self,plate_id):
        self.plates.pop(plate_id)
        self.show_plates()

    def change_node_list(self,plate_id,parameters = {}):
        self.plates[plate_id].change(parameters)
        self.show_plates()

    def show_node_list(self):
        if self.plates == {}:
            print('No plates')
        else:
            for p in self.plates:
                print('\n',p,':')
                for key in self.plates[p].parameters:
                    print('\t',key,': ',self.plates[p].parameters[key])


    def add_bullet(self,bullet_id,parameters = {}):
        n = len(self.bullets)
        if parameters == {}:
            self.bullets[bullet_id] = bullet()
        else:
            self.bullets[bullet_id] = bullet(parameters= parameters)
        self.show_bullets()

    def remove_bullet(self,bullet_id):
        self.bullets.pop(bullet_id)
        self.show_bullets()

    def change_bullet(self,bullet_id,parameters = {}):
        self.bullets[bullet_id].change(parameters)
        self.show_bullets()

    def show_bullets(self):
        if self.bullets == {}:
            print('No bullets')
        else:
            for p in self.bullets:
                print('\n',p,':')
                for key in self.bullets[p].parameters:
                    print('\t',key,': ',self.bullets[p].parameters[key])

    def create(self):
        data = []
        for p in self.plates:
            print(p)
            self.plates[p].create()
            data.append(self.plates[p].data)

        for b in self.bullets:
            print(b)
            self.bullets[b].create()
            data.append(self.bullets[b].data)

        self.data = np.concatenate(data,axis=0)

    def write(self,file_name):
        self.write_mesh(file_name)
        self.write_node_lists()

    def write_node_lists(self):
        for nl in self.node_lists:
            print('Writting node_list: ',nl)
            np.savetxt(nl,self.node_lists[nl],fmt=['%d'])


    def write_mesh(self,file_name):
        print('Total rows: ',len(self.data))
        np.savetxt(file_name,self.data,fmt=['%e','%e','%e','%d','%e'])


    def plot(self):
        import matplotlib
        matplotlib.rcParams.update({'font.size': 5})

        colors = [i for i in matplotlib.colors.cnames]
        random.shuffle(colors)
        color = [colors[int(i)] for i in self.data[:,3]]

        fig = plt.figure(figsize=(5, 3.5), dpi=200)
        fig.tight_layout()
        plt.subplot(221)
        plt.scatter(self.data[:,0],self.data[:,1],color=color)
        plt.gca().set_aspect('equal')
        plt.title('Top',fontsize=6)
        plt.subplot(222)
        plt.scatter(self.data[:,1],self.data[:,2],color=color)
        plt.gca().set_aspect('equal')
        plt.title('Right',fontsize=6)
        plt.subplot(223)
        plt.scatter(self.data[:,0],self.data[:,2],color=color)
        plt.gca().set_aspect('equal')
        plt.title('Front',fontsize=6)
        plt.subplot(224)
        plt.scatter(self.data[:,0],self.data[:,1],color=color)
        plt.gca().set_aspect('equal')
        plt.title('Top',fontsize=6)
        plt.show()
