import ovito
from ovito.vis import * #Allows us to manipulate camera
from ovito import dataset #Tells Ovito to use whatever dataset is currently open
from math import * #MUST include this to access pi, sin, and cos
import numpy as np

class animation:
    def __init__(self,pos=(0,0,100),dir=(0,0,-1),fov=35,filepath=None):
        #Viewport object is from ovito.vis - This is what we will use to directly manipulate camera
        self.vp = Viewport()
        self.vp.type = Viewport.Type.PERSPECTIVE
        #Camera setup info found in "Adjust View" dialogue box
        self.vp.camera_dir = dir
        self.vp.camera_pos = pos
        self.vp.fov = fov*pi/180
        self.filepath = filepath
        self.frame = 0

    def translation(self,  d_x=0,d_y=0,d_z=0,steps=10,with_focus=None,trajectory=None,rend=True):
        pos = [i for i in self.vp.camera_pos]
        for i in range(steps+1):
            if trajectory==None:
                self.vp.camera_pos = (pos[0]+i*d_x/steps, pos[1]+i*d_y/steps, pos[2]+i*d_z/steps)
            else:
                self.vp.camera_pos = trajectory[i]

            if with_focus==None:
                pass
            else:
                dir0 = np.array((with_focus)) - np.array((self.vp.camera_pos))
                dir0 /= (dir0**2).sum()
                self.vp.camera_dir = (dir0[0],dir0[1],dir0[2])

            if rend:
                self.render()

    def zoom(self,s,steps=10,rend=True):
        for i in range(steps):
            self.vp.fov -= s*pi/180
            if rend:
                self.render()

    def rotation(self,theta,steps=36,rend=True):
        dir0 = np.array([i for i in self.vp.camera_dir])
        t = angle(dir0[1],dir0[0])
        mag = sqrt(dir0[0]**2+dir0[1]**2)
        #about z-axis
        for i in range(steps+1):
            new_t = t + i*theta/steps
            dir1 = dir0.copy()
            dir1[:2] = [mag*cos(new_t),mag*sin(new_t)]
            self.vp.camera_dir = (dir1[0],dir1[1],dir1[2])
            if rend:
                self.render()

    def render(self):
        filename = self.filepath + 'image_'+str(self.frame)+'.png'
        rs = RenderSettings(size=(400,400), renderer = TachyonRenderer(),filename = filename)
        self.vp.render(rs)
        self.frame += 1

    def combination(self,d_x,d_y,d_z,s_zoom,steps=10):
        for i in range(steps):
            self.translation(d_x/steps,d_y/steps,d_z/steps,steps=2,rend=False)
            self.zoom(self,s_zoom/steps,steps=1,rend=False)
            self.render()

filepath = '/Users/ravinder/Desktop/Ovito_movies/'
sess1 = animation(filepath=filepath,pos = (52,18,22), dir = (-1,0,0),fov=35)

R = 150
trajectory = [(-49+R*cos(i*6*np.pi/180),-49+R*sin(i*6*np.pi/180),70) for i in range(0,60)]

#sess1.translation(steps=len(trajectory)-1,trajectory=trajectory,with_focus=[-49,-49,0])
sess1.translation(steps=40,d_x = -50)
#sess1.zoom(-0.2,steps=50)
#sess1.rotation(2*np.pi)
