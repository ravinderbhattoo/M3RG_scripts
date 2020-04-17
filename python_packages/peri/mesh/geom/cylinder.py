r"""
This is python package for Peridigm Mesh Generation(disk).

Parameters:
------------------
r1,r2     : Radii of cylinder, float
thickness : Thickness of disk, float
e_size    : Size of element, float
ids       : Block id of disk as in Peridigm Input file, int
center    : center of disk
------------------
"""
import numpy as np
import math
import random

def sample():
    return {'r1':10,'r2':20,'length':40,'res': 1 , 'id': 2,'center':[0,0,70],'rand':0.03}

def create(prop={}):

    prop_ =  {'r1':10,'r2':20,'length':40,'res': 1 , 'id': 2,'center':[0,0,0],'rand':0.03}

    prop_.update(prop)
    r1 = prop_['r1']
    r2 = prop_['r2']
    thickness = prop_['length']
    e_size = prop_['res']
    ids = prop_['id']
    center = prop_['center']
    rand = prop_['rand']


    r1=float(r1)
    r2=float(r2)
    thickness=float(thickness)
    e_size=float(e_size)
    ids=int(ids)

    nm_radial_ele = int( ( r2 - r1 )/e_size + 0.5 )
    er_size = ( r2 - r1 ) / nm_radial_ele

    nm_thickness = int( thickness/e_size + 0.5 )
    et_size = (thickness/nm_thickness)

    data = []
    total_vol = 0

    for i in range(nm_thickness):

        z=(i+0.5)*et_size + random.uniform(-rand,rand)*et_size
        for j in range(nm_radial_ele):
            r_in = r1 + (j) * er_size
            r_out = r1 + (j+1) * er_size

            nm_cir_ele = int(( math.pi * (r_in+r_out) ) / e_size + 0.5 )
            ec_size = ( math.pi * (r_in+r_out) ) / nm_cir_ele

            for k in range(nm_cir_ele):

                theta = 2*math.pi / nm_cir_ele
                e_area = math.pi* (r_out**2 - r_in**2) / nm_cir_ele
                r_g = (2 * math.pi * (r_out**3 - r_in**3) * math.sin(theta)/ 3 / theta ) / e_area / nm_cir_ele

                x = r_g * math.cos(k*theta) + random.uniform(-rand,rand)*e_size
                y = r_g * math.sin(k*theta) + random.uniform(-rand,rand)*e_size
                vol = e_area * et_size

                data.append([x+center[0],y+center[1],z+center[2],ids,vol])
                total_vol += vol

    print('Volume:',total_vol,'(',thickness*math.pi*(r2**2-r1**2),')')
    return np.array(data)

#################### END