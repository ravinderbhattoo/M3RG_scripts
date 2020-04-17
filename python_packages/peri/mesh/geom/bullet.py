r"""
This is python package for Peridigm Mesh Generation(bullet).
"""


import math
import random
import numpy as np
def sample():
    return {'radius':7.6,'length':40,'res': 1 , 'id': 2,'center':[0,0,0],'rand':0.03}

def create(prop={}):

    r"""
    Parameters:
    ------------------
    diameter  : Diameter of bullet, float
    length    : length of bullet, float
    e_size    : Size of element, float
    ids       : Block id of bullet as in Peridigm Input file, int
    center    : center of bullet
    ------------------
    """
    prop_ = {'radius':7.6,'length':40,'res': 1 , 'id': 2,'center':[0,0,0],'rand':0.03}

    prop_.update(prop)
    radius = prop_['radius']
    length = prop_['length']
    e_size = prop_['res']
    ids = prop_['id']
    center = prop_['center']
    rand = prop_['rand']

    nm_radial_ele = int( (radius-e_size/2)/e_size + 0.5 )
    er_size = (radius-e_size/2) / nm_radial_ele

    nm_length = int( length/e_size + 0.5 )
    et_size = (length/nm_length)

    data = []
    total_vol = 0

    for i in range(nm_length):

        z=(i+0.5)*et_size + random.uniform(-rand,rand)*et_size
        total_vol += et_size*math.pi*e_size**2
        data.append([0+center[0],0+center[1],z+center[2],ids,et_size*math.pi*e_size**2])

        for j in range(nm_radial_ele):
            r_in = e_size/2 + (j) * er_size
            r_out = e_size/2 + (j+1) * er_size

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

    mask = []
    for i in data:
        if (i[0]-center[0])**2+(i[1]-center[1])**2 < (i[2]-center[2])**2 /(length/3)**2 * radius**2:
            mask.append(True)
        else:
            mask.append(False)

    data = np.array(data)
    data = data[mask]
    print('Volume:',total_vol,'(',length*math.pi*radius**2,')')

    return data

###################### END