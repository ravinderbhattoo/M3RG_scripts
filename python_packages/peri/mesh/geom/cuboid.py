r"""
This is python package for Peridigm Mesh Generation(cuboid).

"""
import numpy as np

def sample():
    return {'Lx': 100, 'Ly': 100, 'Lz': 100, 'res': 1, 'id': 1, 'center': [0, 0, 0], 'rand': 0.03}

def create(prop={}):
    r"""
    Parameters:
    ------------------
    Lx      : Output data file (mesh) for Peridigm, /path_to_dir/filename
    radius    : Radius of cuboid, float
    thickness : Thickness of cuboid, float
    e_size    : Size of element, float
    ids       : Block id of cuboid as in Peridigm Input file, int
    center    : center of cuboid
    ------------------
    """
    prop_ = {'Lx': 100, 'Ly': 100, 'Lz': 100, 'res': 1, 'id': 1, 'center': [0, 0, 0], 'rand': 0.03}

    prop_.update(prop)

    Lx = prop_['Lx']
    Ly = prop_['Ly']
    Lz = prop_['Lz']
    center = prop_['center']
    res = prop_['res']
    id = prop_['id']
    rand = prop_['rand']

    size_x = int(Lx/res+1)
    res_x = Lx/size_x
    size_y = int(Ly/res+1)
    res_y = Ly/size_y
    size_z = int(Lz/res+1)
    res_z = Lz/size_z

    coord_x = np.arange(0.5*res_x,Lx,res_x)  #[(0.5+i)*res_x for i in range(size_x)]
    coord_y = np.arange(0.5*res_y,Ly,res_y)  #[(0.5+i)*res_y for i in range(size_y)]
    coord_z = np.arange(0.5*res_z,Lz,res_z)  #[(0.5+i)*res_z for i in range(size_z)]

    data = np.zeros((size_x*size_y*size_z,5))
    vol = res_x*res_y*res_z

    data[:,3] = id
    data[:,4] = vol
    total_vol = vol*len(data)

    for ind1,i in enumerate(coord_x):
        for ind2,j in enumerate(coord_y):
            for ind3,k in enumerate(coord_z):
                ind = ind1*size_y*size_z + ind2*size_z + ind3
                data[ind,0:3] = [i+center[0]-Lx/2, j+center[1]-Ly/2, k+center[2]-Lz/2]

    print('Volume:',total_vol,'(',Lx*Ly*Lz,')')
    data[:,:3]  = data[:,:3] + np.random.random((len(data),3))*rand*np.array([res_x,res_y,res_z])

    return data

#################### END
