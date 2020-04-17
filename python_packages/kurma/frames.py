"""
Get frames from lammps trajectory file.
Dependency:
numpy
pandas
math
"""
import pandas as pd
import numpy as np
import math

def get_frames(inpf,fn=[1],variable_size=False):
    """
    Get specific frame from lammps traj file.

    Parameters
    inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)
    fn     : List of frame number
    natoms : Number of total atoms

    """

    if fn=='all':
            fn = range(1,get_frame_total(inpf)+1)

    fn_copy = []
    frame = 1
    data_frames = []
    box = []
    names = []
    count = 0
    size_of_frame = get_natoms(inpf) + 9
    max_frame = max(fn)
    print('Reading file: ',inpf)

    with open(inpf) as f:
        for line in f:
            if (frame in fn) or variable_size:
                if count < 3 or count == 4:
                    count += 1
                elif count == 3:
                    natoms = int(line.strip())
                    count += 1
                elif count == 5:
                    if frame in fn:
                        box.append([])
                        box[-1].append([float(x) for x in line.strip().split(' ')])
                    count += 1
                elif 5 < count < 8:
                    if frame in fn:
                        box[-1].append([float(x) for x in line.strip().split(' ')])
                    count += 1
                elif count == 8:
                    if frame in fn:
                        box[-1] = np.array(box[-1])
                        names =  [x for x in line.strip().split(' ')][2:]
                        df = {}
                        for n in names:
                            df[n] = []
                    count += 1
                elif 8 < count < 9+natoms-1:
                    if frame in fn:
                        v_list = [float(x) for x in line.strip().split(' ')]
                        for ind,n in enumerate(names):
                            df[n].append(v_list[ind])
                    else:
                        pass
                    count +=1
                else:
                    if frame in fn:
                        print('frame: ',frame)
                        v_list = [float(x) for x in line.strip().split(' ')]
                        for ind,n in enumerate(names):
                            df[n].append(v_list[ind])
                        fn_copy.append(frame)
                        if 'id' in names:
                            data_frames.append(pd.DataFrame(df).sort_values('id'))
                        else:
                            data_frames.append(pd.DataFrame(df))
                        if sorted(fn_copy) == sorted(fn):
                            break
                        if 'type' in data_frames[-1].columns:
                            pass
                        else:
                            data_frames[-1]['type'] = data_frames[-1]['id']*0+1 
                    else:
                        pass
                    count = 0
                    frame += 1
                    if frame>max_frame:
                        break

            else:
                count += 1
                if count == size_of_frame:
                    frame += 1
                    if frame>max_frame:
                        break
                    count = 0

    return data_frames,box,names


def get_frame_headers(inpf):
    """
    Get headers from lammps traj file.

    Parameters
    ----------
    inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)
    """
    import numpy as np
    head=np.genfromtxt(inpf,dtype='str',skip_header=8,max_rows=1)
    return head[2:]

def get_frame_total(inpf):
    """
    Get total number of frames in lammps traj file.

    Parameters
    ----------
    inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)

    """

    def filelen(fname):
        r"""
        Get length of file.
        """
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    natoms = get_natoms(inpf)
    return int((filelen(inpf))/(9+natoms))

def get_box(inpf,fn,natoms):
    """
    Get box from lammps traj file.

    Parameters
    ----------
    inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)
    fn     : Frame number, [int]
    natoms : Number of atoms, [int]

    """

    return np.genfromtxt(inpf,dtype='float',skip_header=5+(fn-1)*(natoms+9),max_rows=3)

def get_natoms(inpf):
    """
    Get natoms from lammps traj file.

    Parameters
    ----------
    inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)
    """
    return np.genfromtxt(inpf,dtype='int',skip_header=3,max_rows=1)
