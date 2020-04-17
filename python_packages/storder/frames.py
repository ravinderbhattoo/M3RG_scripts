r"""
    Get frames from lammps traj file.
    """
import pandas as pd
import numpy as np
import math

def getframe(inpf,fn):
    r"""
        Get specific frame from lammps traj file.

        Parameters
        ----------
        inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)
        fn     : frame number
        natoms : number of total atoms
        ----------
        """
    fn_copy = []
    frame = 1
    data_frames = []
    box = []
    names = []
    count = 0
    with open(inpf) as f:
        for line in f:
            if count < 3 or count == 4:
                count += 1
            elif count == 3:
                natoms = int(line.strip())
                count += 1
            elif count == 5:
                box.append([])
                box[-1].append([float(x) for x in line.strip().split(' ')])
                count += 1
            elif 5 < count < 8:
                box[-1].append([float(x) for x in line.strip().split(' ')])
                count += 1
            elif count == 8:
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
                else:
                    pass
                count = 0
                frame += 1

    return data_frames,box,names


def getframe_h(inpf):
    r"""
        Get specific frame from lammps traj file.

        Parameters
        ----------
        inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)
        ----------
        """

    import numpy as np
    head=np.genfromtxt(inpf,dtype='str',skip_header=8,max_rows=1)
    return head[2:]

def getframe_total(inpf,natoms):
    r"""
        Get specific frame from lammps traj file.

        Parameters
        ----------
        inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)
        natoms : Number of atoms, [int]
        ----------
        """
    import misc.file_len as fl
    return int((fl.len(inpf))/(9+natoms))

def box(inpf,fn,natoms):
    return np.genfromtxt(inpf,dtype='float',skip_header=5+(fn-1)*(natoms+9),max_rows=3)
