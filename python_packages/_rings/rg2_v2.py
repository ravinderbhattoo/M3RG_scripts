r"""
This is python package for 2D rings calculations.
yet to be completed
"""

import math
import os
import numpy as np
import pandas as pd
from lmps import frames as fr
import matplotlib.pyplot as plt


def get_rings(inpf,fn=[1],cut=6,bond_length=1.7,max_ring_size=6):
    os.makedirs(os.path.dirname(inpf), exist_ok=True)

    data_frames,box,names =  fr.get_frames(inpf,fn)
    data = data_frames[0]
    natoms = len(data)
    print('Number of atoms: ',natoms)

    lx_max = box[0][0][1] - box[0][0][0]
    ly_max = box[0][1][1] - box[0][1][0]
    lz_max = box[0][2][1] - box[0][2][0]

    if 'xs' in names:
        normalised = True
        print('Trajectory file is normalised.')
    else:
        normalised = False

    if normalised:
        data['x'] = data['xs']*lx_max
        data['y'] = data['ys']*ly_max
        data['z'] = data['zs']*lz_max


    print('Lx = ',lx_max)
    print('Ly = ',ly_max)

    #maintain periodicity
    data['x'] -= data['x'].min()
    data['y'] -= data['y'].min()

    df1 = data[data['x']<cut].copy()
    df1['id'] += 1.0e6
    df1['x'] += lx_max

    data = pd.concat([data,df1])

    df2 = data[data['y']<cut].copy()
    df2['id'] += 1.0e6
    df2['y'] += ly_max


    data = pd.concat([data,df2])

    df3 = data[(data['x']>(lx_max-cut)) & (data['x']<(lx_max))].copy()
    df3['id'] += 1.0e6
    df3['x'] -= lx_max

    data = pd.concat([data,df3])

    df4 = data[(data['y']>(ly_max-cut)) & (data['y']<(ly_max))].copy()
    df4['id'] += 1.0e6
    df4['y'] -= ly_max

    data = pd.concat([data,df4])


    print('New atoms in extra padding: ',len(data)-natoms)
    neighbour = {}
    for i in range(-1,int(lx_max/cut+1)+1):
        for j in range(-1,int(ly_max/cut+1)+1):
            mask_inner = (data['x']<=(i+1)*cut) & (data['x']>=i*cut) & (data['y']<=(j+1)*cut) & (data['y']>=j*cut) & (data['y']<=ly_max) & (data['x']<=lx_max) & (data['x']>=0) & (data['y']>=0)
            mask_outer = (data['x']<(i+2)*cut) & (data['x']>(i-1)*cut) & (data['y']<(j+2)*cut) & (data['y']>(j-1)*cut)

            iter_obj = data[mask_inner]
            comp_obj = data[mask_outer]

            for index, row in iter_obj.iterrows():
                r = np.sqrt((comp_obj['x']-row['x'])**2+(comp_obj['y']-row['y'])**2)
                mask = (r<bond_length) & (r>0.1*bond_length)
                neighbour[row['id']] = {}
                neighbour[row['id']]['id'] = comp_obj[mask]['id'].values
                neighbour[row['id']]['bond_length'] = r[mask]

    print('Calculating rings.')
    for i in neighbour:
        print('starting node: ',i)
        counter = 1
        trail = move_check(i,i,neighbour,counter,max_ring_size+3)
        print('trail: ',trail)
        break


def move_check(node,check,neighbour,counter,max_ring_size,offset='',parent=None):
    print(offset,'node: ',node)
    try:
        nb = neighbour[node]['id']
        all = []
        bool1 = False
        for i in nb:
            if (i == parent):
                print(offset,i,'parent')
                pass
            else:
                print('|',offset,i)
                if i == check:
                    print(offset,'its a ring!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    bool1 = True
                    all.append([i])
                else:
                    if counter>max_ring_size:
                        print(offset,'end')
                        pass
                    else:
                        a,b = move_check(i,check,neighbour,counter+1,max_ring_size,offset+'\t',parent=node)
                        print('else: ',a,b)
                        if b:
                            bool1 = b
                            print('adasdada')
                            for j in a:
                                all.append(j.append(i))
                        else:
                            pass

        print(all)
        return all,bool1
    except:
        return []

get_rings('./crystalline.dat')







def plot_rings(r_inpf,l,b,size=''):
    data = np.loadtxt(r_inpf,skiprows=1)
    data[:,2]=data[:,2]*float(l)
    data[:,3]=data[:,3]*float(b)

    rgs3 = np.load(r_inpf+'_rgs3.npy')

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if size =='':
        for ind,i in enumerate(rgs3):
            ax.plot([data[x,2] for x in i]+[data[i[0],2]]   ,[data[x,3] for x in i]+[data[i[0],3]],'black')
    else:
        for ind,i in enumerate(rgs3):
            if len(i)==size:
                ax.plot([data[x,2] for x in i]+[data[i[0],2]]   ,[data[x,3] for x in i]+[data[i[0],3]],'black')

    plt.axis('equal')
    plt.savefig(r_inpf+'.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)
    plt.show()


def plot(r_inpf):
    data = np.loadtxt(r_inpf+'_rings',skiprows=1)
    head = np.genfromtxt(r_inpf+'_rings',max_rows=1)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(data[:,0] ,data[:,1])
    plt.xlabel(head[0])
    plt.ylabel(head[1])
    plt.savefig(r_inpf+'_plot.png', dpi=600, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)
    plt.show()
