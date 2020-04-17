r"""
This is python package for 2D rings calculations.
yet to be completed
"""

import math
import os
import numpy as np
global anb
anb=[]
def get_rings(inpf,l,b,natom,cut,CN_cut,size_ring):
    global anb
    os.makedirs(os.path.dirname(inpf), exist_ok=True)

    #load data
    data2=np.loadtxt(inpf,skiprows=1)

    #scalling of coordinates i.e xs->x
    data=data2
    data[:,2]=data[:,2]*l
    data[:,3]=data[:,3]*b

    if max(abs(data[:,2]))>l:
        l=max(data[:,2])-min(data[:,2])
        b=max(data[:,3])-min(data[:,3])

    #Number of cells in all region(x and y)
    n=math.floor(l/cut+0.000001)
    m=math.floor(b/cut+0.000001)
    cuts=CN_cut**2

    print(n,m)
    #size of each cell(x and y)
    lbin=l/n
    bbin=b/m

    #extra cells around system for periodicity
    n=n+2
    m=m+2

    #neighbourlist for each cell

    nb=[]
    for i in range(0,m):
        for j in range(0,n):
            if j>0 and j<n-1 and i>0 and i<m-1:
                nb.append([(i-1)*n+j-1,(i-1)*n+j,(i-1)*n+j+1,i*n+j-1,i*n+j,i*n+j+1,(i+1)*n+j-1,(i+1)*n+j,(i+1)*n+j+1])
            else:
                nb.append(['edge',-1*(i==0)+1*(i==m-1),-1*(j==0)+1*(j==n-1)])

    #Empty cells as list
    cell=[]
    for i in range(0,n*m):
        cell.append([])
    print('Total Cells : ',len(cell))

    #Outer cell to leave
    leaveit=[]
    leaveit2=[]
    for x in range(0,n):
        leaveit.append(x)
        leaveit.append(n*(m-1)+x)
        leaveit2.append(n+x)
        leaveit2.append(n*(m-2)+x)
    for x in range(0,m):
        leaveit.append(n*x)
        leaveit.append(n*x+n-1)
        leaveit2.append(n*x+1)
        leaveit2.append(n*x+n-2)
    #insert index of each atom in every cell
    for atom_index,atom_coord in enumerate(data):
        k=(math.floor(atom_coord[3]/bbin)+1)*n + math.floor(atom_coord[2]/lbin+1)
        cell[k].append(atom_index)

    #filling extracells
    cell[0]=cell[(m-2)*n+n-1]
    cell[n-1]=cell[(m-2)*n+1]
    cell[(m-1)*n]=cell[2*n-2]
    cell[(m-1)*n+n-1]=cell[n+1]
    for i in range(1,n-1):
        cell[i]=cell[(m-2)*n+i]
        cell[(m-1)*n+i]=cell[n+i]
    for j in range(1,m-1):
        cell[j*n]=cell[(j+1)*n-2]
        cell[(j+1)*n-1]=cell[j*n+1]




    #initialise zeros array for every atom (to insert no. of atoms at r distance)
    anb=[]
    snb=[]
    for i in data:
        anb.append([])
        snb.append([])

    #nested loop over all atoms
    for ind,i in enumerate(cell):

        ids1=[]
        ids2=[]
        if (ind in leaveit):
            print('Cell no. : ',ind,' has been left')
        else:
            print('Cell no. : ',ind)

            for ns in nb[ind]:
                if (nb[ns][0]=='edge'):
                    ids1.append(ns)
                else:
                    ids2.append(cell[ns])

            for ind2,atom_ind in enumerate(i):
                x0=data[atom_ind,2]
                y0=data[atom_ind,3]
                for k in ids2:
                    x1=data[k,2]
                    y1=data[k,3]
                    dx=x1-x0
                    dy=y1-y0
                    theta=np.arctan2(dy,dx)

                    theta[np.logical_and(theta<0,1)] += 2*math.pi
                    r=dx**2+dy**2
                    for inx,x in enumerate(r):
                        if x<cuts and x>0.1:
                            anb[atom_ind].append(k[inx])
                            snb[atom_ind].append(theta[inx])

            for ind2,atom_ind in enumerate(i):
                x0=data[atom_ind,2]
                y0=data[atom_ind,3]
                for ns in ids1:
                    for k in cell[ns]:
                        x1=data[k,2]+nb[ns][2]*l
                        y1=data[k,3]+nb[ns][1]*b
                        dx=x1-x0
                        dy=y1-y0
                        theta=math.atan2(dy,dx)
                        if theta<0:
                            theta += 2*math.pi
                        r=dx**2+dy**2
                        if r<cuts and x>0.1:
                            anb[atom_ind].append(k)
                            snb[atom_ind].append(theta)



    for ind,i in enumerate(anb[0:]):
       if i==[]:
           pass
       else:
           anb[ind]=[list(x) for x in zip(*sorted(zip(snb[ind], anb[ind]), key=lambda pair: pair[0]))][1]

   #ring calculations
    cy=int(size_ring)
    rings=[]
    rings2=[]
    getr=[]
    for ind,i in enumerate(cell):
        if (ind in leaveit) or (ind in leaveit2):
            pass
        else:
            for a_id in i:

                nbs=anb[a_id]
                rg=1*recurse(nbs,a_id,cy,a_id,[-1])
                getr.append([len(rings),a_id])
                #refine rings

                cr=[]
                s=len(nbs)
                l=99
                for ind2 in range(0,len(nbs)):
                    for k in rg:

                        if (nbs[s-ind2-1]==k[2] and nbs[s-ind2-2]==k[-1] ):
                            if len(k)<=l:
                                l=len(k)
                    for k in rg:
                        if (nbs[s-ind2-1]==k[2] and nbs[s-ind2-2]==k[-1] ):
                            if len(k)==l:
                                if k in cr:
                                    pass
                                else:
                                    cr.append(k)
                rings.append(cr)

    rgs=[]
    pt=[]
    for r in rings:
        for i in r:
            s=([int(j) for j in i[1:]])
            m=min(s)
            rgs.append(s)
            pt.append(s.index(m))

    rgs2=[]
    for ind2,j in enumerate(rgs):
        s=[]
        for ind,k in enumerate(j):
            s.append(j[-ind+pt[ind2]])
        rgs2.append(s)

    rgs3=[]
    for i in rgs2:
        if i in rgs3:
            pass
        else:
            rgs3.append(i)

    np.save((inpf+'_rgs3'),rgs3)

    bld=[]
    bdi=[]
    bd=[]
    for i in rgs3:
        bl=[]
        for ind,j in enumerate(i):
            l=len(i)-1
            dx=data[i[l-ind],2]-data[i[l-ind-1],2]
            dy=data[i[l-ind],3]-data[i[l-ind-1],3]
            bl.append(math.sqrt(dx**2+dy**2))
            ij=sorted([i[l-ind],i[l-ind-1]])
            if ij in bdi:
                bd[bdi.index(ij)][1].append(len(i))
            else:
                bdi.append(ij)
                bd.append([bl[-1],[len(i)]])

        bld.append(bl[::-1])

    np.save((inpf+'_bld'),bld)
    r_size = [len(x) for x in rgs3]
    unique, counts = np.unique(r_size, return_counts=True)
    np.savetxt((inpf+'_rings'),np.transpose(np.vstack((unique,counts))),header='Ring_size Number',comments='')
    return rgs3



def recurse(i,ch,cy,gc,chb):
    global anb
    lt=[]
    if cy>0:
        cy=cy-1
        for j in i:
            if j==gc:
                lt.append(['R',j])
            else:
                ii=1*anb[j]
                if ch in ii:
                    ii.remove(ch)
                for p in chb:
                    if p in ii:
                        ii.remove(p)
                chb2=1*chb
                chb2.append(j)
                r=recurse(ii,j,cy,gc,chb2)
                for k in r:
                    if r==[]:
                        pass
                    else:
                        k.append(j)
                        lt.append(k)
    else:
        pass
    return lt


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
