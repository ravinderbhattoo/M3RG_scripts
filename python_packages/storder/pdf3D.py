r"""
    Calculate PDF(g(r)) and Cordination Number (CN)

    How to calculate the pair correlation function g(r)
    This explanation is for three-dimensional data. To calculate g(r), do the following:
    1.Pick a value of dr
    2.Loop over all values of r that you care about:
    3.Consider each particle you have in turn. Count all particles that are a distance between r and r + dr away from the particle you're considering. You can think of this as all particles in a spherical shell surrounding the reference particle. The shell has a thickness dr.
    4.Divide your total count by N, the number of reference particles you considered -- probably the total number of particles in your data.
    5.Divide this number by 4 pi r^2 dr, the volume of the spherical shell (the surface area 4 pi r^2, multiplied by the small thickness dr). This accounts for the fact that as r gets larger, for trivial reasons you find more particles with the given separation.
    6.Divide this by the particle number density. This ensures that g(r)=1 for data with no structure. In other words, if you just had an arbitrarily placed spherical shell of inner radius r and outer radius r+dr, you'd expect to find about rho * V particles inside, where rho is the number density and V is the volume of that shell.
    7.In 2D, follow the algorithm as above but divide by 2 pi r dr instead of step #3 above.

    Parameters
    ----------

    inpf   : Input file name, [String] (format: id, atom type, xs, ys, zs)
    l      : Length of region, Angstrom [float]
    b      : Bredth of region, Angstrom [float]
    natom  : Number of atoms [int]
    cut    : Cut off  distance, Angstrom [float]
    dr     : Step size, Angstrom [float]
    CN_cut : Cutoff for CN, Angstrom [float]
    ----------
    """
import numpy as np
import math
from collections import Counter

def grn_cal(inpf,lx,ly,lz,cut,dr,CN_cut):

    #load data
    data2 = np.loadtxt(inpf,skiprows=1)

    #scalling of coordinates i.e xs->x
    data = data2
    data[:,2] = data[:,2]*lx
    data[:,3] = data[:,3]*ly
    data[:,4] = data[:,4]*lz
    types = Counter(data[:,1]).keys()
    types = [x for x in types]
    freq = Counter(data[:,1]).values()
    permutations = int(len(types)*(len(types)+1)/2)
    lup = [i for i in range(permutations)]
    count=0
    lookup=np.zeros((len(types),len(types)),dtype=int)
    for i in range(len(types)):
        for j in range(len(types)-i):
            lookup[i,i+j] = lup[count]
            lookup[i+j,i] = lup[count]
            count += 1
                            
    natom = len(data)

    #Number of cells in all region(x and y)
    nx = math.floor(lx/cut-0.001)+1
    ny = math.floor(ly/cut-0.001)+1
    nz = math.floor(lz/cut-0.001)+1
    cut_squared = cut**2

    #size of each cell(x and y)
    xbin = lx/nx
    ybin = ly/ny
    zbin = lz/nz

    print(nx,ny,nz)




    #neighbourlist for each cell
    nbx = [-1,0,1]
    nby = [-1,0,1]
    nbz = [-1,0,1]

    print(nx,ny,nz)
    #Empty cells as list
    cell=[[[[] for i in range(nz)] for j in range(ny)] for k in range(nx)]

    circle_x = [i for i in range(nx)] ; circle_x += [0,nx-1]
    circle_y = [i for i in range(ny)] ; circle_y += [0,ny-1]
    circle_z = [i for i in range(nz)] ; circle_z += [0,nz-1]

    flagx = [1] + [0]*(nx-2) + [1]
    flagy = [1] + [0]*(ny-2) + [1]
    flagz = [1] + [0]*(nz-2) + [1]
    
    print('Total Cells : ',len(cell)*len(cell[0])*len(cell[0][0]))

    #insert index of each atom in cell
    for atom_index,atom_coord in enumerate(data):
        ix=math.floor(atom_coord[2]/xbin*0.999)
        iy=math.floor(atom_coord[3]/ybin*0.999)
        iz=math.floor(atom_coord[4]/zbin*0.999)
        cell[ix][iy][iz].append(atom_index)

    #initialise zeros array for every atom (to insert no. of atoms at r distance)
    grn=np.zeros((permutations,natom,math.floor(cut/dr)+1))

    for i,plane in enumerate(cell):
        print('Plane : ',i)
        for j,line in enumerate(plane):
            for k,box in enumerate(line):
                for nbi in nbx:
                    for nbj in nby:
                        for nbk in nbz:
                            nblist=cell[circle_x[i+nbi]][circle_y[j+nbj]][circle_z[k+nbk]]
                            
                            for ind,atom in enumerate(box):
                                ii=types.index(data[atom,1])
                                x0=data[atom,2]
                                y0=data[atom,3]
                                z0=data[atom,4]
                                r1=lx*flagx[i]*nbi + data[nblist,2] - x0
                                r2=ly*flagy[j]*nbj + data[nblist,3] - y0
                                r3=lz*flagz[k]*nbk + data[nblist,4] - z0
                                r=r1**2+r2**2+r3**2
                                for ind1,dist in enumerate(r):
                                    if dist<cut_squared:
                                        jj=types.index(data[nblist[ind1],1])
                                        loc=math.floor(math.sqrt(dist)/dr)
                                        grn[lookup[ii,jj]][atom][loc] += 1
                                        
    return grn

def get_pdf(inpf,lx,ly,lz,cut,dr,CN_cut,grn=None):
    if grn == None:
        grn = grn_cal(inpf,lx,ly,lz,cut,dr,CN_cut)
    else:
        pass
    
    for i in grn:
        natoms = len(i)
        vol = lx*ly*lz
        g = np.sum(i,0)/natoms
        r = [dr+dr*(x+0.5) for x in range(len(g)-1)]
        pdf = 1*g[1:]
        for ind,j in enumerate(r):
            d_vol = 4*math.pi* j**2 *dr
            pdf[ind] = (g[ind]/d_vol)/(natoms/vol)

    return pdf,r     

                                        
def get_CN(inpf,lx,ly,lz,cut,dr,CN_cut,grn=None):
    if grn == None:
        grn = grn_cal(inpf,lx,ly,lz,cut,dr,CN_cut)
    else:
        pass

    
