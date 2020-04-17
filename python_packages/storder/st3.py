r"""
    Only for 2D.
    Calculate structure parameters like PDF(g(r)), Cordination Number (CN), Structure Factor S(Q) and .....

    """

import numpy as np

def neigh(i,j,k,dim=3):
    nb=[]
    list1 = [-1,0,1]
    for x in list1:
        for y in list1:
            if dim==3:
                for z in list1:
                    nb.append('x'+str(i+x)+'y'+str(j+y)+'z'+str(k+z))
            else:
                nb.append('x'+str(i+x)+'y'+str(j+y)+'z1')
    return nb

def pdf(inpf,outf,xs,ys,zs,natom,cut,dr,CN_cut,dim=3):
    r"""
        Only for 2D
        Calculate PDF(g(r)), Cordination Number (CN)

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

        inpf   : Input file name, [String] (format: id, atom_type, xs, ys, zs)
        outf   : Output file name, [String]
        xs     : x-scaling,  [float]
        ys     : y-scaling,  [float]
        zs     : z-scaling,  [float]
        natom  : Number of atoms [int]
        cut    : Cut off  distance, Angstrom [float]
        dr     : Step size, Angstrom [float]
        CN_cut : Cutoff for CN, Angstrom [float]
        ----------
        """

    import math
    import os
    os.makedirs(os.path.dirname(outf), exist_ok=True)


    #load data
    data=np.loadtxt(inpf,skiprows=1)

    #scalling of coordinates i.e xs->x
    data[:,2]=data[:,2]*xs
    data[:,3]=data[:,3]*ys
    data[:,4]=data[:,4]*zs

    data[:,2] -= min(data[:,2])
    data[:,3] -= min(data[:,3])
    data[:,4] -= min(data[:,4])

    Lx=max(data[:,2])
    Ly=max(data[:,3])
    Lz=max(data[:,4])

    #Number of cells in all region(x and y)
    if Lx>cut:
        nx=math.floor(Lx/cut+0.000001)
        xbin=Lx/nx

    else:
        nx=1
        xbin=cut
        print('Lx < cut')
    if Ly>cut:
        ny=math.floor(Ly/cut+0.000001)
        ybin=Ly/ny
    else:
        ny=1
        ybin=cut
        print('Ly < cut')
    if Lz>cut:
        nz=math.floor(Lz/cut+0.000001)
        zbin=Lz/nz
    else:
        nz=1
        zbin=cut
        print('Lz < cut')

    x_c = 1*data[:,2]
    y_c = 1*data[:,3]
    z_c = 1*data[:,4]
    ai = 1*data[:,0]
    m1 = x_c<xbin
    m2 = x_c>Lx-xbin
    x_c = np.concatenate((x_c,x_c[m1]+Lx,x_c[m2]-Lx))
    y_c = np.concatenate((y_c,y_c[m1],y_c[m2]))
    z_c = np.concatenate((z_c,z_c[m1],z_c[m2]))
    ai = np.concatenate((ai,ai[m1],ai[m2]))
    m1 = y_c<ybin
    m2 = y_c>Ly-ybin
    x_c = np.concatenate((x_c,x_c[m1],x_c[m2]))
    y_c = np.concatenate((y_c,y_c[m1]+Ly,y_c[m2]-Ly))
    z_c = np.concatenate((z_c,z_c[m1],z_c[m2]))
    ai = np.concatenate((ai,ai[m1],ai[m2]))

    cuts=cut**2

    #size of each cell(x and y)
    print('nx,ny,nz',nx,ny,nz)

    cells = {}
    #insert index of each atom in every cell
    for atom_index,x in enumerate(x_c):
        kx= math.floor(x/xbin+1)
        ky= math.floor(y_c[atom_index]/ybin+1)
        kz= math.floor(z_c[atom_index]/zbin+1)
        id = 'x'+str(kx)+'y'+str(ky)+'z'+str(kz)
        try:
            cells[id].append(atom_index)
        except:
            cells[id] = [atom_index]

    #initialise zeros array for every atom (to insert no. of atoms at r distance)
    gr=np.zeros((len(x_c),math.floor(cut/dr)+1))

    #nested loop over all atoms to calculate CN and PDF
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                name = 'x'+str(i+1)+'y'+str(j+1)+'z'+str(k+1)
                neighbour = neigh(i+1,j+1,k+1,dim)
                for atomid in cells[name]:
                    x0 = x_c[atomid]
                    y0 = y_c[atomid]
                    z0 = z_c[atomid]
                    for nb in neighbour:
                        x1 = x_c[cells[nb]]
                        y1 = y_c[cells[nb]]
                        z1 = z_c[cells[nb]]
                        r2 = (x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2
                        val=gr[atomid]
                        for x in r2:
                            if x<cuts:
                                x=math.sqrt(x)
                                rr=math.floor(x/dr)
                                val[rr]=val[rr]+1
                        gr[atomid]=val

    mask=[]
    gr=gr[:,:-2]
    for ind,row in enumerate(gr):
        if row[0]==0:
            mask.append(ind)
    gr=np.delete(gr,mask,0)

    if CN_cut<cut:
        ind=math.floor(CN_cut/dr)
        CN_all=np.sum(gr[:,1:ind],axis=1)
        unique, counts = np.unique(CN_all, return_counts=True)
        CN_dist=dict(zip(unique, counts))
        CN_ar=np.average(CN_all)

        bond_freq=np.sum(gr[:,1:ind],axis=0)/len(gr)
        bond_length=0
        for ind3,i in enumerate(bond_freq):
            bond_length=bond_length+(ind3+1.5)*i
        bond_length=bond_length*dr/CN_ar
        print('Avg. Bond Length : ',bond_length)
        print('CN : ',CN_dist)
        print('Average CN : ',CN_ar)
    else:
        print('General cutoff is less than cutoff for cordination no.!!!')




    #Average gr for all atom
    s=gr[0,1:]*0
    for x in gr:
        s=x[1:]+s
    grn=s/len(gr)


    #calculating PDF
    V2=Lx*Ly
    V3=Lx*Ly*Lz
    pdf=[]
    rd=[]
    tr=[]
    for j in range(0,len(grn)):
        rd.append((j+1)*dr+dr/2)
        if dim==3:
            pass
        else:
            dv=2*3.14*rd[j]*dr
            rho = natom/V2
        pdf.append((grn[j]/dv)/rho)
        tr.append(pdf[j]*rd[j]*4*3.14*rho)

    unique=np.asarray(unique)
    counts=np.asarray(counts)

    #np.savetxt(('grn_'+inpf),np.transpose(np.vstack((rd,grn))),header='r(A) average_atoms',comments='')
    np.savetxt((outf+'_pdf'),np.transpose(np.vstack((rd,pdf))),header='r(A) g(r)',comments='')
    #np.savetxt((outf+'_tr'),np.transpose(np.vstack((rd,tr))),header='r(A) t(r)',comments='')
    np.savetxt((outf+'_CNh'),np.transpose(np.vstack((gr[:,0],CN_all))),header='atom_id CN',comments='')
    np.savetxt((outf+'_CN'),np.transpose(np.vstack((unique,counts))),header='#Average Bond Length(A) : '+str(bond_length)+'#Average CN : '+str(CN_ar)+'    #CN : '+str(CN_dist),comments='')
    np.savetxt((outf+'_gr'),gr)


def SQ(inpf,outf,R,Q,dq,rho,drop=0):
    r"""
        Only for 2D
        Structure Factor S(Q)

        Parameters
        ----------

        inpf   : Input file name for pdf, [String]
        outf   : Output file name for SQ, [String]
        rho    : Density of atoms [String]
        Q      : Cut off  distance, Angstrom [float]
        dq     : Step size, Angstrom [float]
        R      : Half of box size, Angstrom [float]
        drop  : Drop first few points [float]
        ----------
        """
    import math
    import os
    os.makedirs(os.path.dirname(outf), exist_ok=True)
    #rho = eval(rho)
    p=np.loadtxt(inpf,skiprows=1)
    pdf=p[:,1]
    rd=p[:,0]
    dr=rd[2]-rd[1]
    sq=[]
    q=[]
    for j in range(0,int(Q/dq)):
        sigma=0
        q.append(dq*(j+1))
        for i in range(0,len(pdf)-1):
            pdfv=(pdf[i]+pdf[i+1])/2
            rdv=(rd[i]+rd[i+1])/2
            F=math.sin(math.pi*rdv/R) / (math.pi*rdv/R)
            sig=2*math.pi*rdv * (pdfv-1) * (math.sin(q[j]*rdv)/(q[j]*rdv)) * F * dr
            sigma=sigma+sig

        sq.append(1+rho*sigma)

    mask = 0
    for ind,i in enumerate(sq):
        if i>drop:
            mask = ind
            break

    q=q[mask:]
    sq=sq[mask:]
    np.savetxt((outf+'_SQ'),np.transpose(np.vstack((q,sq))),header='Q S(Q)',comments='')






def plot(list1,shift=0,ymaximum='auto'):
    r"""
        Plot calculated structure parameters PDF(g(r)) and Structure Factor S(Q).

        Parameters
        ----------
        inpf   : Input file name for pdf/SQ, [String]
        ----------
        """
    import matplotlib.pyplot as plt
    try:
        for ind,inpf in enumerate(list1):
            p=np.loadtxt(inpf ,skiprows=1)
            head=np.genfromtxt(inpf,dtype='str',max_rows=1)
            plt.plot(p[4:,0],[shift*ind+x for x in p[4:,1]])
            if ymaximum=='auto':
                pass
            else:
                plt.ylim(-1,ymaximum)
            plt.xlabel(head[0], fontsize=18)
            plt.ylabel(head[1], fontsize=18)
            plt.savefig(inpf+'.png', dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None)

        plt.show()

    except:
        print('failed!!!')
