r"""
    Only for 2D.
    Calculate structure parameters like PDF(g(r)), Cordination Number (CN), Structure Factor S(Q) and .....

    """

import numpy as np

def pdf(inpf,outf,ls,bs,natom,cut,dr,CN_cut):
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
        ls      : Length scaling,  [float]
        bs      : Bredth scaling,  [float]
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
    data2=np.loadtxt(inpf,skiprows=1)

    #scalling of coordinates i.e xs->x
    data=data2
    data[:,2]=data[:,2]*ls
    data[:,3]=data[:,3]*bs

    l=(max(data[:,2])-min(data[:,2]))
    b=(max(data[:,3])-min(data[:,3]))

    #Number of cells in all region(x and y)
    n=math.floor(l/cut+0.000001)
    m=math.floor(b/cut+0.000001)
    cuts=cut**2

    #size of each cell(x and y)
    print('n,m',n,m)
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
    for x in range(0,n):
        leaveit.append(x)
        leaveit.append(n*(m-1)+x)
    for x in range(0,m):
        leaveit.append(n*x)
        leaveit.append(n*x+n-1)


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
    gr=np.zeros((natom,math.floor(cut/dr)+1))

    #nested loop over all atoms to calculate CN and PDF
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
                    r=dx**2+dy**2
                    val=gr[atom_ind]
                    for x in r:
                        if x<cuts:
                            x=math.sqrt(x)
                            rr=math.floor(x/dr)
                            val[rr]=val[rr]+1
                    gr[atom_ind]=val
                gr[atom_ind,0]=data[atom_ind,0]

            for ind2,atom_ind in enumerate(i):
                x0=data[atom_ind,2]
                y0=data[atom_ind,3]
                for ns in ids1:
                    for k in cell[ns]:
                        x1=data[k,2]+nb[ns][2]*l
                        y1=data[k,3]+nb[ns][1]*b
                        dx=x1-x0
                        dy=y1-y0
                        r=dx**2+dy**2
                        val=gr[atom_ind]

                        if r<cuts:
                            r=math.sqrt(r)
                            rr=math.floor(r/dr)
                            val[rr]=val[rr]+1
                        gr[atom_ind]=val
                gr[atom_ind,0]=data[atom_ind,0]

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
    V=l*b
    pdf=[]
    rd=[]
    tr=[]
    for j in range(0,len(grn)):
        rd.append((j+1)*dr+dr/2)
        dv=2*3.14*rd[j]*dr
        pdf.append((grn[j]/dv)/(natom/V))
        #tr.append(pdf[j]*rd[j]*4*3.14*natom/V)

    unique=np.asarray(unique)
    counts=np.asarray(counts)

    #np.savetxt(('grn_'+inpf),np.transpose(np.vstack((rd,grn))),header='r(A) average_atoms',comments='')
    np.savetxt((outf+'_pdf'),np.transpose(np.vstack((rd,pdf))),header='r(A) g(r)',comments='')
    #np.savetxt((outf+'_tr'),np.transpose(np.vstack((rd,tr))),header='r(A) t(r)',comments='')
    np.savetxt((outf+'_CNh'),np.transpose(np.vstack((gr[:,0],CN_all))),header='atom_id CN',comments='')
    np.savetxt((outf+'_CN'),np.transpose(np.vstack((unique,counts))),header='#Average Bond Length(A) : '+str(bond_length)+'#Average CN : '+str(CN_ar)+'    #CN : '+str(CN_dist),comments='')
    np.savetxt((outf+'_gr'),gr)


def SQ(inpf,outf,R,Q,dq,rho,c_int):
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
        c_int  : Drop first few points [float]
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
    q=q[c_int:]
    sq=sq[c_int:]
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
