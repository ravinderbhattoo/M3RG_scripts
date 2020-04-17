r"""
    Only for 2D.
    Calculate structure parameters like PDF(g(r)), Cordination Number (CN), Structure Factor S(Q) and .....

    """
import math
import os
import numpy as np
import pandas as pd
from dump import frames as fr
import matplotlib.pyplot as plt
from shadow.plot import *
from scipy.optimize import leastsq


def pdf(inpf,outf,cut,dr,CN_cut,fn=[1]):
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
        cut    : Cut off  distance, Angstrom [float]
        dr     : Step size, Angstrom [float]
        CN_cut : Cutoff for CN, Angstrom [float]
        fn     : Frame Number from trajectory file [int]
        ----------
        """
    if os.path.basename(outf)=='':
        outf = './'+outf
    full_path = os.path.abspath('./'+outf)
    os.makedirs(os.path.dirname(full_path), exist_ok=True)
    dir_,file = os.path.dirname(full_path),os.path.basename(full_path)

    print(full_path)

    #load data
    try:
        data_frames,boxs,names =  fr.get_frames(inpf,fn)
    except:
        print('Errors')
        Error

    RD = []
    pdf_ = []
    tr_ = []

    for fns in range(len(fn)):
        print('Calculating PDF for frame: ',fns)
        data = data_frames[fns]
        box = boxs[fns]

        natoms = len(data)
        print('Number of atoms: ',natoms)

        lx_max = box[0][1] - box[0][0]
        ly_max = box[1][1] - box[1][0]
        lz_max = box[2][1] - box[2][0]

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
        r = []
        for i in range(int(lx_max/cut+1)):
            for j in range(int(ly_max/cut+1)):
                mask_inner = (data['x']<=(i+1)*cut) & (data['x']>=i*cut) & (data['y']<=(j+1)*cut) & (data['y']>=j*cut) & (data['y']<=ly_max) & (data['x']<=lx_max) & (data['x']>=0) & (data['y']>=0)
                mask_outer = (data['x']<(i+2)*cut) & (data['x']>(i-1)*cut) & (data['y']<(j+2)*cut) & (data['y']>(j-1)*cut)

                iter_obj = data[mask_inner]
                comp_obj = data[mask_outer]

                for index, row in iter_obj.iterrows():
                    r.append(np.sqrt((comp_obj['x']-row['x'])**2+(comp_obj['y']-row['y'])**2))

        CN = []
        size = int(cut/dr)+1
        pdf = np.zeros((size))
        tr = np.zeros((size))
        for atom in r:
            CN.append(np.sum((atom<CN_cut) & (atom>0.01)))
            for i in atom:
                try:
                    pdf[int(i/dr)] += 1
                except:
                    pass
        pdf /= len(r)
        rho = natoms/lx_max/ly_max
        for ind,val in enumerate(pdf):
            pdf[ind] = val/(2*np.pi*((ind+0.5)*dr)*dr*rho)
            tr[ind] = pdf[ind]*(ind+0.5)*dr

        pdf = pdf[1:]
        tr = tr[1:]
        rd  = [dr*(i+1+0.5) for i in range(0,len(pdf))]
        values, counts = np.unique(CN, return_counts=True)
        print('CN: ',values)
        print('Counts: ',counts)

        pdf_.append(pdf.copy())
        tr_.append(tr)
        RD = np.array(rd)

    if (len(fn)>1):
        np.savetxt((dir_+'/pdf_'+file),np.hstack((RD.reshape(-1,1),np.array(pdf_).T)), header='r(\AA) g(r) '+' '.join(['g(r)-{}'.format(i+1) for i in range(len(fn)-1)]),comments='')

        np.savetxt((dir_+'/tr_'+file),np.hstack((RD.reshape(-1,1),np.array(pdf_).T)),header='r(\AA) t(r) '+' '.join(['t(r)-{}'.format(i+1) for i in range(len(fn)-1)]),comments='')
    else:
        np.savetxt((dir_+'/pdf_'+file), np.transpose(np.vstack((rd,pdf))),header='r(\AA) g(r)',comments='')
        np.savetxt((dir_+'/tr_'+file), np.transpose(np.vstack((rd,tr))),header='r(\AA) t(r)',comments='')
        np.savetxt((dir_+'/CN_'+file), np.transpose(np.vstack((values,counts))),header='Values Counts',comments='')


def SQ(inpf,outf,R,Q,dq,rho=None,c_int=0):
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
    if os.path.basename(outf)=='':
        outf = './'+outf
    full_path = os.path.abspath('./'+outf)
    os.makedirs(os.path.dirname(full_path), exist_ok=True)
    dir_,file = os.path.dirname(full_path),os.path.basename(full_path)

    print(full_path)

    p=np.loadtxt(inpf,skiprows=1)
    pdf=p[:,1:].T
    if rho==None:
        rho = [1]*len(pdf)
        print('rho is not provided. rho=1')
    if len(rho)!=len(pdf):
        rho = rho*len(pdf)
    print('Size of pdf: ',len(pdf))
    rd=p[:,0]
    dr=rd[2]-rd[1]
    sq_=[]
    ind = 0
    for pdf_ in pdf:
        rho_ = rho[ind]
        ind += 1
        sq = []
        q = []
        for j in range(0,int(Q/dq)):
            sigma=0
            q.append(dq*(j+1))
            for i in range(0,len(pdf_)-1):
                pdfv=(pdf_[i]+pdf_[i+1])/2
                rdv=(rd[i]+rd[i+1])/2
                F=math.sin(math.pi*rdv/R) / (math.pi*rdv/R)
                sig=2*math.pi*rdv * (pdfv-1) * (math.sin(q[j]*rdv)/(q[j]*rdv)) * F * dr
                sigma=sigma+sig

            sq.append(1+rho_*sigma)
        q=np.array(q[c_int:])
        sq_.append(sq[c_int:])

    np.savetxt((dir_+'/SQ_'+file), np.transpose(np.vstack((q,np.array(sq_)))),header='q(1/\AA) S(q) '+' '.join(['S(q)-{}'.format(i+1) for i in range(len(pdf)-1)]),comments='')






def plot(list1,which='avg',ax=None,base=0,shift=0,shift2=0,legends=None,legends2=None,xlim_=None,ylim_=None,**kwargs):
    r"""
        Plot calculated structure parameters PDF(g(r)) and Structure Factor S(Q).

        Parameters
        ----------
        inpf   : Input file name for pdf/SQ, [String]
        ----------
        """
    if legends==None:
        legends = ['']*len(list1)
    do_show = True
    if ax==None:
        fig, [ax] = panel(1,1)
    else:
        do_show = False

    if 1:
        lines = []
        for ind,inpf in enumerate(list1):
            p=np.loadtxt(inpf ,skiprows=1)
            head=np.genfromtxt(inpf,dtype='str',max_rows=1)
            if which=='avg':
                l1 = ax.plot(p[:,0],[base+shift*ind+x for x in p[:,1:].mean(axis=1)],label=legends[ind],**kwargs)
                lines.append(l1[0])
            elif which=='all':
                for i in range(p.shape[1]-1):
                    if legends2==None:
                        label = legends[ind]
                    else:
                        label = legends[ind]+legends2[i]

                    ax.plot(p[:,0],[base+shift*ind+x for x in p[:,i+1]],label=label,**kwargs)
            elif type(which)==type([]):
                for ind2,j in enumerate(which):
                    if legends2==None:
                        label = legends[ind]
                    else:
                        label = legends[ind]+str(j)
                    l1 = ax.plot(p[:,0],[base+shift*ind+ind2*shift2+x for x in p[:,1]],label=label,**kwargs)
                    lines.append(l1[0])

            else:
                l1 = ax.plot(p[:,0],[base+shift*ind+x for x in p[:,1]],label=legends[ind],**kwargs)
                lines.append(l1[0])

        xlabel(r'{}'.format(head[0]),ax=ax)
        ylabel(r'{}'.format(head[1]),ax=ax)
        xlim(xlim_,ax=ax)
        ylim(ylim_,ax=ax)
        legend_on(ax=ax)
        if do_show:
            plt.savefig(inpf+'.png', dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format=None,
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None)
            plt.show()

    else:
        print('failed!!!')

    return lines


def lorentz(p,x):
    if len(p)==3:
        return p[2] * (p[0]/2)/((p[0]/2)**2+(x-p[1])**2)
    else:
        return p[3] + p[2] * (p[0]/2)/((p[0]/2)**2+(x-p[1])**2)

def errorfunc(p,x,z):
        return lorentz(p,x)-z

def fit_lorentz(x,z,p0=[1,1,1]):
    solp, ier = leastsq(errorfunc,
                    p0,
                    args=(x,z),
                    Dfun=None,
                    full_output=False,
                    ftol=1e-9,
                    xtol=1e-9,
                    maxfev=100000,
                    epsfcn=1e-10,
                    factor=0.1)
    return solp, ier












# end
