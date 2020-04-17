r"""
    Only for 2D
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

    filename   : Input file name, [String] (format: id, atom type, xs, ys, zs)
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
import pandas as pd
import matplotlib.pyplot as plt
import lmps.frames as fr
import os

class  storder:
    def __init__(self,filename,type=None):
        self.filename = filename
        if type==None:
            print('\nFile: ',filename)
            self.data, self.box, self.names = fr.get_frames(filename,[1])
            self.atoms = len(self.data[0])
            print('\nAtoms: ', self.atoms)
        self.pdfs = {}
        self.sqs = {}

    def pdf_cal(self,cut=25,dr=0.01,CN_cut=1.92,verbose=False,draw=False,frame=0):
        filename = self.filename+'_'
        data = self.data[frame].values
        box = self.box[frame]

        #scalling of coordinates i.e xs->x
        data[:,2:5] -= data[:,2:5].min(axis=0)

        l = box[0][1] - box[0][0]
        b = box[1][1] - box[1][0]

        if max(data[:,2])-min(data[:,2])<1.001:
            data[:,2]=data[:,2]*l
            data[:,3]=data[:,3]*b

        #pad atoms +x
        mask = data[:,2]<cut
        data = np.concatenate([data,data[mask,:] + np.array([[100000,0,l,0,0]])])
        #pad atoms -x
        mask = (data[:,2] > l-cut) & ((data[:,2] < l))
        data = np.concatenate([data,data[mask,:] + np.array([[100000,0,-l,0,0]])])

        #pad atoms +y
        mask = data[:,3]<cut
        data = np.concatenate([data,data[mask,:] + np.array([[100000,0,0,b,0]])])
        #pad atoms -y
        mask = (data[:,3] > b-cut) & ((data[:,3] < b))
        data = np.concatenate([data,data[mask,:] + np.array([[100000,0,0,-b,0]])])

        l = l+2*cut
        b = b+2*cut

        natom = len(data)

        #Number of cells in all region(x and y)
        n=math.floor(l/cut+0.000001)
        m=math.floor(b/cut+0.000001)
        cuts=cut**2

        #size of each cell(x and y)
        lbin=l/n
        bbin=b/m

        #neighbourlist for each cell
        nb=[]
        for i in range(0,n*m):
            nb.append([i-n,i-n+1,i-n-1,i-1,i,i+1,i+n,i+n+1,i+n-1])

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
            k=(math.floor(atom_coord[3]/bbin))*n + math.floor(atom_coord[2]/lbin)
            cell[k].append(atom_index)


        #initialise zeros array for every atom (to insert no. of atoms at r distance)
        gr=np.zeros((natom,math.floor(cut/dr)+1))

        #nested loop over all atoms to calculate CN
        for ind,i in enumerate(cell):
             ids=[]
             if (ind in leaveit):
                 if verbose:
                     print('Cell no. : ',ind,' has been left')
             else:
                 if verbose:
                     print('Cell no. : ',ind)
                 for ns in nb[ind]:
                         if ns>-0.1 and ns<n*m:
                             ids.append(cell[ns])

                 for ind2,atom_ind in enumerate(i):

                     x0=data[atom_ind,2]
                     y0=data[atom_ind,3]

                     for k in ids:
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
                     gr[atom_ind,0]=int(data[atom_ind,0])

        mask=[]
        for ind,row in enumerate(gr):
            if row[0]==0:
                mask.append(ind)

        gr=np.delete(gr,mask,0)
        gr=gr[0:len(gr)-1,:]

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
        for j in range(0,len(grn)):
            rd.append((j+1)*dr+dr/2)
            dv=2*3.14*((j+1)*dr+dr/2)*dr
            pdf.append((grn[j]/dv)/(natom/V))

        #calculating PDF/atom
        V=l*b
        pdfperatom = np.zeros((len(gr),len(rd)))
        for ind,i in enumerate(gr):
            for ind2,j in enumerate(i[1:]):
                dv=2*3.14*((ind2+1)*dr+dr/2)*dr
                pdfperatom[ind][ind2]=((j/dv)/(natom/V))


        df = pd.DataFrame(np.transpose(pdfperatom),columns=gr[:,0])
        df['rd']=rd

        #print(df)
        df.to_csv(filename+str(frame)+'.pdfperatom.txt')

        unique=np.asarray(unique)
        counts=np.asarray(counts)



        np.savetxt((filename+str(frame)+'_grn.txt'),np.transpose(np.vstack((rd,grn))),header='r(A) average_atoms',comments='')
        np.savetxt((filename+str(frame)+'_pdf.txt'),np.transpose(np.vstack((rd,pdf))),header='r(A) g(r)',comments='')
        np.savetxt((filename+str(frame)+'CNh_.txt'),np.transpose(np.vstack((gr[:,0],CN_all))),header='#Average CN : '+str(CN_ar)+'    #CN : '+str(CN_dist),comments='')
        np.savetxt((filename+str(frame)+'CN_.txt'),np.transpose(np.vstack((unique,counts))),header='#Average Bond Length(A) : '+str(bond_length)+'#Average CN : '+str(CN_ar)+'    #CN : '+str(CN_dist),comments='')

        self.pdfs['frame'+str(frame)] = np.transpose(np.vstack((rd,pdf)))


        #plot PDF
        if draw:
            fig,(ax1,ax2) = plt.subplots(2,1,figsize=[8,6],dpi=200)
            ax1.plot(rd,grn,label='grn')
            ax1.set_title(filename+str(frame))
            ax1.set_ylabel('grn')
            ax2.plot(rd,pdf,label='pdf')
            ax2.set_ylabel('pdf')
            ax2.set_xlabel('distance')
            plt.savefig(filename+str(frame)+'_pdf.png')


    def sq_cal(self,R=50,Q=20,dq=0.01,rho=1,c_int=0,frame=0,draw=False):
        r"""
            Only for 2D
            Structure Factor S(Q)

            Parameters
            ----------

            filename   : Input file name for pdf, [String]
            outf   : Output file name for SQ, [String]
            rho    : Density of atoms [String]
            Q      : Cut off  distance, Angstrom [float]
            dq     : Step size, Angstrom [float]
            R      : Half of box size, Angstrom [float]
            c_int  : Drop first few points [float]
            ----------
            """
        filename = self.filename+'_'
        os.makedirs(os.path.dirname(filename+str(frame)+'_SQ.txt'), exist_ok=True)
        its_ok = True

        try:
            p = self.pdfs['frame'+str(frame)]
        except:
            try:
                p=np.loadtxt(filename+str(frame)+'_pdf.txt',skiprows=1)
            except:
                print('Please make PDF first.\n Use pdf_cal() function.')
                its_ok = False
        if its_ok:
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
            np.savetxt((filename+str(frame)+'_SQ.txt'),np.transpose(np.vstack((q,sq))),header='Q S(Q)',comments='')
            self.sqs['frame'+str(frame)] = np.transpose(np.vstack((q,sq)))

        if draw:
            fig,(ax) = plt.subplots(1,1,figsize=[8,6],dpi=200)
            ax.plot(q,sq,label='S(q)')
            plt.savefig(filename+str(frame)+'_SQ.png')
