"""
Import trajectory and perform different operations.

Dependency:
------------
        numpy
        pandas
        math
        kurma.frames

"""


import pandas as pd
import numpy as np
import math
from kurma import frames as fr
from kurma.structure_order.pdf import PDF as PDF

import matplotlib.pyplot as plt
from tnks import static, dynamic

class trajectory:
    def __init__(self,inpf='',fn=[]):
        if inpf == '':
            pass
        else:
            if fn == []:
                fn = range(1,fr.get_frame_total(inpf)+1)
            else:
                pass

            self.natoms = fr.get_natoms(inpf)
            self.nframes = fr.get_frame_total(inpf)
            self.frames,self.box,self.names = fr.get_frames(inpf,fn)
            self.fn = fn[0:len(self.frames)]
            self.compute = {}
            self.compute_dict = {}
            self.chuncks = {}
            self.xchuncks = {}
            self.figures = {}


    def show_frame(self,index = None,frame = None):
        if index == None:
            index = self.fn.index(frame)

        if index + 1 > self.nframes:
            print("index ({}) excedes total frames ({}) !!!".format(index,self.nframes))
        else:
            print(self.frames[index])


    def make_chuncks(self,name='along_z',dir='z',step=None,start=None,end=None,n=2, logical = 'True' ):
        chuncks = []
        if end<start:
            print('end < start !!')
        else:
            if step == None:
                step = abs(end-start)/n
            else:
                n = int(abs(end-start)/step)+1
            for ind,f in enumerate(self.frames):
                c = []
                for i in range(n):
                    c.append((f[dir] > start+i*step) & (f[dir] < start+(i+1)*step) & eval(logical))
                chuncks.append(c)
            self.xchuncks[name] = [start+(i+0.5)*step for i in range(n)]
            self.chuncks[name] = chuncks

    def compute_per_chunck(self,cols=[[]],ids=[],chuncks_id = '',op = 'sum',f_args=[],f_kwargs={},mvavg=1,trickle_self=False):
        if chuncks_id == '':
            print('No chuncks_id provided')
        else:
            for indx,col in enumerate(cols):
                w_sum = []
                if trickle_self:
                    op(col,ids[indx],self=self,*f_args,**f_kwargs)
                else:
                    for ind,f in enumerate(self.frames):
                        sum = []
                        for c in self.chuncks[chuncks_id][ind]:
                            sum.append(op(f[col][c],*f_args,**f_kwargs))
                        w_sum.append(movavg(sum,mvavg))
                    self.compute[ids[indx]]=np.array((w_sum))
                    self.compute_dict[ids[indx]] = chuncks_id

    def compute_full_frame(self,cols=[[]],ids=[],op='average',f_args=[],f_kwargs={},trickle_self=False):
        for indx,col in enumerate(cols):
            w_sum = []
            if trickle_self:
                op(col,ids[indx],self=self,*f_args,**f_kwargs)
            else:
                for ind,f in enumerate(self.frames):
                    f_kwargs.update({'box':self.box[ind]})
                    w_sum.append(op(f[col],*f_args,**f_kwargs))
                self.compute[ids[indx]]=np.array((w_sum))
                self.compute_dict[ids[indx]] = 'full_frame'

    def count_per_chunck(self,ids=[],chuncks_id = ''):
        if chuncks_id == '':
            print('No chuncks_id provided')
        else:
            w_sum = []
            for ind,f in enumerate(self.frames):
                sum = []
                for c in self.chuncks[chuncks_id][ind]:
                    sum.append(f['z'][c].count())
                w_sum.append(sum)
            self.compute[ids[0]]=np.array((w_sum))
            self.compute_dict[ids[0]] = chuncks_id


    def plot_figures(self,x=None,y=None,names=None,fix=None,flip=False,scale=True):
        if (y == None) or (x == None) or (names == None):
            raise ValueError('x, y, names are neccessory.')
        else:
            min_x,max_x,abs_x = cal_abs(x)
            min_y,max_y,abs_y = cal_abs(y)

            if scale==False:
                abs_x = [1]*len(abs_x)
                abs_y = [1]*len(abs_y)
        return  make_frames(time = time,x = x,y = y,names = names,group = names,abs_y = abs_y)


    def plot_computes(self,flip=False,names=[],fix=None, mode='lines',type='scatter',marker={},line={},scale=True):
        x = []
        y = []
        for n in names:
            if flip:
                y.append(np.transpose(self.compute[n]))
                x.append([self.fn]*len(self.xchuncks[self.compute_dict[n]]))
            else:
                x.append([self.xchuncks[self.compute_dict[n]]]*len(self.fn))
                y.append(self.compute[n])
        return self.plot_figures(x = x, y = y, names = names,
                        fix=fix,flip=flip,scale=scale)

    def plot(self,fig=dict(),name=''):
        if name == '':
            return plt.plot(**fig)
        else:
            return plt.plot(**self.figures[name])

def make_data(x=[],y=[],names=[],group=[],abs_y=[],time=0,):
    ret = []
    for ind,i in enumerate(y):
        ret.append(dict(x=x[ind][time],
                     y=y[ind][time]/abs_y[ind],
                     label=names[ind]+"/{:.2f}".format(abs_y[ind])))
    return ret

def make_layout(time = [],range_x = [],range_y = [-1,1], xlabel='', slabel='Value'):
    def fun_layout(ax, range_x=range_x, range_y=range_y, xlabel=xlabel):
        ax.set_xlim(range_x)
        ax.set_ylim(range_y)
        ax.set_xlabel(xlabel)
        ax.set_title('Overview')
    return fun_layout

def make_frames(time = [],x=[],y=[],names=[],group=[],abs_y=[]):
    return [{'data':make_data(x = x,y = y,names = names,group=names,abs_y = abs_y, time=ind), 'name': str(t)} for ind,t in enumerate(time)]

def cal_abs(p):
    abs = []
    min_i = []
    max_i = []
    for i in p:
        min_i += [np.min(i)]
        max_i += [np.max(i)]
        abs.append(np.max(np.abs([min_i,max_i])))
    return np.min(min_i),np.max(max_i),abs


def update_dict(a,b,i):
    a.update(b)
    return a

def movavg(x,n):
    x1 = x.copy()
    if int(n/2)-n/2 == 0:
        pass
    else:
        for i in range(len(x)-n+1):
            x1[i+int((n)/2)] = sum(x[i:i+n])/n
    return x1


def pdfs(df,id,self=None,ind=None,dim=3,*args,**kwargs):
    pdf = PDF(info = [self.frames, self.box, self.names], dim = dim)
    pdf.make_pdf(frame_index=range(len(self.fn)),*args,**kwargs)
    self.compute[id+'_'+list(pdf.PDF.keys())[0].split('_')[0]] = list(pdf.PDF.values())
    if not((pdf.TCF==None) or (pdf.TCF=={})):
        self.compute[id+'_tcf_'+list(pdf.PDF.keys())[0].split('_')[0]] = list(pdf.TCF.values())
    if not((pdf.SQ==None) or (pdf.SQ=={})):
        self.compute[id+'_sq_'+list(pdf.PDF.keys())[0].split('_')[0]] = list(pdf.SQ.values())

def gr_n(df,id,self=None,ind=None,dim=3,*args,**kwargs):
    pdf = PDF(info = [self.frames, self.box, self.names], dim = dim)
    pdf.make_gr_n(frame_index=range(len(self.fn)),*args,**kwargs)
    self.compute[id+'_'+list(pdf.gr_n.keys())[0].split('_')[0]] = list(pdf.gr_n.values())

def shi_n(df=None,pos=None,box=None,n=3,rm=2,):
    if pos==None:
        pos = df.values
    shi_n = np.zeros((pos.shape[0]),dtype=type(1j))
    l = box[:, 1] - box[:, 0]
    l = l[:2]
    for ind,i in enumerate(pos):
        delta = (pos-i.reshape(1,-1))[:,:2]
        dist = (np.minimum(delta, np.abs(delta - l))**2).sum(axis=1)
        mask = (0.001<dist)*(dist<=rm**2)
        delta = delta[mask]
        delta = delta/l
        delta = np.where(delta>1/2,delta-1,0) + np.where(delta<-1/2,1+delta,0) + np.where((delta<-1/2)+(delta>1/2),0,delta)
        thetas = np.arctan2(*delta.T[::-1])
        try:
            thetas -= thetas[0]
        except:
            thetas = [np.pi/2*n]
        shi_n[ind] = __shi_n(thetas,n=n)
    return shi_n

def __shi_n(thetas,n=3):
    return np.mean(np.exp(1j*n*np.array(thetas)))
