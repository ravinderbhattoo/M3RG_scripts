r"""
    Import trajectory and perform different operations.
    Dependency:
        numpy
        pandas
        math
        lmps.frames
    """
import pandas as pd
import numpy as np
import math
from lmps import frames as fr
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

p_color = ['#1f77b4',
    '#ff7f0e',
    '#2ca02c',
    '#d62728',
    '#9467bd',
    '#8c564b',
    '#e377c2',
    '#7f7f7f',
    '#bcbd22',
    '#17becf',]

class trajectory:
    frames = None
    names = None
    box = None
    natoms = None
    nframes = None
    fn = None
    compute = {}
    compute_dict = {}
    chuncks = {}
    xchuncks = {}
    figures = {}
    func = {'sum': np.sum,
            'average': np.average,
            'std': np.std,
            'max': np.max,
            'min': np.min,}

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

    def compute_per_chunck(self,cols=[],ids=[],chuncks_id = '',op = 'sum',mvavg=1):
        if chuncks_id == '':
            print('No chuncks_id provided')
        else:
            for indx,col in enumerate(cols):
                w_sum = []
                for ind,f in enumerate(self.frames):
                    sum = []
                    for c in self.chuncks[chuncks_id][ind]:
                        sum.append(self.func[op](f[col][c].values))
                    w_sum.append(movavg(sum,mvavg))
                self.compute[ids[indx]]=np.array((w_sum))
                self.compute_dict[ids[indx]] = chuncks_id

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


    def iplot_figure(self,x=None,y=None,names=None,fix=None,mode='lines',type='scatter',marker={},line={},flip=False,scale=True):
        line0 = line.copy()
        line0['dash']='dot'

        if (y == None) or (x == None) or (names == None):
            raise ValueError('x, y, names are neccessory.')
        else:
            min_x,max_x,abs_x = cal_abs(x)
            min_y,max_y,abs_y = cal_abs(y)

            if scale==False:
                abs_x = [1]*len(abs_x)
                abs_y = [1]*len(abs_y)

            figure = {}

            fix_data = []
            if fix==None:
                pass
            else:
                fix_data = make_data(x=x,y=y,names=names,
                        group=names,abs_y=abs_y,time=fix,mode='lines',marker=marker,line=line0)
            figure['data'] = make_data(x=x,y=y,names=names,group=names,abs_y=abs_y,time=0,
                                mode=mode,marker=marker,line=line) + fix_data
            if flip:
                time = [i for i in range(0,len(x[0]))]
                xlabel = 'time/frame'
                slabel = 'bin/position'

            else:
                time = self.fn
                xlabel = 'bin/position'
                slabel = 'time/frame'

            figure['layout'] = make_layout(time = time,range_x = [min_x,max_x],range_y = [min_y,max_y],xlabel=xlabel,slabel=slabel)

            figure['frames'] = make_frames(time = time,x = x,y = y,names = names,group = names,abs_y = abs_y, mode=mode,marker=marker,line=line)

        return figure

    def iplot_computes(self,flip=False,names=[],fix=None,mode='lines',type='scatter',marker={},line={},scale=True):
        x = []
        y = []
        for n in names:
            if flip:
                y.append(np.transpose(self.compute[n]))
                x.append([self.fn]*len(self.xchuncks[self.compute_dict[n]]))
            else:
                x.append([self.xchuncks[self.compute_dict[n]]]*len(self.fn))
                y.append(self.compute[n])
        return self.iplot_figure(x = x, y = y, names = names,
                        fix=fix,mode=mode,type=type,marker=marker,line=line,flip=flip,scale=scale)


        if name == '':
            return iplot(fig)
        else:
            return iplot(self.figures[name])


    def plot(self,fig=dict(),name=''):
        if name == '':
            return iplot(fig)
        else:
            return iplot(self.figures[name])



def make_data(x=[],y=[],names=[],group=[],abs_y=[],time=0,mode='',marker={},line={}):
    ret = []
    global p_color
    if len(x)>len(p_color):
        p_color *= 1+int(len(x)/len(p_color))

    for ind,i in enumerate(y):
        marker['color'] = p_color[ind]
        line['color'] = p_color[ind]
        ret.append(dict(x=x[ind][time],
                     y=y[ind][time]/abs_y[ind],
                     mode=mode,
                     marker=marker.copy(),
                     line=line.copy(),
                     legendgroup=group[ind]+' inits',
                     name=names[ind]+"/{:.2f}".format(abs_y[ind]) ))
    return ret

def make_layout(time = [],range_x = [],range_y = [],xlabel='',slabel='Value'):
    myaxis = default_axis({'autorange': False,'showgrid':False,})
    return { 'xaxis': myaxis.copy().add({'range': range_x,'title':xlabel}),
             'yaxis': myaxis.copy().add({'range': [-1,1] }),
             'title': 'Overview',
             'updatemenus': [{'type': 'buttons',
                              'buttons': [{'label': 'Play',
                                           'method': 'animate',
                                           'args': [None, {'frame': {'duration': '', 'redraw': True},
                                                     'fromcurrent': True,
                                                     'transition':  {'duration': '', 'easing': 'quadratic-in-out'}
                                                          }],
                                              }]}],
             'sliders':[dict(currentvalue = {"prefix": slabel+": "},pad={'t':50},steps=[{'args':[[str(t)],{'frame':{'duration':'', 'redraw': True},'mode':'immediate','transition':{'duration':''}}],
                                        'label':str(t),'method':'animate',} for t in time])
                        ],
            }


def make_frames(time = [],x=[],y=[],names=[],group=[],abs_y=[],mode='',marker={},line={}):
    return [{'data':make_data(x = x,y = y,names = names,group=names,abs_y = abs_y, time=ind,mode=mode,marker=marker,line=line), 'name': str(t)} for ind,t in enumerate(time)]


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

class default_axis:
    def __init__(self,a_dict):
        self.d_axis = dict(
            showgrid=False,
            showline=True,
            mirror='ticks',
            gridcolor='#eeeeee',
            gridwidth=1,
            zeroline=False,
            linewidth=1,
            ticks='inside',
            ticklen=8,
            tickwidth=1,
            tickformat='',
            exponentformat = 'none',
            showticklabels=True,
            tickprefix = '<span style="font-weight:bold">')
        for i in a_dict:
            self.d_axis[i] = a_dict[i]
    def add(self,a_dict):
        new_dict = self.d_axis.copy()
        for i in a_dict:
            new_dict[i] = a_dict[i]
        return new_dict
    def copy(self):
        return self
