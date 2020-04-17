"""
LAMMPS log file module.

"""

import numpy as np
import pandas as pd
import sys

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

def reader(filename):
    """
    Read LAMMPS log files.

    :param filename: Path of log file

    """
    data = []
    tables = dict(df=[],label=[])
    label = {'a':'Active Fix:'}
    with open(filename,'r') as f:
        rd = 0
        for line in f:
            l = line.strip()
            if rd:
                if l[:4]=='Loop':
                    rd = 0
                    has2end = 0
                    tables['label'].append('\n'.join(label.values()))
                    tables['df'].append(pd.DataFrame(np.array(data),columns=names))
                    data = []
                    label.pop('min',None)
                else:
                    try:
                        data.append([float(i) for i in l.split()])
                    except:
                        pass
            else:
                if l[:4]=='Step':
                    rd = 1
                    has2end = 1
                    names = l.split()
                elif l[:3] == 'fix':
                    label[l.split()[1]] = ' '.join(l.split()[2:])
                elif l[:5] == 'unfix':
                    label.pop(l.split()[1])
                elif l[:8] == 'minimize':
                    label['min'] = 'minimize'
                else:
                    pass
        if has2end:
            rd = 0
            has2end = 0
            tables['label'].append('\n'.join(label.values()))
            tables['df'].append(pd.DataFrame(np.array(data[:-1]),columns=names))
            data = []
            label.pop('min',None)
    return tables


class plotter:
    global reader
    def __init__(self,filename,func = None,figsize=[8,6]):
        self.filename = filename

        if func==None:
            self.tables = reader(filename)
        else:
            self.tables = func(filename)

        self.run = 0
        self.cum = False
        self.prop = None
        self.x = None
        self.axcolor = 'lightgoldenrodyellow'
        self.make_fig()
        self.redraw_all()

    def show(self):
        plt.show()

    def make_fig(self):
        self.fig, (self.ax) = plt.subplots(figsize=[10,7])
        plt.subplots_adjust(left=0.25, right=0.80, bottom=0.25)

    def redraw_all(self):
        self.fig.clear()

        all_names = []
        for i in self.tables['df']:
            all_names += list(i.columns)
        self.names = list(set(all_names))
        self.prop = self.names[0]
        self.x = self.names[0]

        self.draw_plot()
        self.draw_radio()
        self.draw_slider()
        self.draw_reset()
        self.draw_cum()
        self.draw_save()
        self.draw_export()

        self.fig.canvas.draw_idle()


    def draw_plot(self):
        x = self.tables['df'][self.run][self.names[0]]
        y = self.tables['df'][self.run][self.names[0]]
        self.l, = plt.plot(x, y,'b')
        self.ax = plt.gca()
        self.ax.set_xlabel(self.x)
        self.ax.set_ylabel(self.prop)
        self.ax.set_xlim([x.min(),x.max(),])
        self.ax.set_ylim([y.min(),y.max(),])
        self.ax.legend([self.tables['label'][self.run]])

    def update_plot(self):
        if self.cum:
            x = []
            y = []
            for df in self.tables['df']:
                x += list(df[self.x])
                y += list(df[self.prop])
            x = np.array(x)
            y = np.array(y)
        else:
            x = self.tables['df'][self.run][self.x]
            y = self.tables['df'][self.run][self.prop]
        self.ax.set_xlabel(self.x)
        self.ax.set_ylabel(self.prop)
        self.ax.set_xlim([x.min(),x.max(),])
        self.ax.set_ylim([y.min(),y.max(),])
        self.ax.legend([self.tables['label'][self.run]])
        self.l.set_ydata(y)
        self.l.set_xdata(x)
        self.fig.canvas.draw_idle()


    def draw_reset(self):
        resetax = plt.axes([0.7, 0.015, 0.1, 0.04])
        self.button = Button(resetax, 'Reset', color=self.axcolor, hovercolor='0.975')
        self.button.on_clicked(self.reset)

    def reset(self,event):
        self.srun.reset()
        self.update_plot()

    def draw_cum(self):
        resetax = plt.axes([0.55, 0.015, 0.1, 0.04])
        self.button2 = Button(resetax, 'Toggle All', color=self.axcolor, hovercolor='0.975')
        self.button2.on_clicked(self.cumm)

    def cumm(self,event):
        if self.cum:
            self.cum = False
        else:
            self.cum = True
        self.update_plot()

    def draw_save(self):
        resetax = plt.axes([0.4, 0.015, 0.1, 0.04])
        self.button3 = Button(resetax, 'Save Fig', color=self.axcolor, hovercolor='0.975')
        self.button3.on_clicked(self.save_it)

    def save_it(self,event):
        fig, ax = plt.subplots(figsize=[6,6])
        x = self.tables['df'][self.run][self.x]
        y = self.tables['df'][self.run][self.prop]
        ax.plot(x,y)
        ax.set_xlabel(self.x)
        ax.set_ylabel(self.prop)
        plt.savefig('./Figure_eps.eps')

    def draw_export(self):
        resetax = plt.axes([0.25, 0.015, 0.1, 0.04])
        self.button4 = Button(resetax, 'Export Data', color=self.axcolor, hovercolor='0.975')
        self.button4.on_clicked(self.export_it)

    def export_it(self,event):
        x = np.array(self.l.get_xdata()).reshape(-1,1)
        y = np.array(self.l.get_ydata()).reshape(-1,1)
        data = np.hstack([x,y])
        np.savetxt('./export.txt',data)

    def draw_slider(self):
        n = len(self.tables['df'])
        self.maxrun = n-1
        run = plt.axes([0.25, 0.1, 0.55, 0.03], facecolor=self.axcolor)
        self.srun = Slider(run, 'Run block', 0, n,valfmt="%i",valinit=self.run)
        self.srun.on_changed(self.update_slider)

    def update_slider(self,val):
        self.run  = min(int(val),self.maxrun)
        self.update_plot()

    def draw_radio(self):
        rax = plt.axes([0, 0, 0.15, 1], facecolor=self.axcolor,) #aspect='equal')
        self.radio = RadioButtons(rax, self.names, active=np.argmax(self.prop==self.names))
        self.radio.on_clicked(self.radiofunc)

        rax = plt.axes([0.85, 0, 0.15, 1], facecolor=self.axcolor,) #aspect='equal')
        self.radio2 = RadioButtons(rax, self.names, active=np.argmax(self.prop==self.names))
        self.radio2.on_clicked(self.radiofunc2)

    def radiofunc(self,label):
        self.prop = label
        self.update_plot()

    def radiofunc2(self,label):
        self.x = label
        self.update_plot()


if __name__ == '__main__':
    lp = plotter(sys.argv[1])
    lp.show()
