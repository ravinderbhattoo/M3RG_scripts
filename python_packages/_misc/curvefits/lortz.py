r"""
    Lorentzian fit
    ----------
    """
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

def lortz(x,c,f,base=0):
    ret=(1/math.pi)*0.5*f/((x-c)**2+(0.5*f)**2)+base
    return ret


def fit(xdata,ydata,p0,c='none',f='none',b=0):
    if (c=='none' and f=='none'):
        popt, pcov = curve_fit(lambda x,c,f: lortz(x, c, f, b),xdata,ydata,p0)

    elif (c=='none'):
        popt, pcov = curve_fit(lambda x,c: lortz(x, c, f, b),xdata,ydata,p0)

    else:
        popt, pcov = curve_fit(lambda x,f: lortz(x, c, f, b),xdata,ydata,p0)

    return [popt,pcov]


def get(xdata,ydata,p0,x1,x2,step,center='none',fwhm='none',base=0):
    popt,pcov = fit(xdata,ydata,p0,c=center,f=fwhm,b=base)
    xaxis = np.linspace(x1,x2,step)
    if center!='none':
        c=center
        f=popt
    elif fwhm!='none':
        f=fwhm
        c=popt
    else:
        c=popt[0]
        f=popt[1]


    yaxis=lortz(xaxis,c,f,base)
    err = pcov

    plt.scatter(xdata,ydata,color='r')
    plt.plot(xaxis,yaxis,color='b')


    return [xaxis,yaxis,[c,f],err]





%matplotlib qt5
import pandas as pd
def lfit(i,a,b,center=0,base=0,hz=0):
    filenames = ['glass.dat_SQ','crystalline.dat_SQ','prinstine.dat_SQ']
    for f in filenames[i:i+1]:
        df = pd.read_csv('/Users/ravinder/Main/Graphene/graphene1/'+f,sep=' ')
        x = df['Q']
        y = df['S(Q)']

        plt.plot(x[1000:5000],y[1000:5000])

        mask = (x>a) & (x<b)

        xdata = x[mask]
        ydata = y[mask]

        x_1 = 2
        x_2 = 4

        xy = get(xdata,ydata,[1],x_1,x_2,500,center=center,base=base)

        x1 = xy[0]
        y1 = xy[1]

        yp = max(y1)
        yp2 = (yp+hz)/2

        xvals = []
        yvals = []

        sign = 1
        for ind,i in enumerate(y1):
            if (i-yp2)*sign>0:
                xvals.append(x1[ind])
                yvals.append(y1[ind])
                sign = -1*sign
        width = xvals[1]-xvals[0]
        plt.plot(xvals,yvals)
        plt.plot([x_1,x_2],[hz,hz])

        print('width: ',width)
        print('Cov: ',xy[3])


lfit(0,2.9,3.26,center=3.08,base=0.09,hz=0.35) #0.533/0.4969

lfit(1,3.03,3.185,center=3.11,base=-0.1,hz=0.25) #0.340

lfit(2,3.03,3.16,center=3.097,base=-0.1,hz=0.25) #0.340
