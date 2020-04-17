import numpy as np
import scipy as sp
from scipy.optimize import curve_fit

def ISF():
    pass

def Ornstein_Zernike_function(r,E,K):
    r = np.array(r)
    return  K*(r**(-1/4)*np.exp(-r/E))

def fit_Ornstein_Zernike_function(r,y,p0=None,dim=2):

    def check_recusive(lt,y):
        if y>lt[-1]:
            lt.pop()
            return check_recusive(lt,y)
        else:
            return lt+[y]

    xs = [0]
    ys = [np.inf]
    y_last = 0
    for x_,y_ in zip(r,y):
        if y_<y_last:
            if y_last>ys[-1]:
                ys = check_recusive(ys,y_last)
                xs = xs[:len(ys)-1]
                xs.append(x_last)
            else:
                ys.append(y_last)
                xs.append(x_last)
        y_last = y_
        x_last = x_

    if dim==2:
        popt, pcov = curve_fit(Ornstein_Zernike_function, xs[1:], ys[1:],p0=p0)
        return popt,(xs[1:],ys[1:])

def stretched_exp_func(t,ta,B,A):
    return A*np.exp(-(t/ta)**B)

def fit_stretched_exp_func(x,y,p0=None):
    popt, pcov = curve_fit(stretched_exp_func, x, y,p0=p0)
    return popt
