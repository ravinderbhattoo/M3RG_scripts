import numpy as onp
import jax.numpy as np
from jax import vmap, jit
from functools import partial
from jax.ops import index_update, index_add, index_min, index_max


def cal_inds(r,_1,_2,N):
    return np.where((_1<r)*(_2>r),np.arange(N),N)
vmap_cal_inds = vmap(cal_inds,(0,None,None,None))


def cal_intermediate(pos_ca0,shi_ca0,pos_fa,shi_fa,l,hist,_1,_2,N):
    delta = np.abs(pos_fa - pos_ca0)
    r = np.sqrt((np.minimum(delta, np.abs(delta - l))**2).sum(axis=1))
    inds = vmap_cal_inds(r,_1,_2,N)
    val = shi_fa*shi_ca0
    index_add(hist,inds,val)
    index_update(hist,[0,1],0)
    index_update(hist,[N-1],0)
    return hist

vmap_cal_intermediate = vmap(partial(cal_intermediate),(0,0,None,None,None,None,None,None,None))

def cal_gr_n(atoms, pos, type, box, radius=10, bins=100, shi_n=None, key='', dim=3):
    r_step = radius / bins
    if len(atoms) == 2:
        shi = onp.abs(shi_n).reshape(-1)

        ca, fa = atoms

        def masked_data(ca):
            mask_ca = False
            if ca == None:
                ca = ['all']
                return pos, shi
            else:
                for x in ca:
                    mask_ca = mask_ca | (type == x)
                return pos[mask_ca,:], shi[mask_ca]

        pos_ca, shi_ca = masked_data(ca)
        pos_fa, shi_fa = masked_data(fa)

        l = box[:, 1] - box[:, 0]

        hist, _ = onp.histogram(
            [-1], range=[0, radius], bins=bins)
        _1 = _[:-1]
        _2 = _[1:]

        hist = vmap_cal_intermediate(pos_ca,shi_ca,pos_fa,shi_fa,l,hist,_1,_2,len(hist))

        print(hist.shape)

        rs = (_[1:] + _[:-1]) / 2

        if dim == 3:
            div = onp.array(
                [4 * onp.pi * i**2 * r_step for i in rs])
            vol = l.prod()
        elif dim == 2:
            div = onp.array(
                [2 * onp.pi * i * r_step for i in rs])
            vol = l[:-1].prod()
        else:
            raise()

        hist = hist / len(pos_ca) / div / len(pos_fa) * vol
        key += '{}-{}'.format(','.join([str(i) for i in ca]),
                             ','.join([str(i) for i in fa]))
        return hist, rs

def make_pdf(self, atoms=[None, None], radius=30, bins=300, frame_index=None, cal_sq=True, Q=30, q_bins=300, R=50, rho=[1]):
    r"""
    Calculates pair distribution function. (https://en.wikipedia.org/wiki/Radial_distribution_function)

    Parameters
    ----------
    :type atoms: list
    :param atoms: [a,b] where a and b are pair of atoms

    :type radius: float, int
    :param radius: max r value PDF (radius=30)

    :type bins: int
    :param bins: Number of bins (bins=100)

    :type frame_index: int
    :param frame_index: Frame number from trajectory (frame_index=None)

    :type cal_sq: bool
    :param cal_sq: if Ture, it will calculate SQ (default: True)

    :type Q: float, int
    :param Q: max Q for SQ (Q=30)

    :type q_bins: int
    :param q_bins: Number of bins for SQ (q_bins=100)

    :type R: float
    :param R: R for calculating SQ, generally half of box size (R=50)

    :type rho: float
    :param rho: density for calculating SQ (rho=1)

    """

    print('Calculating PDF for: ', atoms)
    self.radius = radius
    self.r_bins = bins
    self.r_step = radius / bins
    self.do_cal_sq = cal_sq
    self.sq_prop = {'R': R, 'dq': Q/q_bins, 'Q': Q, 'rho': rho}

    if len(atoms) == 2:
        if frame_index == None:
            frame_index = range(len(self.frames))

        for f, b, i in zip([self.frames[i] for i in frame_index], [self.box[i] for i in frame_index], frame_index):
            print('Frame index: ', i)
            self.c_frame = f
            self.c_box = b
            self.c_fn = i
            self.cal_pdf(atoms)
    else:
        pass

def cal_pdf(self, atoms):
    if atoms == []:
        pass
    else:
        ca, fa = atoms

        mask_ca = False
        if ca == None:
            ca = ['all']
            pos_ca = self.c_frame[['x', 'y', 'z']].values
        else:
            for x in ca:
                mask_ca = mask_ca | (self.c_frame['type'] == x)
            pos_ca = self.c_frame[['x', 'y', 'z']].values[mask_ca,:]

        mask_fa = False
        if fa == None:
            fa = ['all']
            pos_fa = self.c_frame[['x', 'y', 'z']].values
        else:
            for x in fa:
                mask_fa = mask_fa | (self.c_frame['type'] == x)
            pos_fa = self.c_frame[['x', 'y', 'z']].values[mask_fa,:]

        l = self.c_box[:, 1] - self.c_box[:, 0]

        hist, _ = onp.histogram(
            [-1], range=[0, self.radius], bins=self.r_bins)

        for i in pos_ca:
            delta = onp.abs(pos_fa - i)
            r2 = (onp.minimum(delta, onp.abs(delta - l))**2).sum(axis=1)
            r2 = r2[r2<(self.radius)**2]
            h, _ = onp.histogram(
                onp.sqrt(r2), range=[1e-5, self.radius], bins=self.r_bins)
            hist += h


        rs = (_[1:] + _[:-1]) / 2
        if self.dim == 3:
            div = onp.array(
                [4 * onp.pi * i**2 * self.r_step for i in rs])
            vol = l.prod()
        elif self.dim == 2:
            div = onp.array(
                [2 * onp.pi * i * self.r_step for i in rs])
            vol = l[:-1].prod()
        else:
            raise()

        hist = hist / len(pos_ca) / div / len(pos_fa) * vol
        key = '{}-{}'.format(','.join([str(i) for i in ca]),
                             ','.join([str(i) for i in fa])) + '_' + str(self.c_fn)
        self.PDF[key] = [hist.copy(), rs.copy()]
        self.PDF_rho[key] = len(pos_fa) / vol
        self.TCF[key] = [hist * rs, rs.copy()]
        if self.do_cal_sq:
            self.cal_sq([key], **self.sq_prop)

def cal_sq(self, keys, Q=20, dq=0.1, R=50, rho = None):
    for key in keys:
        print('Calculating SQ for: ', key)
        if key not in self.PDF.keys():
            print('No such PDF exists !!!')
        else:
            p, r = self.PDF[key]
            rho_ = self.PDF_rho[key]
            dr = r[1]-r[0]
            sq_list = onp.zeros((int(Q/dq)))
            q_list = onp.zeros((int(Q/dq)))
            for i in range(1,1+int(Q/dq)):
                q_list[i-1] = dq*i
                # F = sin(pi*r/R) / (pi*r/R)
                # sq = integral [2*pi*r*(pdf-1)*sin(q*r)/q/r ] dr

                F = onp.sin(onp.pi*r/R) / (onp.pi*r/R)
                sq_list[i-1] = 1+rho_*(2*onp.pi*r*(p-1)*onp.sin(dq*i*r)/(dq*i*r)*F*dr).sum()

            self.SQ[key] = [sq_list,q_list]
