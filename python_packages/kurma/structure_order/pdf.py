"""
This module is for creating functions related to structure parameters.

"""

from kurma import frames as fr
import numpy as np
import itertools

import matplotlib.pyplot as plt

class PDF:
    r"""

    Example:
    --------

    .. code-block:: python

        from kurma.structure_order import pdf
        import matplotlib.pyplot as plt
        from shadow.plot import xlabel, ylabel, xlim, ylim, panel, xticks
        import numpy as np

        P = pdf.PDF('./frame0_s.trj', dim=2, fn=[1])
        P.make_pdf()
        # P.make_pdf(atoms=[[2],[1]], frame_index=None, bins=200, radius=30)
        leg_prop = {'loc': 1, 'bbox_to_anchor': [0.9, 0.9]}
        P.PDF.keys()

        fig, [ax] = panel(1,1)
        P.plot('all-all', 'PDF', label=r'\textbf{PDF1}', leg_prop=leg_prop, ax=ax, mean=True,shift=5)
        #P.plot('2-1_0', 'PDF', label=r'\textbf{PDF2}', leg_prop=leg_prop, ax=ax)
        xlabel(r'Distance, r (\AA)')
        ylabel(r'Pair distribution function (PDF)')
        plt.legend()
        plt.show()

        fig, [ax] = panel(1,1)
        P.plot('all-all', 'TCF', label=r'\textbf{TCF1}', leg_prop=leg_prop, ax=ax, mean=True,shift=80)
        #P.plot('2-1_0', 'TCF', label=r'\textbf{TCF2}', leg_prop=leg_prop,ax=ax)
        xlabel(r'Distance, r (\AA)')
        ylabel(r'Total correlation function (TCF)')
        plt.legend()
        plt.show()

        fig, [ax] = panel(1,1)
        P.plot('all-all', 'SQ', label=r'\textbf{SQ1}', leg_prop=leg_prop, mean=True, ax=ax, shift=10)
        #P.plot('2-1_0', 'SQ', label=r'\textbf{SQ2}', leg_prop=leg_prop, ax=ax)
        xlabel(r'Wave-vector, q (k)')
        ylabel(r'Structure function (SQ)')
        plt.legend()
        plt.show()

        cutoff={'2-1':1.6}
        c = P.CN(cutoff,frame_index=[0],default=1.6)

        shift = 0
        fig, [ax] = panel(1,1)
        for key in P.bond_angles[0].keys():
            if P.bond_angles[0][key][0] != []:
                a = str(key).split('0')
                P.plot_bond_angles([[int(i) for i in a]],ax=ax,shift=shift,label=r'\textbf{{{}}}'.format(key))
                shift += 10

        ylabel('Frequency')
        xlabel('Bond angle distribution (\AA)')
        xticks(range(0,190,30))
        ylim([-10,50])
        plt.legend()
        plt.show()

    """

    def __init__(self, filename=None, dim=3, fn=[1], info=None):
        r"""
        Parameters:
        -----------

        param: filename: Path of trajectory file

        param: dim: Dimensions of system (2D or 3D) (dim=3)

        param: fn: list of frame numbers (see dump.get_frames for more info.) (fn=[1])

        """

        print('Calculates only for orthogonal box !!!')
        if filename!=None:
            self.frames, self.box, self.names = fr.get_frames(filename, fn=fn)
        else:
            self.frames, self.box, self.names = info

        for i in range(len(self.frames)):
            if 'type' not in self.names:
                self.frames[i]['type'] = np.ones(len(self.frames[i]))
            self.frames[i][['id', 'type']] = self.frames[i][[
                'id', 'type']].astype(int)
            self.frames[i].index = self.frames[i]['id']
        self.dim = dim
        self.PDF = {}
        self.PDF_rho = {}
        self.gr_n = {}
        self.TCF = {}
        self.SQ = {}
        self.nb_list = {}
        self.bond_angles = {}

    def make_gr_n(self, atoms=[None, None], radius=30, bins=300, frame_index=None, shi_n=None, key=''):
        print('Calculating g{}(r) for: '.format(key), atoms)
        r_step = radius / bins

        if len(atoms) == 2:
            if frame_index == None:
                frame_index = range(len(self.frames))

            for f, b, i in zip([self.frames[i] for i in frame_index], [self.box[i] for i in frame_index], frame_index):
                print('Frame index: ', i)
                self.c_frame = f
                self.c_box = b
                self.c_fn = i
                shi = np.abs(shi_n[i]).reshape(-1)
                self.cal_gr_n(atoms,shi,radius, bins, r_step, key)
        else:
            pass

    def cal_gr_n(self, atoms, shi, radius, bins, r_step, key):
        if atoms == []:
            pass
        else:
            ca, fa = atoms

            mask_ca = False
            if ca == None:
                ca = ['all']
                pos_ca = self.c_frame[['x', 'y', 'z']].values
                shi_ca = shi
            else:
                for x in ca:
                    mask_ca = mask_ca | (self.c_frame['type'] == x)
                pos_ca = self.c_frame[['x', 'y', 'z']].values[mask_ca,:]
                shi_ca = shi[mask_ca]

            mask_fa = False
            if fa == None:
                fa = ['all']
                pos_fa = self.c_frame[['x', 'y', 'z']].values
                shi_fa = shi
            else:
                for x in fa:
                    mask_fa = mask_fa | (self.c_frame['type'] == x)
                pos_fa = self.c_frame[['x', 'y', 'z']].values[mask_fa,:]
                shi_fa = shi[mask_fa]

            l = self.c_box[:, 1] - self.c_box[:, 0]

            hist, _ = np.histogram(
                [-1], range=[0, radius], bins=bins)
            g3 = np.zeros(hist.shape)
            for i,shi_ in zip(pos_ca,shi_ca):
                delta = np.abs(pos_fa - i)
                r2 = (np.minimum(delta, np.abs(delta - l))**2).sum(axis=1)
                mask = r2<(radius)**2
                inds = np.digitize(np.sqrt(r2[mask]), _)
                r2 = shi_fa[mask]*shi_
                for e_ind,r2_ in zip(inds,r2):
                    if (e_ind<len(hist)) and (e_ind>1):
                        g3[e_ind] +=  r2_

            rs = (_[1:] + _[:-1]) / 2
            if self.dim == 3:
                div = np.array(
                    [4 * np.pi * i**2 * r_step for i in rs])
                vol = l.prod()
            elif self.dim == 2:
                div = np.array(
                    [2 * np.pi * i * r_step for i in rs])
                vol = l[:-1].prod()
            else:
                raise()

            g3 = g3 / len(pos_ca) / div / len(pos_fa) * vol
            key += '{}-{}'.format(','.join([str(i) for i in ca]),
                                 ','.join([str(i) for i in fa])) + '_' + str(self.c_fn)
            self.gr_n[key] = [g3.copy(), rs.copy()]


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

            hist, _ = np.histogram(
                [-1], range=[0, self.radius], bins=self.r_bins)

            for i in pos_ca:
                delta = np.abs(pos_fa - i)
                r2 = (np.minimum(delta, np.abs(delta - l))**2).sum(axis=1)
                r2 = r2[r2<(self.radius)**2]
                h, _ = np.histogram(
                    np.sqrt(r2), range=[1e-5, self.radius], bins=self.r_bins)
                hist += h


            rs = (_[1:] + _[:-1]) / 2
            if self.dim == 3:
                div = np.array(
                    [4 * np.pi * i**2 * self.r_step for i in rs])
                vol = l.prod()
            elif self.dim == 2:
                div = np.array(
                    [2 * np.pi * i * self.r_step for i in rs])
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
                sq_list = np.zeros((int(Q/dq)))
                q_list = np.zeros((int(Q/dq)))
                for i in range(1,1+int(Q/dq)):
                    q_list[i-1] = dq*i
                    # F = sin(pi*r/R) / (pi*r/R)
                    # sq = integral [2*pi*r*(pdf-1)*sin(q*r)/q/r ] dr
                    F = np.sin(np.pi*r/R) / (np.pi*r/R)
                    sq_list[i-1] = 1+rho_*(2*np.pi*r*(p-1)*np.sin(dq*i*r)/(dq*i*r)*F*dr).sum()

                self.SQ[key] = [sq_list,q_list]

    def plot(self, key, which, *args, leg_prop={}, mean=False, lz_fit={}, ax=None, shift=0, **kwargs):
        which = self.__getattribute__(which)
        if key in which.keys() or mean:
            if ax is None:
                fig, ax = plt.panel(1, 1)
            if mean:
                p = 0
                ind = 0
                for a in which.keys():
                    if key in a:
                        p0, p1 = which[a]
                        p += p0
                        ind += 1
                p = p / ind
                x = list(which.values())[0][1]
                ax.plot(x, p+shift, *args, **kwargs)
            else:
                p, x = which[key]
                ax.plot(x, p+shift, *args, **kwargs)
        return p, x

    def CN(self, cutoff={}, default=2.0, frame_index=[0]):
        for f, box, fi in zip([self.frames[i] for i in frame_index],
                            [self.box[i] for i in frame_index], frame_index):
            print('Calculating nb_list for frame index ', fi)
            nb_list = {}
            L = box[:, 1] - box[:, 0]
            pos = f[['x', 'y', 'z']].values
            uniques = np.unique(f.type)
            self.bond_angles[fi] = {}
            for i in uniques:
                for j in uniques:
                    for k in uniques:
                        self.bond_angles[fi][i + 100 *
                                             j + 10000 * k] = [[], []]
            comb = list(itertools.combinations(uniques, 2))
            for u in uniques:
                comb.append((u, u))
            cut = {}
            for c in comb:
                cut['{}-{}'.format(*c)] = default
            cut.update(cutoff)
            types_list = []
            ids_list = []
            for ca, t, id in zip(pos, f.type, f.id):
                delta = np.abs(pos - ca)
                r2 = ((np.minimum(delta, np.abs(delta - L))**2).sum(axis=1))
                mask = (r2 < np.max(list(cut.values()) +
                                   [default])**2) & (r2 > 1.0e-3)
                r2 = r2[mask]
                ids = f.id[mask]
                types = f.type[mask]
                data = pos[mask]
                comb = ['{}-{}'.format(*np.sort([t, i])) for i in types]
                l1 = []
                l2 = []
                for x, y, z, tp in zip(r2, ids, comb, types):
                    if x < cut[z]**2:
                        l1.append(tp)
                        l2.append(y)
                types_list.append(l1)
                ids_list.append(l2)
                for s_atoms in itertools.combinations(zip(l1, l2), 2):
                    a = (data[ids == s_atoms[0][1]] - ca)[0]
                    mask = (a>L/2)
                    a[mask] -= L[mask]
                    mask = (a<-L/2)
                    a[mask] += L[mask]

                    b = (data[ids == s_atoms[1][1]] - ca)[0]
                    mask = (b>L/2)
                    b[mask] -= L[mask]
                    mask = (b<-L/2)
                    b[mask] += L[mask]

                    ct = (a * b.T).sum() / \
                        np.sqrt((a * a.T).sum() * (b * b.T).sum())
                    c2 = np.arccos(ct)
                    self.bond_angles[fi][s_atoms[0][0] + 100 *
                                         t + 10000 * s_atoms[1][0]][0] += [c2]
                    self.bond_angles[fi][s_atoms[0][0] + 100 * t + 10000 *
                                         s_atoms[1][0]][1] += [[s_atoms[0][1], id, s_atoms[1][1]]]
            self.nb_list[fi] = [f.id.values, f.type.values,
                                types_list, ids_list, len(ids_list)]
        return plt.gcf()

    def plot_bond_angles(self, key_list, *args, bins=360, ax=None, frame_index=[0], shift=0, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        hist = []
        for fi in frame_index:
            for key in key_list:
                y = self.bond_angles[fi][key[0] +
                                         key[1] * 100 + key[2] * 10000][0]
                y = 180 / np.pi * np.array(y)
                h, _ = np.histogram(y, range=[0, 180], bins=bins)
                h = h / h.sum() * 100
                x = (_[1:] + _[:-1]) / 2
                ax.plot(x, h+shift, *args, **kwargs)
                hist += [h]
        return hist, x
