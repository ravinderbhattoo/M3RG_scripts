from kurma.structure_order import pdf
import numpy as np
import importlib
import matplotlib.pyplot as plt
from shadow.plot import xlabel, ylabel, xlim, ylim, panel, xticks, set_things
set_things()

P = pdf.PDF('./structure_order/0.1_300.dump', dim=2, fn=[1])

P.make_pdf(atoms=[[1], [2]])
P.make_pdf(atoms=[None, None])

leg_prop = {'loc': 1, 'bbox_to_anchor': [0.9, 0.9]}
fig, [ax] = panel(1,1)
P.plot('1-2', 'PDF', 'b-', label='1-2', leg_prop=leg_prop, ax=ax, mean=True, shift=5)
P.plot('all-all', 'PDF', 'r-', label='all', leg_prop=leg_prop, ax=ax, mean=True,)

xlabel(r'Distance, r (\AA)')
ylabel(r'Pair distribution function (PDF)')
xlim([0,30])
plt.legend()


leg_prop = {'loc': 1, 'bbox_to_anchor': [0.9, 0.9]}
fig, [ax] = panel(1,1)
P.plot('1-2', 'TCF', 'b-', label='1-2', leg_prop=leg_prop, ax=ax, mean=True, shift=5)
P.plot('all-all', 'TCF', 'r-', label='all', leg_prop=leg_prop, ax=ax, mean=True,)

xlabel(r'Distance, r (\AA)')
ylabel(r'Total correlation function (TCF)')
xlim([0,30])
plt.legend()


leg_prop = {'loc': 1, 'bbox_to_anchor': [0.9, 0.9]}
fig, [ax] = panel(1,1)
P.plot('1-2', 'SQ', 'b-', label='1-2', leg_prop=leg_prop, ax=ax, mean=True, shift=5)
P.plot('all-all', 'SQ', 'r-', label='all', leg_prop=leg_prop, ax=ax, mean=True,)

xlabel(r'Wave-vector, k (\AA$_{-1}$)')
ylabel(r'Structure factor (SQ)')
xlim([0,30])
plt.legend()



cutoff={'1-2':1.8}
c = P.CN(cutoff, frame_index=[0], default=1.8)

leg_prop = {'loc': 1, 'bbox_to_anchor': [0.9, 0.9]}
fig, [ax] = panel(1,1)
P.plot_bond_angles([[1,2,1],[2,1,2]],label=r'\textbf{1-2-1}')
xticks(range(0,190,30))
xlabel(r'Bond angle, (\AA)')
ylabel(r'Frequency')
plt.legend()
