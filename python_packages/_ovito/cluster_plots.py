from ovito.data import *
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
fig2 = plt.figure()
ax = fig.add_subplot(1,1,1)
ax2 = fig2.add_subplot(1,1,1)
fig.show()
ax.clear()
fig2.show()
ax2.clear()
counter = 0
max_size = 0
max_num = 0

num_clus = np.zeros((100,1))

def modify(frame, input, output):
    if frame == 0:
        max_num = 0
        max_size = 0
    global counter,max_size,max_num
    counter += 1
    if counter==10:
        counter = 0
        ax.clear()
    summ = []
    length =  max(output['Cluster'].marray)
    num_clus[frame] = length
    for i in range(1,length+1):
        summ.append(sum(output['Cluster'].marray==i))
    ax.plot(range(0,length),summ)
    max_size = max(max_size,max(summ))
    ax.set_ylim([0,max_size])
    fig.canvas.draw()
    ax2.clear()
    ax2.plot(range(0,100),num_clus)
    ax2.plot(frame,length,'-r*')
    max_num = max(max_num,length)
    ax2.set_ylim([0,max_num])
    fig2.canvas.draw()