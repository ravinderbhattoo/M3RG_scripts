import ovito
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import numpy as np


# This user-defined function is called by OVITO to let it draw arbitrary graphics on top of the viewport.
# It is passed a QPainter (see http://qt-project.org/doc/qt-5/qpainter.html).

def render(painter, **args):
	ovito.dataset.selected_node.compute()
	output = ovito.dataset.selected_node.output
	heightWindow = painter.window().height()
	widthWindow = painter.window().width()
	x0 = widthWindow/2
	y0 = heightWindow/2

	painter.setPen(QPen(Qt.green,  4, Qt.DashLine))
	r= 200		
	painter.drawEllipse(x0-r,y0-r,2*r,2*r)
	#painter.drawEllipse(x0-10,y0-10,20,20)
		
	ids = output['Particle Identifier'].marray
	pos = output['Position'].marray
	bonds = output['Bonds'].array
	
	mask = ids==3091
	id = np.argmax(mask)
	bs = bonds[bonds[:,0]==id,:]
	
	for b in bs:
		c = pos[b[0]]
		o = pos[b[1]]
		dx = o[0]-c[0]
		dy = o[1]-c[1]
		l = np.sqrt(dx**2 + dy**2)
		sl = l*60
		theta = np.arctan2(dy,dx)
		painter.drawLine(x0,y0,x0+sl*np.cos(theta),y0-sl*np.sin(theta))
	
	





	
	# This demo code prints the current animation frame into the upper left corner of the viewport.
	text1 = "Frame {}".format(ovito.dataset.anim.current_frame)
	#painter.drawText(10, 10 + painter.fontMetrics().ascent(), text1)

	# Also print the current number of particles into the lower left corner of the viewport.
	node = ovito.dataset.selected_node
	num_particles = (node.compute().number_of_particles if node else 0)
	text2 = "{} particles".format(num_particles)
	#painter.drawText(10, painter.window().height() - 10, text2)
	