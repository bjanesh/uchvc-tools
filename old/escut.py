#!/usr/local/python
"""
use matplotlib, mouse events, and Paths to define and match an extended source cut in one or more filters
after defining input files for each filter, this script pops open a matplotlib window which the user
can draw a region around the sources that should be included. the script then matches the two samples and writes
the combined sources to an output file.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
def escut(fwhmlog_g='getfwhm_g.log', fwhmlog_i='getfwhm_i.log'):
	# do an extended source cut to get rid of background galaxies/etc. subject to change
	if not os.path.isfile('escut_i.pos') :
		def onclicki(event):
			global fig, ax,i_mag, i_fwhm, ix, iy
			if event.button ==1 :
				# print event.xdata, event.ydata
				ax.cla()
				ax.scatter(i_mag,i_fwhm, edgecolors='none')
				ax.set_ylim(0,15)
				ax.set_xlabel('magnitude')
				ax.set_ylabel('FWHM')
				ix.append(event.xdata)
				iy.append(event.ydata)
			# print ix,iy
				ax.plot(ix, iy ,'r-', lw=3, alpha=0.6)
			elif event.button == 3 :
				#print 'Removedlast point from shape'
				ix.pop()
				iy.pop()
				ax.cla()
				ax.scatter(i_mag,i_fwhm, edgecolors='none')
				ax.set_ylim(0,15)
				ax.set_xlabel('magnitude')
				ax.set_ylabel('FWHM')
				ax.plot(ix, iy ,'r-', lw=3, alpha=0.6)
			elif event.key == 'q' :
				print 'Finished recording points'
				fig.canvas.mpl_disconnect(cid)
			fig.canvas.draw()
			return True
	
		def onclickg(event):
			global fig, ax, g_mag, g_fwhm, gx, gy
			if event.button == 1 :
			# print event.xdata, event.ydata
				ax.cla()
				ax.scatter(g_mag,g_fwhm, edgecolors='none')
				ax.set_ylim(0,15)
				ax.set_xlabel('magnitude')
				ax.set_ylabel('FWHM')
				gx.append(event.xdata)
				gy.append(event.ydata)
				ax.plot(gx, gy ,'r-', lw=3, alpha=0.6)
			elif event.button == 3 :
				#print 'Removedlast point from shape'
				gx.pop()
				gy.pop()
				ax.cla()
				ax.scatter(g_mag,g_fwhm, edgecolors='none')
				ax.set_ylim(0,15)
				ax.set_xlabel('magnitude')
				ax.set_ylabel('FWHM')
				ax.plot(gx, gy ,'r-', lw=3, alpha=0.6)
			elif event.key == 'q' :
				print 'Finished recording points'
				fig.canvas.mpl_disconnect(cid)
			fig.canvas.draw()
			return True
	
		# get preliminary magnitude data which is close enough for this purpose
		i0l = []
		g0l = []
		ifl = []
		gfl = []
		ixl = []
		gxl = []
		iyl =[]
		gyl = []
	
		with open(fwhmlog_i,'r') as fin :
			lines = fin.readlines()
			for i, line in enumerate(lines):
				items =line.split()
				if items[0] != '#':
					if items[12] != 'INDEF':
						ixl.append(float(items[0]))
						iyl.append(float(items[1]))
						i0l.append(float(items[5]))
						ifl.append(float(items[12]))
					else:
						ixl.append(float(items[0]))
						iyl.append(float(items[1]))
						i0l.append(float(items[5]))
						ifl.append(999.999)
				else:
				    print 'Ignored line',i,'in', fwhmlog_i, 'during read'
	
		with open(fwhmlog_g,'r') as fin :
			lines = fin.readlines()
			for i, line in enumerate(lines):
				items = line.split()
				if items[0] != '#':
					if items[12] != 'INDEF':
						gxl.append(float(items[0]))
						gyl.append(float(items[1]))
						g0l.append(float(items[5]))
						gfl.append(float(items[12]))
					else:
						gxl.append(float(items[0]))
						gyl.append(float(items[1]))
						g0l.append(float(items[5]))
						gfl.append(999.999)
				else:
				    print 'Ignored line',i,'in', fwhmlog_g, 'during read'
	
		i_i = np.array(i0l)
		g_i = np.array(g0l)
		i_fwhm = np.array(ifl)
		g_fwhm = np.array(gfl)
	
		i_mag = i_i
		g_mag = g_i
	
		if not os.path.isfile('escutRegion_i.txt'):
			print 'Recording region for I-BAND extended source cut.'
			print 'left click to select vertices, right click to delete last point'
			print "press 'q' when finished"
			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.scatter(i_mag,i_fwhm,edgecolors='none')
			ax.set_ylim(0,15)
			ax.set_xlabel('magnitude')
			ax.set_ylabel('FWHM')
	
			ix = []
			iy = []
			cid = fig.canvas.mpl_connect('button_press_event', onclicki)
	
			plt.show()
	
			with open('escutRegion_i.txt','w+') as fout:
				for i in range(len(ix)):
					print >> fout, ix[i], iy[i]
	
		else :
			print 'Using prerecorded escut file for i-band'
			ix, iy = np.loadtxt('escutRegion_i.txt',usecols=(0,1),unpack=True)
	
	
		if not os.path.isfile('escutRegion_g.txt'):
			print 'Recording region for G-BAND extended source cut.'
			print 'left click to select vertices, right click to delete last point'
			print "press 'q' when finished"
			fig =plt.figure()
			ax = fig.add_subplot(111)
			ax.scatter(g_mag,g_fwhm,edgecolors='none')
			ax.set_ylim(0,15)
			ax.set_xlabel('magnitude')
			ax.set_ylabel('FWHM')
	
			gx = []
			gy = []
			cid = fig.canvas.mpl_connect('button_press_event', onclickg)
	
			plt.show()
	
			with open('escutRegion_g.txt','w+') as fout:
				for i in range(len(gx)):
					print >> fout, gx[i], gy[i]
		else :
			print 'Using prerecorded escut file for g-band'
			gx, gy = np.loadtxt('escutRegion_g.txt',usecols=(0,1),unpack=True)
	
		iverts = zip(ix,iy)        # set up the Path necessary for testing membership
		ipoints = zip(i_mag,i_fwhm)
		i_filter = Path(iverts)
		gverts = zip(gx,gy)        # set up the Path necessary for testing membership
		g_filter = Path(gverts)
		gpoints = zip(g_mag,g_fwhm)
	
		stars_f_i = i_filter.contains_points(ipoints)
		stars_f_g = g_filter.contains_points(gpoints)
	
		# escutFlag2 = np.genfromtxt('escutFlag2.txt', dtype=bool, usecols=(0,),unpack=True)
	
	
		match_f = list(stars_f_i)
		match_f[i] = stars_f_i[i] and stars_f_g[i] #and escutFlag2[i]
		test = [match_f[i] for i in range(len(match_f)) if (match_f[i])]
	
		print len(test), 'stars matched'

		escut_pos_file_g = open("escut_g.pos", 'w+')
		escut_pos_file_i = open("escut_i.pos", 'w+')
		for i in range(len(match_f)) :
			if match_f[i] == True :
				print >> escut_pos_file_g, gxl[i], gyl[i]
				print >> escut_pos_file_i, ixl[i], iyl[i]
		escut_pos_file_g.close()
		escut_pos_file_i.close()        

