"""
use matplotlib, mouse events, and Paths to define and match an extended source cut in one or more filters
after defining input files for each filter, this script pops open a matplotlib window which the user
can draw a region around the sources that should be included. the script then matches the two samples and writes
the combined sources to an output file.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path

def onclicki(event):
	global i_mag, i_fwhm, ix, iy
	if event.button == 1 :
		# print event.xdata, event.ydata
		ax.cla()
		ax.scatter(i0,i_fwhm, edgecolors='none')
		ax.set_ylim(0,15)
		ax.set_xlabel('magnitude')
		ax.set_ylabel('FWHM')
		ix.append(event.xdata)
		iy.append(event.ydata)
		# print ix,iy
		ax.plot(ix, iy ,'r-', lw=3, alpha=0.6)
	elif event.button == 3 :
		# print 'Removed last point from shape'
		ix.pop()
		iy.pop()
		ax.cla()
		ax.scatter(i0,i_fwhm, edgecolors='none')
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
	global g_mag, g_fwhm, gx, gy
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
		# print 'Removed last point from shape'
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
	
mkobsfile = 'ifirst_tol7.out'
ifwhm = 'getfwhm_i.log'
gfwhm = 'getfwhm_g.log'
title_string = 'AGC198606'

# get preliminary magnitude data which is close enough for this purpose
i0l = []
g0l = []
ifl = []
gfl = []
ixl = []
gxl = []
iyl = []
gyl = []

with open(mkobsfile,'r') as fin:
	lines = fin.readlines()
	for i, line in enumerate(lines):
		items = line.split()
		if len(items) > 0 :
			if items[1] == 'odi_i':
				i0l.append(float(items[6]))
				ixl.append(float(items[4]))
				iyl.append(float(items[5]))
			elif items[1] == 'odi_g':
				g0l.append(float(items[6]))
				gxl.append(float(items[4]))
				gyl.append(float(items[5]))
			else:
				print 'Ignored line',i,'in',mkobsfile,'during read'
		# else:
		# 	print 'Ignored line',i,'in',mkobsfile,'during read'

with open(ifwhm,'r') as fin :
	lines = fin.readlines()
	for i, line in enumerate(lines):
		items = line.split()
		if items[0] != '#':
			if items[12] != 'INDEF':
				ifl.append(float(items[12]))
			else:
				ifl.append(999.999)
		# else:
		# 	print 'Ignored line',i,'in', ifwhm, 'during read'
			
with open(gfwhm,'r') as fin :
	lines = fin.readlines()
	for i, line in enumerate(lines):
		items = line.split()
		if items[0] != '#':
			if items[12] != 'INDEF':
				gfl.append(float(items[12]))
			else:
				gfl.append(999.999)
		# else:
		# 	print 'Ignored line',i,'in', gfwhm, 'during read'

i0 = np.array(i0l)
g0 = np.array(g0l)
i_fwhm = np.array(ifl)
g_fwhm = np.array(gfl)

f_zp_g = open(title_string+'_g_phot.zp')
data_g = f_zp_g.read()
fl_g = data_g.split('\n', 1)[0]
zp_vals_g = fl_g.split()
cal_zp_g = float(zp_vals_g[2])
cal_color_g = float(zp_vals_g[0])
cal_zp_ge = float(zp_vals_g[3])
cal_color_ge = float(zp_vals_g[1])

f_zp_i = open(title_string+'_i_phot.zp')
data = f_zp_i.read()
fl_i = data.split('\n', 1)[1]
zp_vals_i = fl_i.split()
cal_zp_i = float(zp_vals_i[0])
cal_zp_ie = float(zp_vals_i[1])

# cal_zp_i = 25.0
# cal_zp_g = 25.0
# cal_color_g = 0.0

g_mag = g0 + cal_zp_g + cal_color_g*(g0-i0)
i_mag = i0 + cal_zp_i

print 'Recording region for extended source cut.'
print 'left click to select vertices, right click to delete last point'
print "press 'q' when finished"
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(i0,i_fwhm,edgecolors='none')
ax.set_ylim(0,15)
ax.set_xlabel('magnitude')
ax.set_ylabel('FWHM')

ix = []
iy = []
cid = fig.canvas.mpl_connect('button_press_event', onclicki)

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(g_mag,g_fwhm,edgecolors='none')
ax.set_ylim(0,15)
ax.set_xlabel('magnitude')
ax.set_ylabel('FWHM')

gx = []
gy = []
cid = fig.canvas.mpl_connect('button_press_event', onclickg)

plt.show()

iverts = zip(ix,iy)		# set up the Path necessary for testing membership
ipoints = zip(i0,i_fwhm)
i_filter = Path(iverts)
gverts = zip(gx,gy)		# set up the Path necessary for testing membership
g_filter = Path(gverts)
gpoints = zip(g_mag,g_fwhm)

stars_f_i = i_filter.contains_points(ipoints)
stars_f_g = g_filter.contains_points(gpoints)

print len([stars_f_i[i] for i in range(len(stars_f_i)) if (stars_f_i[i])]), 'stars in i'

match_f = list(stars_f_i)
for i in range(len(stars_f_i)) :
	match_f[i] = stars_f_i[i] & stars_f_g[i]

test = [match_f[i] for i in range(len(match_f)) if match_f[i]]

print len(test), 'stars matched'

escut_pos_file_g = open("escut_g.pos", 'w+')
escut_pos_file_i = open("escut_i.pos", 'w+')
for i in range(len(match_f)) :
	if match_f[i] == True :
		print >> escut_pos_file_g, gxl[i], gyl[i]
		print >> escut_pos_file_i, ixl[i], iyl[i]
escut_pos_file_g.close()
escut_pos_file_i.close()
