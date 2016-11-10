#! /usr/local/bin/python
import os, sys, getopt, warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.path import Path
from matplotlib import cm
from astropy import wcs
from astropy.io import fits
from pyraf import iraf
try :
	from scipy import ndimage
except ImportError :
	from stsci import ndimage

def main(argv):
	fwhm = 2.0
	filter_string = 'psmap'
	fwhm_string = '2.0'
	iraf.tv(_doprint=0)
	home_root = os.environ['HOME']
	imexam_flag = False

	# print "Getting fits files..."
	# Load the FITS header using astropy.io.fits
	for file_ in os.listdir("./"):
		if file_.endswith("i_sh.fits"):
			fits_file_i = file_

 	for file_ in os.listdir("./"):
 		if file_.endswith("g_sh.fits"):
 			fits_file_g = file_

	fits_i = fits.open(fits_file_i)
	fits_g = fits.open(fits_file_g)
	# print "Opened fits files:",fits_file_g,"&",fits_file_i

	# objid = fits_i[0].header['OBJECT']
	title_string = fits_file_i[0:9]		# get the first 9 characters of the filename.

	# set up some filenames
	out_file = filter_string + '_' + fwhm_string + '_' + title_string + '.pdf'
	mag_file = 'calibrated_mags.dat'
	mark_file = 'f_list_' + filter_string + '_' + fwhm_string + '_' + title_string + '.dat'
	circ_file = 'c_list_' + filter_string + '_' + fwhm_string + '_' + title_string + '.dat'


	# read in magnitudes, colors, and positions(x,y)
	gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
	print len(gxr), "total stars"

	# filter out the things with crappy color errors
	color_error_cut = np.sqrt(2.0)*0.2
	mag_error_cut = 0.2

	gmi_errr = [np.sqrt(g_ierrr[i]**2 + i_ierrr[i]**2) for i in range(len(gxr))]
	gx = [gxr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	gy = [gyr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	g_mag = [g_magr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	g_ierr = [g_ierrr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	ix = [ixr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	iy = [iyr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	i_mag = [i_magr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	i_ierr = [i_ierrr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	gmi = [gmir[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
	gmi_err = [np.sqrt(g_ierrr[i]**2 + i_ierrr[i]**2) for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]

	# nid = np.loadtxt(mag_file,usecols=(0,),dtype=int,unpack=True)
	pixcrd = zip(ix,iy)


	# print "Reading WCS info from image header..."
	# Parse the WCS keywords in the primary HDU
	warnings.filterwarnings('ignore', category=UserWarning, append=True)
	w = wcs.WCS(fits_i[0].header)

	# Print out the "name" of the WCS, as defined in the FITS header
	# print w.wcs.name

	# Print out all of the settings that were parsed from the header
	# w.wcs.print_contents()

	# Convert pixel coordinates to world coordinates
	# The second argument is "origin" -- in this case we're declaring we
	# have 1-based (Fortran-like) coordinates.
	world = w.all_pix2world(pixcrd, 1)
	ra_c, dec_c = w.all_pix2world(0,0,1)
	ra_c_d,dec_c_d = deg2HMS(ra=ra_c, dec=dec_c, round=True)
	print 'Corner RA:',ra_c_d,':: Corner Dec:',dec_c_d

	fwhm_i = fits_i[0].header['FWHMPSF']
	fwhm_g = fits_g[0].header['FWHMPSF']

	print 'Image FWHM :: g = {0:5.3f} : i = {1:5.3f}'.format(fwhm_g,fwhm_i)

	fits_i.close()
	fits_g.close()

	# split the ra and dec out into individual arrays and transform to arcmin from the corner
	i_ra = [abs((world[i,0]-ra_c)*60) for i in range(len(world[:,0]))]
	i_dec = [abs((world[i,1]-dec_c)*60) for i in range(len(world[:,1]))]
	# also preserve the decimal degrees for reference
	i_rad = [world[i,0] for i in range(len(world[:,0]))]
	i_decd = [world[i,1] for i in range(len(world[:,1]))]

	# bin the filtered stars into a grid with pixel size XXX
	# print "Binning for m-M =",dm
	bins = 165
	width = 22 

	grid, xedges, yedges = np.histogram2d(i_dec, i_ra, bins=[bins,bins], range=[[0,width],[0,width]])
	hist_points = zip(xedges,yedges)

	sig = ((bins/width)*fwhm)/2.355
	pltsig = fwhm/2.0

	# convolve the grid with a gaussian
	# print "Convolving for m-M =",dm
	grid_gaus = ndimage.filters.gaussian_filter(grid, sig, mode='constant', cval=0)
	S = np.array(grid_gaus*0)
	S_th = 3.0

	# print grid_gaus[0:44][0:44]

	grid_mean = np.mean(grid_gaus)
	grid_sigma = np.std(grid_gaus)
	S = (grid_gaus-grid_mean)/grid_sigma

	above_th = [(int(i),int(j)) for i in range(len(S)) for j in range(len(S[i])) if (S[i][j] >= S_th)]

	# for i in range(143) :
	# 	for j in range(143) :
	# 		ie = i + 22
	# 		je = j + 22
	# 		grid_mean = np.mean(grid_gaus[i:ie,j:je])
	# 		grid_sigma = np.std(grid_gaus[i:ie,j:je])
	# 		S[(11+i),(11+j)] = (grid_gaus[(11+i),(11+j)]-grid_mean)/grid_sigma
	# 		# print i,j,grid_mean,grid_sigma,S[(22+i),(22+j)]

	# find the maximum point in the grid and center the circle there
	x_cent, y_cent = np.unravel_index(grid_gaus.argmax(),grid_gaus.shape)
	print 'Max of gaussian convolved grid located at:','('+'{0:6.3f}'.format(yedges[y_cent])+','+'{0:6.3f}'.format(xedges[x_cent])+')'
	# print grid_gaus[x_cent][y_cent]
	print 'Value of S at above:','{0:6.3f}'.format(S[x_cent][y_cent])

	x_cent_S, y_cent_S = np.unravel_index(S.argmax(),S.shape)
	print 'Max of S located at:','('+'{0:6.3f}'.format(yedges[y_cent_S])+','+'{0:6.3f}'.format(xedges[x_cent_S])+')'
	# print grid_gaus[x_cent][y_cent]
	print 'Value of S at above:','{0:6.3f}'.format(S[x_cent_S][y_cent_S])

	print 'Number of bins above S_th: {0:4d}'.format(len(above_th))

	
	# make a circle to highlight a certain region
	cosd = lambda x : np.cos(np.deg2rad(x))
	sind = lambda x : np.sin(np.deg2rad(x))
	x_circ = [yedges[y_cent] + pltsig*cosd(t) for t in range(0,359,1)]
	y_circ = [xedges[x_cent] + pltsig*sind(t) for t in range(0,359,1)]
	xy_points = zip(i_ra,i_dec)
	verts_circ = zip(x_circ,y_circ)
	circ_filter = Path(verts_circ)

	stars_circ = circ_filter.contains_points(xy_points)	

	i_mag_c = [i_mag[i] for i in range(len(i_mag)) if (stars_circ[i])]
	gmi_c = [gmi[i] for i in range(len(i_mag)) if (stars_circ[i])]
	i_ra_c = [i_ra[i] for i in range(len(i_mag)) if (stars_circ[i])]
	i_dec_c = [i_dec[i] for i in range(len(i_mag)) if (stars_circ[i])]
	i_rad_c = [i_rad[i] for i in range(len(i_mag)) if (stars_circ[i])]
	i_decd_c = [i_decd[i] for i in range(len(i_mag)) if (stars_circ[i])]
	i_x_c = [ix[i] for i in range(len(i_mag)) if (stars_circ[i])]
	i_y_c = [iy[i] for i in range(len(i_mag)) if (stars_circ[i])]

	f2 = open(circ_file, 'w+')
	for i in range(len(i_x_c)) :
		print >> f2, '{0:8.2f} {1:8.2f} {2:12.8f} {3:12.8f} {4:8.2f} '.format(i_x_c[i],i_y_c[i],i_rad_c[i],i_decd_c[i],i_mag_c[i],gmi_c[i])
	f2.close()

	# print "{0:3d} stars in filter, {1:3d} stars in circle, {2:3d} stars in both.".format(len(i_mag_f),len(i_mag_c),len(i_mag_fc))
	# for i in range(len(i_x_fc)) :
	# 	print (i_x_fc[i],i_y_fc[i])

	circ_c_x = (yedges[y_cent]/60.)+ra_c
	circ_c_y = (xedges[x_cent]/60.)+dec_c
	circ_pix_x, circ_pix_y = w.wcs_world2pix(circ_c_x,circ_c_y,1)

	hi_c_ra, hi_c_dec = 142.5104167, 16.6355556
	hi_c_x, hi_c_y = abs((hi_c_ra-ra_c)*60), abs((hi_c_dec-dec_c)*60)
	hi_x_circ = [hi_c_x + pltsig*cosd(t) for t in range(0,359,1)]
	hi_y_circ = [hi_c_y + pltsig*sind(t) for t in range(0,359,1)]

	hi_pix_x,hi_pix_y = w.wcs_world2pix(hi_c_ra,hi_c_dec,1)

	center_file = 'im_cens_'+fwhm_string+'.dat'
	im_cens = open(center_file,'w+')
	print >> im_cens, -circ_pix_x, circ_pix_y, 'circle_center'
	print >> im_cens, hi_pix_x, hi_pix_y, 'HI_centroid'
	im_cens.close()

	mark_radius = int(pltsig*60/0.11)
	mark_radii = str(repr(mark_radius-2)+','+repr(mark_radius-1)+','+repr(mark_radius)+','+repr(mark_radius+1)+','+repr(mark_radius+2))

	# if disp_flag :
	# 	iraf.unlearn(iraf.tv.tvmark)
	# 	iraf.tv.tvmark.setParam('label',"no")
	# 	iraf.tv.tvmark.setParam('pointsize',7)
	# 	iraf.tv.tvmark.setParam('mark',"circle")
	# 	iraf.tv.display(image=fits_file_g, frame=1)
	# 	iraf.tv.display(image=fits_file_i, frame=2)
	# 	iraf.tv.tvmark(frame=1, coords=mark_file, radii="14,15,16", color=207)
	# 	iraf.tv.tvmark(frame=1, coords=circ_file, radii="20,21,22", color=209)
	# 	iraf.tv.tvmark(frame=1, coords=center_file, txsize=4, mark="plus", color=208, label="yes")
	# 	iraf.tv.tvmark(frame=1, coords=center_file, radii=mark_radii, color=208)
	# 	iraf.tv.tvmark(frame=2, coords=mark_file, radii="14,15,16", color=207, label="yes")
	# 	iraf.tv.tvmark(frame=2, coords=circ_file, radii="20,21,22", color=209)
	# 	iraf.tv.tvmark(frame=2, coords=center_file, txsize=4, mark="plus", color=208, label="yes")
	# 	iraf.tv.tvmark(frame=2, coords=center_file, radii=mark_radii, color=208)
	# 
	# setup the pdf output
	pp = PdfPages('psmap_'+ out_file)
	# plot
	# print "Plotting for m-M = ",dm
	ax0 = plt.subplot(2,2,1)
	plt.scatter(i_ra, i_dec,  color=gmi, marker='o', s=(27-np.array(i_mag)), edgecolors='none', cmap=cm.rainbow)
	plt.plot(x_circ,y_circ,linestyle='-', color='magenta')
	plt.plot(hi_x_circ,hi_y_circ,linestyle='-', color='black')
	plt.colorbar()
	# plt.scatter(i_ra_c, i_dec_c,  color='red', marker='o', s=3, edgecolors='none')
	plt.ylabel('Dec (arcmin)')
	plt.xlim(0,max(i_ra))
	plt.ylim(0,max(i_dec))
	plt.title('sky positions')
	ax0.set_aspect('equal')

	ax1 = plt.subplot(2,2,2)
	plt.scatter(gmi, i_mag,  color='black', marker='o', s=1, edgecolors='none')
	plt.scatter(gmi_c, i_mag_c,  color='red', marker='o', s=3, edgecolors='none')
	plt.tick_params(axis='y',left='on',right='off',labelleft='on',labelright='off')
	ax1.yaxis.set_label_position('left')
	plt.ylabel('$i$')
	plt.ylim(25,15)
	plt.xlim(-1,4)
	ax1.set_aspect(0.5)

	ax2 = plt.subplot(2,2,3)

	extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
	plt.imshow(S, extent=extent, interpolation='nearest', cmap=cm.cubehelix)
	cbar_S = plt.colorbar()
	# cbar_S.tick_params(labelsize=10)
	plt.plot(x_circ,y_circ,linestyle='-', color='magenta')
	plt.plot(hi_x_circ,hi_y_circ,linestyle='-', color='black')
	# X, Y = np.meshgrid(xedges,yedges)
	# ax3.pcolormesh(X,Y,grid_gaus)
	plt.xlabel('RA (arcmin)')
	plt.ylabel('Dec (arcmin)')

	# plt.ylabel('Dec (arcmin)')
	plt.xlim(0,max(i_ra))
	plt.ylim(0,max(i_dec))
	ax2.set_aspect('equal')

	ax3 = plt.subplot(2,2,4)
	# plt.imshow(imagedata)
	# extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
	# plt.imshow(grid, extent=extent, interpolation='nearest')
	# cbar_gaus = plt.colorbar()
	# # cbar_gaus.tick_params(labelsize=10)
	# plt.plot(x_circ,y_circ,linestyle='-', color='magenta')
	# plt.plot(hi_x_circ,hi_y_circ,linestyle='-', color='black')
	# plt.xlabel('RA (arcmin)')
	# plt.ylabel('Dec (arcmin)')
	# plt.xlim(0,max(i_ra))
	# plt.ylim(0,max(i_dec))
	plt.scatter(gmi_c, i_mag_c,  color='red', marker='o', s=3, edgecolors='none')
	plt.tick_params(axis='y',left='on',right='off',labelleft='on',labelright='off')
	ax0.yaxis.set_label_position('left')
	plt.xlabel('$(g-i)$')
	plt.ylabel('$i$')
	plt.ylim(25,15)
	plt.xlim(-1,4)
	ax3.set_aspect(0.5)	


	plt.suptitle(title_string+ ' -- ' + fwhm_string + ' arcmin smoothing')
	plt.show

	pp.savefig()
	pp.close()

	# ax0 = plt.subplot(1,1,1)
	# 
	# pp.savefig()
	# 
	# pp.close()

	if imexam_flag :
		iraf.unlearn(iraf.tv.imexamine, iraf.rimexam)
		iraf.tv.rimexam.setParam('radius',int(fwhm_i))
		iraf.tv.rimexam.setParam('rplot',12.)
		iraf.tv.rimexam.setParam('fittype','gaussian')
		iraf.tv.rimexam.setParam('center','yes')
		iraf.tv.imexamine(input=fits_file_i, frame=2)

	# clean up
	if os.path.isfile(center_file) :
		os.remove(center_file)
	if os.path.isfile(mark_file) :
		os.remove(mark_file)
	if os.path.isfile(circ_file) :
		os.remove(circ_file)

def deg2HMS(ra='', dec='', round=False):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
	 if str(dec)[0] == '-':
		ds, dec = '-', abs(dec)
	 deg = int(dec)
	 decM = abs(int((dec-deg)*60))
	 if round:
		decS = int((abs((dec-deg)*60)-decM)*60)
	 else:
		decS = (abs((dec-deg)*60)-decM)*60
	 DEC = '{0}{1}:{2}:{3}'.format(ds, deg, decM, decS)

  if ra:
	 if str(ra)[0] == '-':
		rs, ra = '-', abs(ra)
	 raH = int(ra/15)
	 raM = int(((ra/15)-raH)*60)
	 if round:
		raS = int(((((ra/15)-raH)*60)-raM)*60)
	 else:
		raS = ((((ra/15)-raH)*60)-raM)*60
	 RA = '{0}{1}:{2}:{3}'.format(rs, raH, raM, raS)

  if ra and dec:
	 return (RA, DEC)
  else:
	 return RA or DEC

if __name__ == "__main__":
	main(sys.argv[1:])	
