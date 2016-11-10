#! /usr/local/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.path import Path
from scipy import ndimage
from astropy.io import fits
from astropy import wcs
from photutils import detect_sources, source_properties, properties_table
from photutils.utils import random_cmap
from distfit import distfit
import aplpy
from scipy.stats import binned_statistic_2d
import matplotlib.patches as patches

hdulist = fits.open('agc249525_i_sh.fits')
hdu = hdulist[0]
naxis1, naxis2 = hdu.header['NAXIS1'], hdu.header['NAXIS2']
w = wcs.WCS(hdu.header)
ra_c0, dec_c0 = w.all_pix2world(0,0,1)
ra_cN, dec_cN = w.all_pix2world(naxis1,naxis2,1)
print ra_c0, dec_c0, ra_cN-ra_c0, dec_cN-dec_c0

mag_file = 'calibrated_mags.dat'

# read in magnitudes, colors, and positions(x,y)
gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
print len(gxr), "total stars"

# filter out the things with crappy color errors
color_error_cut = np.sqrt(2.0)*0.2
mag_error_cut = 0.2

gmi_errr = np.array([np.sqrt(g_ierrr[i]**2 + i_ierrr[i]**2) for i in range(len(gxr))])
gx = np.array([gxr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
gy = np.array([gyr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
g_mag = np.array([g_magr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
g_ierr = np.array([g_ierrr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
ix = np.array([ixr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
iy = np.array([iyr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
i_mag = np.array([i_magr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
i_ierr = np.array([i_ierrr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
gmi = np.array([gmir[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
gmi_err = np.array([np.sqrt(g_ierrr[i]**2 + i_ierrr[i]**2) for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])

print len(gx), "after color+mag error cut"
# nid = np.loadtxt(mag_file,usecols=(0,),dtype=int,unpack=True)
pixcrd = zip(ix,iy)

world = w.all_pix2world(pixcrd, 1)
ra_c, dec_c = w.all_pix2world(0,0,1)
# ra_c_d,dec_c_d = deg2HMS(ra=ra_c, dec=dec_c, round=True)
# print 'Corner RA:',ra_c_d,':: Corner Dec:',dec_c_d

# split the ra and dec out into individual arrays and transform to arcmin from the corner
i_ra = [abs((world[i,0]-ra_c)*60) for i in range(len(world[:,0]))]
i_dec = [abs((world[i,1]-dec_c)*60) for i in range(len(world[:,1]))]
# also preserve the decimal degrees for reference
i_rad = [world[i,0] for i in range(len(world[:,0]))]
i_decd = [world[i,1] for i in range(len(world[:,1]))]

bin_count, xedges, yedges, bin_id = binned_statistic_2d(i_rad, i_decd, gmi, statistic='count', bins=[5,5], range=[[ra_cN,ra_c0],[dec_c0,dec_cN]], expand_binnumbers=True)
xcenters= 0.5*(xedges[1:]+xedges[:-1])
ycenters= 0.5*(yedges[1:]+yedges[:-1])
# for i,b in enumerate(bin_id[0]):
#     print bin_id[0][i], bin_id[1][i], xedges[i], yedges[i], bin_count[i]

print bin_count, np.sum(bin_count), np.median(bin_count), np.std(bin_count)

fig=plt.figure(3)#, figsize=(22,17))  
fig.subplots_adjust(hspace=0.0, wspace=0.0)  
plt.clf()
for x in range(5):
    for y in range(5):
        bin_stars = np.where((bin_id[0] == x+1) & (bin_id[1] == y+1))
        # print bin_stars
        i_magBIN = i_mag[bin_stars]
        gmiBIN = gmi[bin_stars]
        
        with open('grid'+repr(x)+repr(y)+'.reg','w+') as grid:
            for i in range(len(ix[bin_stars])):
                # if gmiBIN[i] > 1 and gmiBIN[i] <1.5 and i_magBIN[i] > 22.:
                print >> grid, ix[bin_stars][i], iy[bin_stars][i]
        
        ax = plt.subplot2grid((5,5),(4-y,4-x))
        plt.scatter(gmiBIN, i_magBIN, color='black', marker='o', s=3, edgecolors='none')
        plt.text(3.45, 16.,'$\\alpha$ = {0:6.3f}'.format(xcenters[x]), horizontalalignment='right', verticalalignment='center', fontsize=4)
        plt.text(3.45, 17.,'$\delta$ = {0:6.3f}'.format(ycenters[y]), horizontalalignment='right', verticalalignment='center', fontsize=4)
        plt.text(3.45, 18.,'$N$ = {0:4d}'.format(int(bin_count[x,y])), horizontalalignment='right', verticalalignment='center', fontsize=4)
        if x==4 and y==2:
            plt.ylabel('$i$')
        if x==2 and y==0:    
            plt.xlabel('$(g-i)$')
        plt.ylim(27,15)
        plt.xlim(-1,3.75)
        ax.xaxis.set_ticks([0,1,2,3])
        if y>0:
            ax.xaxis.set_ticklabels([])
        if x <4 and y>0:
            ax.yaxis.set_ticklabels([])
# plt.tight_layout()
# plt.show()
plt.savefig('cmd_grid.pdf')