#! /usr/local/bin/python
# -*- coding: utf-8 -*-
import os, sys, getopt, warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.path import Path
from matplotlib import cm
from astropy import wcs
from astropy.io import fits
from pyraf import iraf
import scipy.stats as ss
from photutils import detect_sources, source_properties, properties_table
from photutils.utils import random_cmap
try :
    from scipy import ndimage
except ImportError :
    print 'bad import'

def getHIcentroid(object):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    uchvcdb = os.environ['HOME']+'/projects/uchvc-db/predblist.sort.csv'
    name, coords = np.loadtxt(uchvcdb, usecols=(1,2), dtype=str, delimiter=',', unpack=True)
    # find the right row
    coord = [this for i,this in enumerate(coords) if object in name[i]][0]

    # parse the coordinate into a better string
    rah = coord[0:2]
    ram = coord[2:4]
    ras = coord[4:8]
    ded = coord[8:11]
    dem = coord[11:13]
    des = coord[13:15]
    
    ra = rah+':'+ram+':'+ras
    dec = ded+':'+dem+':'+des
    coord_hi = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    ra_hi = coord_hi.ra
    dec_hi = coord_hi.dec
    return ra_hi, dec_hi
    
def dist2HIcentroid(ra, dec, ra_hi, dec_hi):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    c_hi = SkyCoord(ra = ra_hi, dec = dec_hi, unit=(u.hourangle, u.deg))
    c_peak = SkyCoord(ra = ra, dec = dec, unit=(u.hourangle, u.deg))
    sep = c_hi.separation(c_peak)
    return sep.arcsecond
    
def grid_smooth(i_ra_f, i_dec_f, fwhm, width, height):
    # bin the filtered stars into a grid with pixel size XXX
    # print "Binning for m-M =",dm
    # bins = 165
    # width = 30
    bins_h = int(height * 60. / 8.)
    bins_w = int(width * 60. / 8.)

    grid, xedges, yedges = np.histogram2d(i_dec_f, i_ra_f, bins=[bins_h,bins_w], range=[[0,height],[0,width]])
    hist_points = zip(xedges,yedges)

    sig = ((bins_w/width)*fwhm)/2.355
    pltsig = fwhm/2.0

    # convolve the grid with a gaussian
    grid_gaus = ndimage.filters.gaussian_filter(grid, sig, mode='constant', cval=0)
    S = np.array(grid_gaus*0)
    S_th = 3.0

    grid_mean = np.mean(grid_gaus)
    grid_sigma = np.std(grid_gaus)
    S = (grid_gaus-grid_mean)/grid_sigma

    above_th = [(int(i),int(j)) for i in range(len(S)) for j in range(len(S[i])) if (S[i][j] >= S_th)]

    segm = detect_sources(S, 2.0, npixels=5)
    props = source_properties(S, segm)
    columns = ['id', 'maxval_xpos', 'maxval_ypos', 'max_value', 'area']
    tbl = properties_table(props, columns=columns)
    # print tbl
    # rand_cmap = random_cmap(segm.max + 1, random_state=12345)

    # find the maximum point in the grid and center the circle there
    x_cent, y_cent = np.unravel_index(grid_gaus.argmax(),grid_gaus.shape)
    x_cent_S, y_cent_S = np.unravel_index(S.argmax(),S.shape)
    # print 'Max of S located at:','('+'{0:6.3f}'.format(y_cent_S)+','+'{0:6.3f}'.format(x_cent_S)+')'
    # print 'Value of S at above:','{0:6.3f}'.format(S[x_cent_S][y_cent_S])
    # print 'Number of bins above S_th: {0:4d}'.format(len(above_th))
    return xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl, segm
    
def distfit(n,dists,title,width,height,fwhm,dm,samples=1000):
    from scipy.stats import lognorm

    bins_h = int(height * 60. / 8.)
    bins_w = int(width * 60. / 8.)
    sig = ((bins_w/width)*fwhm)/2.355
    valsLP = []
    for i in range(samples) :
        random_ra = width*np.random.random_sample((n,))
        random_dec = height*np.random.random_sample((n,))
        random_xy = zip(random_ra,random_dec)
        grid_r, xedges_r, yedges_r = np.histogram2d(random_dec, random_ra, bins=[bins_h,bins_w], range=[[0,height],[0,width]])
        hist_points_r = zip(xedges_r,yedges_r)
        grid_gaus_r = ndimage.filters.gaussian_filter(grid_r, sig, mode='constant', cval=0)
        S_r = np.array(grid_gaus_r*0)

        grid_mean_r = np.mean(grid_gaus_r)
        grid_sigma_r = np.std(grid_gaus_r)
        S_r = (grid_gaus_r-grid_mean_r)/grid_sigma_r

        x_cent_r, y_cent_r = np.unravel_index(grid_gaus_r.argmax(),grid_gaus_r.shape)
        valsLP.append(S_r[x_cent_r][y_cent_r])

    x = np.linspace(2, 22, 4000)

    bins, edges = np.histogram(valsLP, bins=400, range=[2,22], normed=True)
    centers = (edges[:-1] + edges[1:])/2.

    al,loc,beta=lognorm.fit(valsLP)
    pct = 100.0*lognorm.cdf(dists, al, loc=loc, scale=beta)
    print 'Significance of detection:','{0:6.3f}%'.format(pct)
    
def galmap(image):
    # read in the results
    ras,decs,psfMag_u,psfMagErr_u,psfMag_g,psfMagErr_g,psfMag_r,psfMagErr_r,psfMag_i,psfMagErr_i,psfMag_z,psfMagErr_z = np.loadtxt(image[:-5]+'.sdss',usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=2)
    probPSF = np.loadtxt(image[:-5]+'.sdss', usecols=(12,), dtype=int, unpack=True, delimiter=',', skiprows=2)
    # keep things that are actually stars (defined as being psf's) and with the right magnitude range (arbitrary)
    keep_gals = (probPSF == 0)
    print 'keeping', len(np.where(keep_gals)[0]), 'galaxies of', len(psfMag_g), 'sources'
    
    fits_g = fits.open(image)
    w = wcs.WCS(fits_g[0].header)
    footprint = w.calc_footprint()
    se_corner = footprint[0]
    ne_corner = footprint[1]
    nw_corner = footprint[2]
    sw_corner = footprint[3]
    # print se_corner, ne_corner, nw_corner, sw_corner
    width = (ne_corner[0]-nw_corner[0])*60.
    height = (ne_corner[1]-se_corner[1])*60.
    print width, height
    ra_corner, dec_corner = w.all_pix2world(0,0,1)
    
    coords2 = zip(ras,decs)
    
    # split the ra and dec out into individual arrays and transform to arcmin from the corner
    i_ra_f = [abs((ras[i]-ra_corner)*60) for i in range(len(ras)) if keep_gals[i]]
    i_dec_f = [abs((decs[i]-dec_corner)*60) for i in range(len(decs)) if keep_gals[i]]
    
    xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl, segm = grid_smooth(i_ra_f, i_dec_f, 2.0, width, height)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(S, extent=extent, interpolation='nearest',cmap=cm.gray)
    plt.imshow(segm, extent=extent, alpha=0.5)
    cbar_S = plt.colorbar()
    plt.xlabel('RA (arcmin)')
    plt.ylabel('Dec (arcmin)')
    plt.xlim(0,yedges[-1])
    plt.ylim(0,xedges[-1])
    plt.show()
    return len(ras), S[x_cent_S][y_cent_S], image[:-5], ra_corner, dec_corner, 2.0, 0.0

n, dists, title, ra, dec, fwhm, dm = galmap('HI1151+20_g.fits')
distfit(n, dists, title, ra, dec, fwhm, dm)
