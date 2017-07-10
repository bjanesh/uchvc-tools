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
# from pyraf import iraf
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
    print object
    # find the right row
    coord = [this for i,this in enumerate(coords) if object.upper() in name[i]][0]
    
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
    ra_hi = coord_hi.ra.deg
    dec_hi = coord_hi.dec.deg
    return ra_hi, dec_hi
    
def dist2HIcentroid(ra, dec, ra_hi, dec_hi):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    c_hi = SkyCoord(ra = ra_hi, dec = dec_hi, unit=(u.deg, u.deg))
    c_peak = SkyCoord(ra = ra, dec = dec, unit=(u.hourangle, u.deg))
    sep = c_hi.separation(c_peak)
    # print c_hi, c_peak, sep.arcsecond
    return sep.arcsecond

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
        DEC = '{0}{1:02d}:{2:02d}:{3:04.1f}'.format(ds, deg, decM, decS)

    if ra:
        if str(ra)[0] == '-':
            rs, ra = '-', abs(ra)
        raH = int(ra/15)
        raM = int(((ra/15)-raH)*60)
        if round:
            raS = int(((((ra/15)-raH)*60)-raM)*60)
        else:
            raS = ((((ra/15)-raH)*60)-raM)*60
        RA = '{0}{1:02d}:{2:02d}:{3:04.1f}'.format(rs, raH, raM, raS)
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC    
    
def make_filter(dm, filter_file):
    # the color-magnitude filter we're going to use in abs. mag. 
    # based on Girardi ugriz isochrones, see Walsh et al. and Betsey's thesis chapter
    if filter_file.endswith("iso_filter.txt") :
        gi_iso,i_iso = np.loadtxt(filter_file, usecols=(0,1),unpack=True)
        i_m_iso = i_iso + dm
        verts = zip(gi_iso,i_m_iso)        # set up the Path necessary for testing membership
        cm_filter = Path(verts)
        
    elif filter_file == 'none':
        gi_iso = None
        i_m_iso = None
        cm_filter = None
    else :
        g_iso,i_iso = np.loadtxt(filter_file, usecols=(8,10),unpack=True)
        gi_iso = g_iso - i_iso
        i_m_iso = i_iso + dm
        verts = zip(gi_iso,i_m_iso)        # set up the Path necessary for testing membership
        cm_filter = Path(verts)     # scale the filter to the DM: to find the apparent mag. just add the DM
    return cm_filter, gi_iso, i_m_iso

def filter_sources(i_mag, i_ierr, gmi, gmi_err, cm_filter, filter_sig = 1):
    if cm_filter == None:
        stars_f = [True] * len(gmi)
    # first get the stars that are in the filter (points & 1-sig error bars!)
    else:
        stars_f = list(gmi)
        for i in range(len(gmi)) :
            nsteps_color = int(abs((float(filter_sig)*gmi_err[i])//0.001))
            nsteps_mag = int(abs((float(filter_sig)*i_ierr[i])//0.001))
            
            if nsteps_color == 0 :
                nsteps_color = 1
            if nsteps_mag == 0 :
                nsteps_mag = 1
            
            # each error bar, sampled as points
            cm_points_l = [(gmi[i]-0.01*j*gmi_err[i],i_mag[i]) for j in range(nsteps_color)]
            cm_points_r = [(gmi[i]+0.01*j*gmi_err[i],i_mag[i]) for j in range(nsteps_color)]
            cm_points_u = [(gmi[i],i_mag[i]-0.01*j*i_ierr[i]) for j in range(nsteps_mag)]
            cm_points_d = [(gmi[i],i_mag[i]+0.01*j*i_ierr[i]) for j in range(nsteps_mag)]
                    
            # check if the errorbar points fall in the filter        
            stars_f_l = cm_filter.contains_points(cm_points_l)        
            stars_f_r = cm_filter.contains_points(cm_points_r)        
            stars_f_u = cm_filter.contains_points(cm_points_u)        
            stars_f_d = cm_filter.contains_points(cm_points_d)        
            
            # if any part of any error bar is in the filter, it gets selected (True)
            stars_f[i] = any(stars_f_l) | any(stars_f_r) | any(stars_f_u) | any(stars_f_d)  
    
    # figure out how many stars are in the filter
    check = [stars_f[i] for i in range(len(stars_f)) if (stars_f[i])]
    # print len(check), "stars in filter"
    return stars_f

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
    return xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl 

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
    # print 'Significance of detection:','{0:6.3f}%'.format(pct)

    if samples > 10000 and pct > 95.:
        plt.clf()
        plt.figure(figsize=(9,4))
        plt.plot(x, lognorm.pdf(x, al, loc=loc, scale=beta),'r-', lw=2, alpha=0.6, label='lognormal distribution')
        plt.scatter(centers, bins, edgecolors='none', label='histogram of $\sigma$ from '+repr(samples)+' \nuniform random samples')
        ax = plt.subplot(111)
        plt.plot([dists,dists],[-1.0,2.0],'k--', lw=2, alpha=1.0, label='best '+title+' detection') 
        plt.ylim(0,1.1)
        plt.xlim(2,12)
        plt.xlabel('$\sigma$ above local mean')
        plt.ylabel('$P(\sigma = X)$')
        plt.legend(loc='best', frameon=True)
        # ax.set_aspect(3)
        # plt.show()
        plt.savefig('{:s}_{:5.2f}_{:3.1f}_dist.pdf'.format(title,dm,fwhm))
        
    return pct
    
def main(argv):
    home_root = os.environ['HOME']
    # set defaults for command line flags
    imexam_flag = False
    disp_flag = False
    fwhm = 2.0
    fwhm_string = repr(fwhm)
    filter_file = os.path.dirname(os.path.abspath(__file__))+'/filter.txt'
    filter_string = 'old'
    dm2 = None
    
    try:
        opts, args = getopt.getopt(argv,"h",["imexam","disp","fwhm","dm=","dm2=","filter"])
    except getopt.GetoptError:
        print 'magfilter.py --fwhm=<fwhm in arcmin> --dm=<DM in mag> --disp=(yes/no)'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'magfilter.py --fwhm=<fwhm in arcmin> --dm=<DM in mag> --disp=(yes/no)'
            sys.exit()
        elif opt in ("--fwhm"):
            fwhm = float(arg)        # in arcmin (7.5 pixels = 1 arcmin)
            fwhm_string = arg        # this is the smoothing scale, not a stellar profile
        elif opt in ("--dm"):
            dm = float(arg)        # in mag
            dm_string = arg
        elif opt in ("--dm2"):
            dm2 = float(arg)
        elif opt in ("--imexam"):
            imexam_flag = True
        elif opt in ("--disp"):
            disp = arg        # in mag
            if disp == 'yes' :
                disp_flag = True
            else :
                disp_flag = False
        elif opt in ("--filter"):
            filt = arg
            if arg == 'old' :
                filter_file = os.path.dirname(os.path.abspath(__file__))+'/filter.txt'
                filter_string = 'old'
            elif arg == 'none':
                filter_file = 'none'
                filter_string = 'none'
            else: 
                filter_file = os.path.dirname(os.path.abspath(__file__))+'/iso_filter.txt'
                filter_string = 'iso'
    
    fwhm_string = fwhm_string.replace('.','_')
    # print "Getting fits files..."
    # Load the FITS header using astropy.io.fits
    for file_ in os.listdir("./"):
        if file_.endswith("i.fits"):
            fits_file_i = file_
 
    for file_ in os.listdir("./"):
        if file_.endswith("g.fits"):
            fits_file_g = file_
              
    fits_i = fits.open(fits_file_i)
    fits_g = fits.open(fits_file_g)
    # print "Opened fits files:",fits_file_g,"&",fits_file_i
    
    # objid = fits_i[0].header['OBJECT']
    title_string = fits_file_i.split('_')[0]       # get the name part of the filename.

    # set up some filenames
    mag_file = 'calibrated_mags.dat'

    # read in magnitudes, colors, and positions(x,y)
    # gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir,fwhm_sr= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8,11),unpack=True)
    gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
    # print len(gxr), "total stars"
    fwhm_sr = np.zeros_like(gxr)
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
    i_ierr = np.array([i_ierrr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
    gmi = [gmir[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
    fwhm_s = [fwhm_sr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
    gmi_err = np.array([np.sqrt(g_ierrr[i]**2 + i_ierrr[i]**2) for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
    
    i_ierrAVG, bedges, binid = ss.binned_statistic(i_mag,i_ierr,statistic='median',bins=10,range=[15,25])
    gmi_errAVG, bedges, binid = ss.binned_statistic(i_mag,gmi_err,statistic='median',bins=10,range=[15,25])
    
    bcenters = (bedges[:-1] + bedges[1:]) / 2
    bxvals = [3.75,3.75,3.75,3.75,3.75,3.75,3.75,3.75,3.75,3.75]
    # print bcenters
    # print i_ierrAVG
    # print gmi_errAVG
    # print len(gx), "after color+mag error cut"
    # nid = np.loadtxt(mag_file,usecols=(0,),dtype=int,unpack=True)
    pixcrd = zip(ix,iy)
    
    
    # print "Reading WCS info from image header..."
    # Parse the WCS keywords in the primary HDU
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    w = wcs.WCS(fits_i[0].header)
    # print fits_i[0].header['naxis1'], fits_i[0].header['naxis2']
    footprint = w.calc_footprint()
    se_corner = footprint[0]
    ne_corner = footprint[1]
    nw_corner = footprint[2]
    sw_corner = footprint[3]
    # print se_corner, ne_corner, nw_corner, sw_corner
    width = (ne_corner[0]-nw_corner[0])*60.
    height = (ne_corner[1]-se_corner[1])*60.
    # print width, height
    
    # Print out the "name" of the WCS, as defined in the FITS header
    # print w.wcs.name
    
    # Print out all of the settings that were parsed from the header
    # w.wcs.print_contents()
    
    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.all_pix2world(pixcrd, 1)
    ra_corner, dec_corner = w.all_pix2world(0,0,1)
    ra_c_d,dec_c_d = deg2HMS(ra=ra_corner, dec=dec_corner, round=True)
    # print 'Corner RA:',ra_c_d,':: Corner Dec:',dec_c_d
    
    fwhm_i = 12.0 #fits_i[0].header['FWHMPSF']
    fwhm_g = 9.0 # fits_g[0].header['FWHMPSF']
    
    # print 'Image FWHM :: g = {0:5.3f} : i = {1:5.3f}'.format(fwhm_g,fwhm_i)
    
    fits_i.close()
    fits_g.close()
    
    # split the ra and dec out into individual arrays and transform to arcmin from the corner
    i_ra = [abs((world[i,0]-ra_corner)*60) for i in range(len(world[:,0]))]
    i_dec = [abs((world[i,1]-dec_corner)*60) for i in range(len(world[:,1]))]
    # also preserve the decimal degrees for reference
    i_rad = [world[i,0] for i in range(len(world[:,0]))]
    i_decd = [world[i,1] for i in range(len(world[:,1]))]
    
    
    i_magBright = [i_mag[i] for i in range(len(i_mag)) if (i_mag[i] < 22.75)]
    g_magBright = [g_mag[i] for i in range(len(i_mag)) if (i_mag[i] < 22.75)]
    ixBright = [ix[i] for i in range(len(i_mag)) if (i_mag[i] < 22.75)]
    iyBright = [iy[i] for i in range(len(i_mag)) if (i_mag[i] < 22.75)]
    i_radBright = [i_rad[i] for i in range(len(i_mag)) if (i_mag[i] < 22.75)]
    i_decdBright = [i_decd[i] for i in range(len(i_mag)) if (i_mag[i] < 22.75)]
    
    if not os.path.isfile('brightStars2275.reg'):
        f1 = open('brightStars2275.reg', 'w+')
        for i in range(len(i_magBright)) :
            print >> f1, '{0:12.4f} {1:12.4f} {2:10.5f} {3:9.5f} {4:8.2f} {5:8.2f} {6:8.2f}'.format(ixBright[i],iyBright[i], i_radBright[i], i_decdBright[i], g_magBright[i], i_magBright[i], g_magBright[i]-i_magBright[i])
        f1.close()
    
    i_magRed = [i_mag[i] for i in range(len(i_mag)) if (gmi[i] > 1.75)]
    g_magRed = [g_mag[i] for i in range(len(i_mag)) if (gmi[i] > 1.75)]
    ixRed = [ix[i] for i in range(len(i_mag)) if (gmi[i] > 1.75)]
    iyRed = [iy[i] for i in range(len(i_mag)) if (gmi[i] > 1.75)]
    i_radRed = [i_rad[i] for i in range(len(i_mag)) if (gmi[i] > 1.75)]
    i_decdRed = [i_decd[i] for i in range(len(i_mag)) if (gmi[i] > 1.75)]
    
    if not os.path.isfile('redStars175.reg'):
        f1 = open('redStars175.reg', 'w+')
        for i in range(len(i_magRed)) :
            print >> f1, '{0:12.4f} {1:12.4f} {2:10.5f} {3:9.5f} {4:8.2f} {5:8.2f} {6:8.2f}'.format(ixRed[i],iyRed[i], i_radRed[i], i_decdRed[i], g_magRed[i], i_magRed[i], g_magRed[i]-i_magRed[i])
        f1.close()
    
    if dm2 and filter_string != 'none':
        dms = np.arange(dm,dm2,0.02)
    else:
        dms = [dm]

    search = open('search.txt','w+')
    
    for dm in dms:
        mpc = pow(10,((dm + 5.)/5.))/1000000.
        dm_string = '{:5.2f}'.format(dm).replace('.','_')
        
        out_file = filter_string + '_' + fwhm_string + '_' + dm_string + '_' + title_string + '.pdf'
        mark_file = 'f_list_' + filter_string + '_' + fwhm_string + '_' + dm_string + '_' + title_string + '.reg'
        filter_reg = 'f_reg_' + filter_string + '_' + fwhm_string + '_' + dm_string + '_' + title_string + '.reg'
        circ_file = 'c_list_' + filter_string + '_' + fwhm_string + '_' + dm_string + '_' + title_string + '.reg'
        fcirc_file = 'fc_list_' + filter_string + '_' + fwhm_string + '_' + dm_string + '_' + title_string + '.reg'
        ds9_file = 'circles_' + filter_string + '_' + fwhm_string + '_' + dm_string + '_' + title_string + '.reg'
        circles_file = 'region_coords.dat'
        
        cm_filter, gi_iso, i_m_iso = make_filter(dm, filter_file)
        stars_f = filter_sources(i_mag, i_ierr, gmi, gmi_err, cm_filter, filter_sig = 1)
        
        xy_points = zip(i_ra,i_dec)
        
        # make new vectors containing only the filtered points
        
        i_mag_f = [i_mag[i] for i in range(len(i_mag)) if (stars_f[i])]
        g_mag_f = [g_mag[i] for i in range(len(i_mag)) if (stars_f[i])]
        gmi_f = [gmi[i] for i in range(len(i_mag)) if (stars_f[i])]
        i_ra_f = [i_ra[i] for i in range(len(i_mag)) if (stars_f[i])]
        i_dec_f = [i_dec[i] for i in range(len(i_mag)) if (stars_f[i])]
        i_rad_f = [i_rad[i] for i in range(len(i_mag)) if (stars_f[i])]
        i_decd_f = [i_decd[i] for i in range(len(i_mag)) if (stars_f[i])]
        i_x_f = [ix[i] for i in range(len(i_mag)) if (stars_f[i])]
        i_y_f = [iy[i] for i in range(len(i_mag)) if (stars_f[i])]
        fwhm_sf = [fwhm_s[i] for i in range(len(i_mag)) if (stars_f[i])]
        n_in_filter = len(i_mag_f)
        
        xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl = grid_smooth(i_ra_f, i_dec_f, fwhm, width, height)
        pct = distfit(n_in_filter,S[x_cent_S][y_cent_S],title_string,width,height,fwhm,dm)
        if pct > 95 :
            pct = distfit(n_in_filter,S[x_cent_S][y_cent_S],title_string,width,height,fwhm,dm, samples=25000)
        
        # make a circle to highlight a certain region
        cosd = lambda x : np.cos(np.deg2rad(x))
        sind = lambda x : np.sin(np.deg2rad(x))
        x_circ = [yedges[y_cent] + 3.0*cosd(t) for t in range(0,359,1)]
        y_circ = [xedges[x_cent] + 3.0*sind(t) for t in range(0,359,1)]
        
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
        fwhm_sc = [fwhm_s[i] for i in range(len(i_mag)) if (stars_circ[i])]
        
        # make a random reference cmd to compare to
        if not os.path.isfile('refCircle.center'):
            rCentx = 16.0*np.random.random()+2.0
            rCenty = 16.0*np.random.random()+2.0
            with open('refCircle.center','w+') as rc:
                print >> rc, '{:8.4f} {:8.4f}'.format(rCentx, rCenty)
        else :
            rCentx, rCenty = np.loadtxt('refCircle.center', usecols=(0,1), unpack=True)
                
        x_circr = [rCentx + 3.0*cosd(t) for t in range(0,359,1)]
        y_circr = [rCenty + 3.0*sind(t) for t in range(0,359,1)]
        
        verts_circr = zip(x_circr,y_circr)
        rcirc_filter = Path(verts_circr)
        
        stars_circr = rcirc_filter.contains_points(xy_points)    
        
        i_mag_cr = [i_mag[i] for i in range(len(i_mag)) if (stars_circr[i])]
        gmi_cr = [gmi[i] for i in range(len(i_mag)) if (stars_circr[i])]
        i_ra_cr = [i_ra[i] for i in range(len(i_mag)) if (stars_circr[i])]
        i_dec_cr = [i_dec[i] for i in range(len(i_mag)) if (stars_circr[i])]
        i_rad_cr = [i_rad[i] for i in range(len(i_mag)) if (stars_circr[i])]
        i_decd_cr = [i_decd[i] for i in range(len(i_mag)) if (stars_circr[i])]
        i_x_cr = [ix[i] for i in range(len(i_mag)) if (stars_circr[i])]
        i_y_cr = [iy[i] for i in range(len(i_mag)) if (stars_circr[i])]
        fwhm_scr = [fwhm_s[i] for i in range(len(i_mag)) if (stars_circr[i])]
        
        i_mag_fc = [i_mag[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        i_ierr_fc = [i_ierr[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        g_ierr_fc = [g_ierr[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        g_mag_fc = [i_mag[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        gmi_fc = [gmi[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        i_ra_fc = [i_ra[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        i_dec_fc = [i_dec[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        i_rad_fc = [i_rad[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        i_decd_fc = [i_decd[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        i_x_fc = [ix[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        i_y_fc = [iy[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        fwhm_sfc = [fwhm_s[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        
        # print len(i_mag_fc), 'filter stars in circle'
        
        # print 'max i mag in circle = ', min(i_mag_fc)
        
        rs = np.array([45, 55, 65, 75, 85, 90,180])
        for r in rs:    
            x_circ = [yedges[y_cent] + r/60.*cosd(t) for t in range(0,359,1)]
            y_circ = [xedges[x_cent] + r/60.*sind(t) for t in range(0,359,1)]
    
            verts_circ = zip(x_circ,y_circ)
            circ_filter = Path(verts_circ)
    
            stars_circ = circ_filter.contains_points(xy_points)
            i_x_fc = [ix[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
            i_y_fc = [iy[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        
            fcirc_file = 'circle'+repr(r)+'.txt'
            with open(fcirc_file,'w+') as f3:
                for i,x in enumerate(i_x_fc):
                    print >> f3, i_x_fc[i], i_y_fc[i]
        
        i_mag_fcr = [i_mag[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        i_ierr_fcr = [i_ierr[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        g_ierr_fcr = [g_ierr[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        g_mag_fcr = [i_mag[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        gmi_fcr = [gmi[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        i_ra_fcr = [i_ra[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        i_dec_fcr = [i_dec[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        i_rad_fcr = [i_rad[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        i_decd_fcr = [i_decd[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        i_x_fcr = [ix[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        i_y_fcr = [iy[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        fwhm_sfcr = [fwhm_s[i] for i in range(len(i_mag)) if (stars_circr[i] and stars_f[i])]
        # print len(i_mag_fcr), 'filter stars in ref. circle'
        
        # 
        # print "{0:3d} stars in filter, {1:3d} stars in circle, {2:3d} stars in both.".format(len(i_mag_f),len(i_mag_c),len(i_mag_fc))
        # for i in range(len(i_x_fc)) :
        #     print (i_x_fc[i],i_y_fc[i])
        
        rcirc_c_x = ra_corner-(rCentx/60.)
        rcirc_c_y = (rCenty/60.)+dec_corner
        rcirc_pix_x, rcirc_pix_y = w.wcs_world2pix(rcirc_c_x,rcirc_c_y,1)
        ra_cr, dec_cr = w.all_pix2world(rcirc_pix_x, rcirc_pix_y,1)
        ra_cr_d,dec_cr_d = deg2HMS(ra=ra_cr, dec=dec_cr, round=False)
        
        circ_c_x = ra_corner-(yedges[y_cent]/60.)
        circ_c_y = (xedges[x_cent]/60.)+dec_corner
        circ_pix_x, circ_pix_y = w.wcs_world2pix(circ_c_x,circ_c_y,1)
        ra_c, dec_c = w.all_pix2world(circ_pix_x, circ_pix_y,1)
        ra_c_d,dec_c_d = deg2HMS(ra=ra_c, dec=dec_c, round=False)
        # print 'Peak RA:',ra_c_d,':: Peak Dec:',dec_c_d
        
        hi_c_ra, hi_c_dec = getHIcentroid(title_string)
        hi_c_x, hi_c_y = abs((hi_c_ra-ra_corner)*60), abs((hi_c_dec-dec_corner)*60)
        hi_x_circ = [hi_c_x + pltsig*cosd(t) for t in range(0,359,1)]
        hi_y_circ = [hi_c_y + pltsig*sind(t) for t in range(0,359,1)]
        
        hi_pix_x,hi_pix_y = w.wcs_world2pix(hi_c_ra,hi_c_dec,1)
        
        sep = dist2HIcentroid(ra_c_d, dec_c_d, hi_c_ra, hi_c_dec)
        # print hi_pix_x, hi_pix_y
        # print ra_cr, dec_cr
        # print ra_c, dec_c
        # print hi_c_ra, hi_c_dec
        
        print 'm-M = {:5.2f} | d = {:4.2f} Mpc | α = {:s}, δ = {:s}, Δ = {:5.1f}" | N = {:4d} | σ = {:6.3f} | β = {:6.3f}%'.format(dm, mpc, ra_c_d, dec_c_d, sep, n_in_filter, S[x_cent_S][y_cent_S], pct)
        print >> search, '{:5.2f} {:4.2f} {:s} {:s} {:5.1f} {:4d} {:6.3f} {:6.3f}'.format(dm, mpc, ra_c_d, dec_c_d, sep, n_in_filter, S[x_cent_S][y_cent_S], pct)        

        #iraf.imutil.hedit(images=fits_g, fields='PV*', delete='yes', verify='no')
        #iraf.imutil.hedit(images=fits_i, fields='PV*', delete='yes', verify='no') 
        if pct > 95.:
            with open(ds9_file,'w+') as ds9:
                print >> ds9, "fk5;circle({:s},{:s},2') # color=yellow width=2 label=ref".format(ra_cr, dec_cr)
                print >> ds9, "fk5;circle({:s},{:s},2') # color=magenta width=2 label=detection".format(ra_c, dec_c)
                print >> ds9, "fk5;circle({:f},{:f},0.5') # color=black width=2 label=HI".format(hi_c_ra, hi_c_dec)
            
            f1 = open(filter_reg, 'w+')
            for i in range(len(i_x_f)) :
                print >> f1, 'image;circle({0:8.2f},{1:8.2f},10) # color=green width=2 edit=0 move=0 delete=0 highlite=0'.format(i_x_f[i],i_y_f[i])
            f1.close()
                
            f1 = open(mark_file, 'w+')
            for i in range(len(i_x_f)) :
                print >> f1, '{0:8.2f} {1:8.2f} {2:12.8f} {3:12.8f} {4:8.2f} {5:8.2f} {6:8.2f} {7:7.3f}'.format(i_x_f[i],i_y_f[i],i_rad_f[i],i_decd_f[i],i_mag_f[i],g_mag_f[i],gmi_f[i],fwhm_sf[i])
            f1.close()
            
            f2 = open(circ_file, 'w+')
            for i in range(len(i_x_c)) :
                print >> f2, '{0:8.2f} {1:8.2f} {2:12.8f} {3:12.8f} {4:8.2f} {5:8.2f} {6:8.2f} {7:7.3f}'.format(i_x_c[i],i_y_c[i],i_rad_c[i],i_decd_c[i],i_mag_c[i],i_mag_c[i]+gmi_c[i],gmi_c[i],fwhm_sc[i])
            f2.close()
            
            with open(fcirc_file,'w+') as f3:
                for i,x in enumerate(i_x_fc):
                    print >> f3, i_x_fc[i], i_y_fc[i], fwhm_sfc[i]
                    
            with open(circles_file,'w+') as f4:
                print >> f4, circ_pix_x, circ_pix_y
            
            # center_file = 'im_cens_'+dm_string+'_'+fwhm_string+'.dat'
            # im_cens = open(center_file,'w+')
            # for k,xp in enumerate(tbl):
            #     circ_c_x = (yedges[tbl['maxval_xpos'][k]]/60.)+ra_corner
            #     circ_c_y = (xedges[tbl['maxval_ypos'][k]]/60.)+dec_corner
            #     circ_pix_x, circ_pix_y = w.wcs_world2pix(circ_c_x,circ_c_y,1)
            #     # ra_c, dec_c = w.all_pix2world(circ_pix_x, circ_pix_y,1)
            #     # print circ_c_x, circ_c_y
            #     print >> im_cens, -circ_pix_x, circ_pix_y, 'circle_center'+repr(k)
            # print >> im_cens, hi_pix_x, hi_pix_y, 'HI_centroid'
            # im_cens.close()
            
            mark_radius = int(pltsig*60/0.11)
            ann_inn = int(np.ceil(pltsig*60/0.11/100.0))*100
            ann_out = int(np.ceil(pltsig*60/0.11/100.0))*100+100
            mark_radii = str(repr(mark_radius-2)+','+repr(mark_radius-1)+','+repr(mark_radius)+','+repr(mark_radius+1)+','+repr(mark_radius+2))
            anni_radii = str(repr(ann_inn-2)+','+repr(ann_inn-1)+','+repr(ann_inn)+','+repr(ann_inn+1)+','+repr(ann_inn+2))
            anno_radii = str(repr(ann_out-2)+','+repr(ann_out-1)+','+repr(ann_out)+','+repr(ann_out+1)+','+repr(ann_out+2))
            # print mark_radius, ann_inn, ann_out
            
            if disp_flag :
                from pyraf import iraf
                iraf.tv(_doprint=0)
                iraf.unlearn(iraf.tv.tvmark)
                iraf.tv.tvmark.setParam('label',"no")
                iraf.tv.tvmark.setParam('pointsize',7)
                iraf.tv.tvmark.setParam('mark',"circle")
                # iraf.tv.display(image=fits_file_g, frame=1)
                iraf.tv.display(image=fits_file_i, frame=1)
                iraf.tv.tvmark(frame=1, coords=mark_file, radii="14,15,16", color=207)
                iraf.tv.tvmark(frame=1, coords=circ_file, radii="20,21,22", color=209)
                iraf.tv.tvmark(frame=1, coords='brightStars22.reg', radii="26,27,28", color=205)
                iraf.tv.tvmark(frame=1, coords='redStars175.reg', radii="32,33,34", color=204)
                iraf.tv.tvmark(frame=1, coords=center_file, txsize=4, mark="plus", color=208, label="yes")
                iraf.tv.tvmark(frame=1, coords=center_file, radii=mark_radii, color=208)
                iraf.tv.tvmark(frame=1, coords=center_file, radii=anni_radii, color=208)
                iraf.tv.tvmark(frame=1, coords=center_file, radii=anno_radii, color=208)
                # iraf.tv.tvmark(frame=2, coords=mark_file, radii="14,15,16", color=207, label="yes")
                # iraf.tv.tvmark(frame=2, coords=circ_file, radii="20,21,22", color=209)
                # iraf.tv.tvmark(frame=2, coords=center_file, txsize=4, mark="plus", color=208, label="yes")
                # iraf.tv.tvmark(frame=2, coords=center_file, radii=mark_radii, color=208)
        
            # setup the pdf output
            pp = PdfPages('f_'+ out_file)
            plt.clf()
            plt.figure(figsize=(11,8.5))
            # plot
            # print "Plotting for m-M = ",dm
            ax0 = plt.subplot(2,2,1)
            plt.scatter(i_ra, i_dec,  color='black', marker='o', s=1, edgecolors='none')
            plt.plot(x_circ,y_circ,linestyle='-', color='magenta')
            plt.plot(x_circr,y_circr,linestyle='-', color='gold')
            plt.plot(hi_x_circ,hi_y_circ,linestyle='-', color='black')
            # plt.scatter(i_ra_c, i_dec_c,  color='red', marker='o', s=3, edgecolors='none')
            plt.scatter(i_ra_f, i_dec_f,  c=gmi_f, marker='o', s=(30-np.array(i_mag_f)), edgecolors='none', cmap=cm.rainbow)
            plt.clim(0,2)
            plt.colorbar()
            plt.ylabel('Dec (arcmin)')
            plt.xlim(0,max(i_ra))
            plt.ylim(0,max(i_dec))
            plt.title('sky positions')
            ax0.set_aspect('equal')
        
            ax1 = plt.subplot(2,2,2)
            
            if os.path.isfile('i_gmi_compl.gr.out'):
                gmiCompl, iCompl = np.loadtxt('i_gmi_compl.gr.out',usecols=(0,1),unpack=True)
                plt.plot(gmiCompl,iCompl, linestyle='--', color='green')
            if os.path.isfile('i_gmi_compl.gr2.out'):
                gmiCompl, iCompl = np.loadtxt('i_gmi_compl.gr2.out',usecols=(0,1),unpack=True)
                plt.plot(gmiCompl,iCompl, linestyle='--', color='red')
            if os.path.isfile('i_gmi_compl2.out'):
                iCompl,gmiCompl = np.loadtxt('i_gmi_compl2.out',usecols=(0,1),unpack=True)
                plt.plot(gmiCompl,iCompl, linestyle='--', color='blue')
            
            plt.plot(gi_iso,i_m_iso,linestyle='-', color='blue')
            plt.scatter(gmi, i_mag,  color='black', marker='o', s=1, edgecolors='none')
            plt.scatter(gmi_f, i_mag_f,  color='red', marker='o', s=3, edgecolors='none')
            # plt.scatter(gmi_c, i_mag_c,  color='red', marker='o', s=3, edgecolors='none')
            plt.errorbar(bxvals, bcenters, xerr=i_ierrAVG, yerr=gmi_errAVG, linestyle='None', color='black', capsize=0, ms=0)
            plt.tick_params(axis='y',left='on',right='off',labelleft='on',labelright='off')
            ax1.yaxis.set_label_position('left')
            plt.ylabel('$i_0$')
            plt.xlabel('$(g-i)_0$')
            plt.ylim(25,15)
            plt.xlim(-1,4)
            plt.title('m-M = ' + '{:5.2f}'.format(dm) + ' (' + '{0:4.2f}'.format(mpc) +  ' Mpc)')
            ax1.set_aspect(0.5)
        
            ax2 = plt.subplot(2,2,3)
        
            extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
            plt.imshow(S, extent=extent, interpolation='nearest',cmap=cm.gray)
            # plt.imshow(segm, extent=extent, cmap=rand_cmap, alpha=0.5)
            cbar_S = plt.colorbar()
            # cbar_S.tick_params(labelsize=10)
            plt.plot(x_circ,y_circ,linestyle='-', color='magenta')
            plt.plot(x_circr,y_circr,linestyle='-', color='gold')
            plt.plot(hi_x_circ,hi_y_circ,linestyle='-', color='black')
            # X, Y = np.meshgrid(xedges,yedges)
            # ax3.pcolormesh(X,Y,grid_gaus)
            plt.xlabel('RA (arcmin)')
            plt.ylabel('Dec (arcmin)')
            plt.title('smoothed star map')
            # plt.ylabel('Dec (arcmin)')
            plt.xlim(0,max(i_ra))
            plt.ylim(0,max(i_dec))
            ax2.set_aspect('equal')
        
            # ax3 = plt.subplot(2,2,4)
            ax3 = plt.subplot2grid((2,4), (1,2))
            plt.scatter(gmi_c, i_mag_c,  color='black', marker='o', s=3, edgecolors='none')
            plt.scatter(gmi_fc, i_mag_fc,  color='red', marker='o', s=3, edgecolors='none')    
            plt.tick_params(axis='y',left='on',right='on',labelleft='off',labelright='off')
            ax0.yaxis.set_label_position('left')
            plt.title('in detection circle')
            plt.xlabel('$(g-i)_0$')
            plt.ylabel('$i_0$')
            plt.ylim(25,15)
            plt.xlim(-1,4)
            # ax3.set_aspect(0.5)    
        
            ax4 = plt.subplot2grid((2,4), (1,3), sharey=ax3)
            plt.scatter(gmi_cr, i_mag_cr,  color='black', marker='o', s=3, edgecolors='none')
            plt.scatter(gmi_fcr, i_mag_fcr,  color='red', marker='o', s=3, edgecolors='none')    
            plt.tick_params(axis='y',left='on',right='on',labelleft='off',labelright='on')
            plt.title('in ref. circle')
            ax0.yaxis.set_label_position('left')
            plt.xlabel('$(g-i)_0$')
            plt.ylim(25,15)
            plt.xlim(-1,4)
            # ax3.set _aspect(0.5)
            plt.tight_layout()
            # plt.suptitle(title_string+ ' -- ' + fwhm_string + ' arcmin smoothing')
            plt.show
        
            pp.savefig()
            pp.close()
        # else :
            # print 'detection at dm='+repr(dm)+' not significant enough to look at plots'
    search.close()    
    if imexam_flag :
        from pyraf import iraf
        iraf.unlearn(iraf.tv.imexamine, iraf.rimexam)
        iraf.tv.rimexam.setParam('radius',int(fwhm_i))
        iraf.tv.rimexam.setParam('rplot',12.)
        iraf.tv.rimexam.setParam('fittype','gaussian')
        iraf.tv.rimexam.setParam('center','yes')
        iraf.tv.imexamine(input=fits_file_i, frame=2)

if __name__ == "__main__":
    main(sys.argv[1:])    
