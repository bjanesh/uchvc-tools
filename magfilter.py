#! /usr/local/bin/python3
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
from scipy import signal
from odi_calibrate import query, filtercomment, usage, write_header
from photutils import detect_sources, source_properties
from photutils.utils import random_cmap
try :
    from scipy import ndimage
except ImportError :
    print('bad import')

def downloadSDSSgal(img1, img2):
    formats = ['csv','xml','html']
    
    astro_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
    public_url='http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
    
    default_url=public_url
    default_fmt='csv'
    
    try: 
        import sys
        import numpy as np
        from astropy.io import fits
        from astropy import wcs
        import os
        import string
    except ImportError:
        print("You should 'pip install astropy' before you try to run this program") 
    
    print('fetching SDSS data from \n--> '+public_url)
    
    image = img1
    
    # read in the image header and save it to a variable for non-destructive editing
    hdulist = fits.open(image)
    hdr = hdulist[0].header
    hdulist.close()
    # get the image dimensions
    xdim = hdr['NAXIS1']
    ydim = hdr['NAXIS2']
    
    # and find the image center
    xc = xdim/2.0
    yc = ydim/2.0
    
    # get the CD matrix keywords
    cd11 = hdr['CD1_1']
    cd22 = hdr['CD2_2']
    # try to load cd12 and cd21, if they don't exist, set them to zero
    try :
        cd12 = hdr['CD1_2']
    except:
        cd12 = 0.0
    try :
        cd21 = hdr['CD2_1']
    except:
        cd21 = 0.0
    
    # get rid of keywords starting with PV, they don't work with astropy.wcs
    # and steven thinks they are redundant anyway
    pvlist = hdr['PV*']
    for pv in pvlist:
        hdr.remove(pv)
        
    # open the second fits image
    hdulist = fits.open(img2)
    hdr_r = hdulist[0].header
    hdulist.close()
    
    pvlist = hdr_r['PV*']
    for pv in pvlist:
        hdr_r.remove(pv)
    
    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdr)
    w_r = wcs.WCS(hdr_r)
    
    # Some pixel coordinates of interest (these are the image centers)
    pixcrd = np.array([[xc,yc]], np.float_)
    
    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 1-based (Fortran-like) coordinates.
    world = w.wcs_pix2world(pixcrd, 1)
    # print(world)    
    rac = world[0][0]
    decc = world[0][1]
    
    # get the biggest radius of the image in arcminutes
    pixscal1 = 3600*abs(cd11)
    pixscal2 = 3600*abs(cd22)
    xas = pixscal1 * xdim # in arcseconds
    yas = pixscal2 * ydim
    xam = xas/60    # to arcminutes
    yam = yas/60
    #print(xam,yam)
    #radius for query: sqrt2 = 1.414
    sizeam = 1.414*(xam+yam)/4
    # print sizeam
    
    if not os.path.isfile(image[:-5]+'.sdssgal'):
        # build the SDSS query
        qry = "select O.ra, O.dec, O.u, O.err_u, O.g, \nO.err_g, O.r, O.err_r, O.i, \nO.err_i, O.z, O.err_z, O.probPSF \nfrom \ndbo.fGetNearbyObjEq("+repr(rac)+","+repr(decc)+","+repr(sizeam)+") \nas N inner join PhotoObjAll as O on O.objID = N.objID order by N.distance"
    
        # print it to the terminal
        print('with query\n-->', qry)
        url = default_url
        fmt = default_fmt
        writefirst = 1
        verbose = 0
    
        # actually do the query
    
        ofp = open(image[:-5]+'.sdssgal','w+')
        if verbose:
            write_header(ofp,'#',url,qry)
        file_ = query(qry,url,fmt)
        # Output line by line (in case it's big)
        line = file_.readline()
        if line.startswith("ERROR"): # SQL Statement Error -> stderr
            ofp = sys.stderr
        if writefirst:
            ofp.write(string.rstrip(line)+os.linesep)
        line = file_.readline()
        while line:
            ofp.write(string.rstrip(line)+os.linesep)
            line = file_.readline()
        ofp.close()
    
    # read in the results
    ras,decs,u,err_u,g,err_g,r,err_r,i,err_i,z,err_z = np.loadtxt(image[:-5]+'.sdssgal',usecols=(0,1,2,3,4,5,6,7,8,9,10,11), unpack=True, delimiter=',', skiprows=2)
    probPSF = np.loadtxt(image[:-5]+'.sdssgal', usecols=(12,), dtype=int, unpack=True, delimiter=',', skiprows=2)
    
    coords2 = list(zip(ras,decs))
    pixcrd2 = w.wcs_world2pix(coords2, 1)
    pixcrd2_r = w_r.wcs_world2pix(coords2, 1)
    
    # keep things that are actually stars (defined as being psf's) and with the right magnitude range (arbitrary)
    
    keep_stars = ((probPSF == 0))# (psfMag_g < gmaglim) & (psfMagErr_g <0.1) & (psfMag_g > gmagbrlim))
    print('keeping', len(np.where(keep_stars)[0]), 'galaxies of', len(g), 'sources')
    
    # then write out separate files for g and i
    with open(image[:-5]+'.galxy','w+') as f1:
        print("# x_g y_g ra dec u uerr g gerr r rerr i ierr z zerr (all modelmags)", file=f1)
        for j,id in enumerate(np.where(keep_stars)[0]):
            print(pixcrd2[id][0], pixcrd2[id][1], ras[id], decs[id], u[id], err_u[id], g[id], err_g[id], r[id], err_r[id], i[id], err_i[id], z[id], err_z[id], file=f1)
            
    with open(img2[:-5]+'.galxy','w+') as f1:
        print("# x_r y_r ra dec u uerr g gerr r rerr i ierr z zerr (all modelmags)", file=f1)
        for j,id in enumerate(np.where(keep_stars)[0]):
            print(pixcrd2_r[id][0], pixcrd2_r[id][1], ras[id], decs[id], u[id], err_u[id], g[id], err_g[id], r[id], err_r[id], i[id], err_i[id], z[id], err_z[id], file=f1)

    
def galaxyMap(fits_file_i, fwhm, dm, filter_file):
    title_string = fits_file_i.split('_')[0]
    x_r, y_r, ra, dec, u, uerr, g, gerr, r, rerr, i, ierr, z, zerr = np.loadtxt(fits_file_i[:-5]+'.galxy', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13), unpack=True)
    gmi = g-i
    gmierr = [np.sqrt(gerr[j]**2 + ierr[j]**2) for j in range(len(g))]
    
    pixcrd = list(zip(x_r,y_r))
    
    fits_i = fits.open(fits_file_i)
    # Parse the WCS keywords in the primary HDU
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    w = wcs.WCS(fits_i[0].header)

    footprint = w.calc_footprint()
    se_corner = footprint[0]
    ne_corner = footprint[1]
    nw_corner = footprint[2]
    sw_corner = footprint[3]
    # print se_corner, ne_corner, nw_corner, sw_corner
    width = (ne_corner[0]-nw_corner[0])*60.
    height = (ne_corner[1]-se_corner[1])*60.

    world = w.all_pix2world(pixcrd, 1)
    ra_corner, dec_corner = w.all_pix2world(0,0,1)
    ra_c_d,dec_c_d = deg2HMS(ra=ra_corner, dec=dec_corner, round=True)

    fits_i.close()
    
    # split the ra and dec out into individual arrays and transform to arcmin from the corner
    i_ra = [abs((world[j,0]-ra_corner)*60) for j in range(len(world[:,0]))]
    i_dec = [abs((world[j,1]-dec_corner)*60) for j in range(len(world[:,1]))]
    
    if dm > 0:
        cm_filter, gi_iso, i_m_iso = make_filter(dm, filter_file)
        stars_f = filter_sources(i, ierr, gmi, gmierr, cm_filter, filter_sig = 1, gal=True)
    else:
        stars_f = [True] * i
    xy_points = list(zip(i_ra,i_dec))
    
    # make new vectors containing only the filtered points
    
    i_f = [i[j] for j in range(len(x_r)) if (stars_f[j])]
    g_f = [g[j] for j in range(len(x_r)) if (stars_f[j])]
    gmi_f = [gmi[j] for j in range(len(x_r)) if (stars_f[j])]
    ra_f = [i_ra[j] for j in range(len(x_r)) if (stars_f[j])]
    dec_f = [i_dec[j] for j in range(len(x_r)) if (stars_f[j])]
    x_f = [x_r[j] for j in range(len(x_r)) if (stars_f[j])]
    y_f = [y_r[j] for j in range(len(x_r)) if (stars_f[j])]
    
    xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl = grid_smooth(ra_f, dec_f, fwhm, width, height)
    
    plt.figure(figsize=(10,4))
    ax1 = plt.subplot(1,2,1)
    ax1.scatter(gmi, i,  color='black', marker='o', s=1, edgecolors='none')
    plt.scatter(gmi_f, i_f,  color='red', marker='o', s=15, edgecolors='none')
    plt.plot(gi_iso,i_m_iso,linestyle='-', color='blue')
    ax1.tick_params(axis='y',left='on',right='off',labelleft='on',labelright='off')
    ax1.yaxis.set_label_position('left')
    ax1.set_ylabel('$i_0$')
    ax1.set_xlabel('$(g-i)_0$')
    ax1.set_ylim(25,15)
    ax1.set_xlim(-1,4)
    
    ax2 = plt.subplot(1,2,2)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    gr = ax2.imshow(S, extent=extent, interpolation='nearest',cmap=cm.gray)
    # ax2.imshow(segm, extent=extent, cmap=rand_cmap, alpha=0.5)
    cbar_S = plt.colorbar(gr)
    cbar_S.set_label('$\sigma$ from local mean')
    ax2.set_xlabel('RA (arcmin)')
    ax2.set_ylabel('Dec (arcmin)')
    ax2.set_xlim(0,width)#max(i_ra))
    ax2.set_ylim(0,height)#max(i_dec))
    plt.tight_layout()
    plt.savefig('galmap_{:s}_{:5.2f}_{:3.1f}.pdf'.format(title_string, dm, fwhm))
    
    return xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl
        
def getHIellipse(object, ra_corner, dec_corner, centroid=False):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    uchvcdb = os.path.dirname(os.path.abspath(__file__))+'/predblist.sort.csv'
    name, coords, ar, br, par = np.loadtxt(uchvcdb, usecols=(1,2,14,15,16), dtype=str, delimiter=',', unpack=True)
    # print object
    # find the right row
    coord = [this for i,this in enumerate(coords) if object.upper() in name[i]][0]
    a = [this for i,this in enumerate(ar) if object.upper() in name[i]][0]
    b = [this for i,this in enumerate(br) if object.upper() in name[i]][0]
    pa = [this for i,this in enumerate(par) if object.upper() in name[i]][0]

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
    
    cosd = lambda x : np.cos(np.deg2rad(x))
    sind = lambda x : np.sin(np.deg2rad(x))
    
    hi_c_x, hi_c_y = abs((ra_hi-ra_corner)*60), abs((dec_hi-dec_corner)*60)
    
    a, b, pa = a.astype(float)/2., b.astype(float)/2., -pa.astype(float)
    
    t = np.array(list(range(0,359,1)))
    ell = np.array([a*cosd(t) , b*sind(t)])
    rot = np.array([[cosd(pa) , -sind(pa)],[sind(pa) , cosd(pa)]])
    ell_rot = np.zeros((2,ell.shape[1]))
    for i in range(ell.shape[1]):
        ell_rot[:,i] = np.dot(rot,ell[:,i])
    hi_x_circ, hi_y_circ = hi_c_x+ell_rot[0,:], hi_c_y+ell_rot[1,:]
    # hi_x_circ = [hi_c_x + a*cosd(t) for t in range(0,359,1)]
    # hi_y_circ = [hi_c_y + b*sind(t) for t in range(0,359,1)]
    if centroid:
        return ra_hi, dec_hi
    else:
        return hi_x_circ, hi_y_circ

def rotateImage(img, angle, pivot):
    padX = [img.shape[1] - pivot[0], pivot[0]]
    padY = [img.shape[0] - pivot[1], pivot[1]]
    imgP = np.pad(img, [padY, padX], 'constant')
    imgR = ndimage.rotate(imgP, angle, order=5, reshape=False)
    return imgR[padY[0] : -padY[1], padX[0] : -padX[1]]
    
def getHIcoincidence(x, y, object, ra_corner, dec_corner, height, width, dm):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    import matplotlib.colors as colors
    uchvcdb = os.path.dirname(os.path.abspath(__file__))+'/predblist.sort.csv'
    name, coords, ar, br, par = np.loadtxt(uchvcdb, usecols=(1,2,14,15,16), dtype=str, delimiter=',', unpack=True)
    # print object
    # find the right row
    coord = [this for i,this in enumerate(coords) if object.upper() in name[i]][0]
    a = [this for i,this in enumerate(ar) if object.upper() in name[i]][0]
    b = [this for i,this in enumerate(br) if object.upper() in name[i]][0]
    pa = [this for i,this in enumerate(par) if object.upper() in name[i]][0]
    
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
    
    cosd = lambda x : np.cos(np.deg2rad(x))
    sind = lambda x : np.sin(np.deg2rad(x))
    
    hi_c_x, hi_c_y = np.zeros(1), np.zeros(1)
    
    hi_c_x[0], hi_c_y[0] = abs((ra_hi-ra_corner)*60), abs((dec_hi-dec_corner)*60)
    
    a, b, pa = a.astype(float)/2., b.astype(float)/2., -pa.astype(float)
    
    t = np.array(list(range(0,359,1)))
    ell = np.array([a*cosd(t) , b*sind(t)])
    rot = np.array([[cosd(pa) , -sind(pa)],[sind(pa) , cosd(pa)]])
    ell_rot = np.zeros((2,ell.shape[1]))
    for i in range(ell.shape[1]):
        ell_rot[:,i] = np.dot(rot,ell[:,i])
    hi_x_circ, hi_y_circ = hi_c_x[0]+ell_rot[0,:], hi_c_y[0]+ell_rot[1,:]

    # x_dist = np.random.multivariate_normal((hi_c_x, hi_c_y), cov, 10000) 
     
    bins_h = int(height * 60. / 8.)
    bins_w = int(width * 60. / 8.)
    grid, xedges, yedges = np.histogram2d(hi_c_y, hi_c_x, bins=[bins_h,bins_w], range=[[0,height],[0,width]], normed=True) 
    xcenters = (xedges[:-1] + xedges[1:])/2.
    ycenters = (yedges[:-1] + yedges[1:])/2.
    # find the pivot point for the rotation == the HI centroid pixel from the histogram
    pivot = np.unravel_index(grid.argmax(),grid.shape)
    
    # convolve the single point with a 2d gaussian w/ a, b as axis ratios
    grid_gaus = ndimage.filters.gaussian_filter(grid, sigma=((bins_w/width)*2.*a/2.355, (bins_w/width)*2.*b/2.355), mode='nearest', cval=0)
    
    # rotate the image with the correct pivot point
    grid_rot = rotateImage(grid_gaus/np.amax(grid_gaus), pa, [pivot[-1], pivot[0]])
    
    if grid_rot[x][y] > 0.9:
        plt.clf()
        plt.figure(figsize=(5.5,5))
        bounds = np.linspace(0, 1, 11)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        
        extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
        gr = plt.imshow(grid_rot, norm=norm, extent=extent, interpolation='nearest')
        plt.plot(hi_x_circ,hi_y_circ,linestyle='-', color='limegreen')
        plt.scatter(ycenters[y], xcenters[x], c='red')
        plt.xlim(0,20)
        plt.ylim(0,20)
        plt.xlabel('RA (arcmin)')
        plt.ylabel('Dec (arcmin)')
        plt.title('{:s} @ dm = {:5.2f} : {:6.3f}%'.format(object, dm, grid_rot[x][y]*100.))
        plt.colorbar(gr, orientation='vertical')
        plt.savefig('{:s}_{:5.2f}_coinc.pdf'.format(object,dm))
    
    return grid_rot[x][y]    
    
def dist2HIcentroid(ra, dec, ra_hi, dec_hi, distance):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    c_hi = SkyCoord(ra = ra_hi, dec = dec_hi, distance=distance*1000000., unit=(u.deg, u.deg, u.pc))
    c_peak = SkyCoord(ra = ra, dec = dec, distance=distance*1000000., unit=(u.hourangle, u.deg, u.pc))
    sep = c_hi.separation(c_peak)
    sep3d = c_hi.separation_3d(c_peak)
    # print c_hi, c_peak, sep.arcsecond
    return sep.arcsecond, sep3d

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
        verts = list(zip(gi_iso,i_m_iso))        # set up the Path necessary for testing membership
        cm_filter = Path(verts)
        
    elif filter_file == 'none':
        gi_iso = None
        i_m_iso = None
        cm_filter = None
    else :
        g_iso,i_iso = np.loadtxt(filter_file, usecols=(8,10),unpack=True)
        gi_iso = g_iso - i_iso
        i_m_iso = i_iso + dm
        verts = list(zip(gi_iso,i_m_iso))        # set up the Path necessary for testing membership
        cm_filter = Path(verts)     # scale the filter to the DM: to find the apparent mag. just add the DM
    return cm_filter, gi_iso, i_m_iso

def make_youngpop(dm, filter_file):
    g_iso,i_iso = np.loadtxt(filter_file, usecols=(8,10),unpack=True)
    gi_iso = g_iso - i_iso
    i_m_iso = i_iso + dm
    # scale the filter to the DM: to find the apparent mag. just add the DM
    return gi_iso, i_m_iso

def filter_sources(i_mag, i_ierr, gmi, gmi_err, cm_filter, filter_sig = 1, gal = False):
    if cm_filter == None:
        stars_f = [True] * len(gmi)
    # first get the stars that are in the filter (points & 1-sig error bars!)
    elif gal:
        print(min(cm_filter.vertices[:,1]))
        stars_f = (0.75 < gmi) & (1.5 > gmi) & (min(cm_filter.vertices[:,1]) < i_mag)
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
    # print bins_h, bins_w
    density = float(len(i_ra_f))/(float(bins_h)*float(bins_w))
    
    grid, xedges, yedges = np.histogram2d(i_dec_f, i_ra_f, bins=[bins_h,bins_w], range=[[0,height],[0,width]])
    hist_points = list(zip(xedges,yedges))

    sig = ((bins_w/width)*fwhm)/2.355
    sig3 = ((bins_w/width)*3.0)/2.355
    pltsig = fwhm/2.0
    pgrid_mean = np.mean(grid)
    pgrid_sigma = np.std(grid)
    
    # convolve the grid with a gaussian
    grid_gaus = ndimage.filters.gaussian_filter(grid, sig, mode='constant', cval=0)
    S = np.array(grid_gaus*0)
    S_th = 3.0
    
    grid_mean = np.mean(grid_gaus)
    grid_sigma = np.std(grid_gaus)
    S = (grid_gaus-grid_mean)/grid_sigma
    # print len(i_ra_f), grid_mean, grid_sigma
    
    above_th = [(int(i),int(j)) for i in range(len(S)) for j in range(len(S[i])) if (S[i][j] >= S_th)]
    
    try: 
        segm = detect_sources(S, 2.0, npixels=5)
        props = source_properties(S, segm)
        columns = ['id', 'maxval_xpos', 'maxval_ypos', 'max_value', 'area']
        tbl = props.to_table(columns=columns)
    except ValueError:
        tbl = []
    # print tbl
    # rand_cmap = random_cmap(segm.max + 1, random_state=12345)
    
    # find the maximum point in the grid and center the circle there
    x_cent, y_cent = np.unravel_index(grid_gaus.argmax(),grid_gaus.shape)
    x_cent_S, y_cent_S = np.unravel_index(S.argmax(),S.shape)
    # print 'Max of S located at:','('+'{0:6.3f}'.format(y_cent_S)+','+'{0:6.3f}'.format(x_cent_S)+')'
    # print 'Value of S at above:','{0:6.3f}'.format(S[x_cent_S][y_cent_S])
    # print 'Number of bins above S_th: {0:4d}'.format(len(above_th))
    # print density, pgrid_mean, pgrid_sigma, grid_mean, grid_sigma, x_cent, y_cent
    
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
        random_xy = list(zip(random_ra,random_dec))
        grid_r, xedges_r, yedges_r = np.histogram2d(random_dec, random_ra, bins=[bins_h,bins_w], range=[[0,height],[0,width]])
        hist_points_r = list(zip(xedges_r,yedges_r))
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

    if samples > 10001 and pct > 95.:
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
        
    return pct, bins, centers

def dm_sigplot(dms, sig_bins, sig_max, fwhm, title_string):    
        plt.clf()
        plt.figure(figsize=(9,4))
        
        # print sig_bins
        # print sig_max
        # print dms
        
        plt.imshow(np.transpose(sig_bins), cmap=plt.cm.Reds, extent=(22, 27, 22, 2), aspect='auto')#, origin=origin)
        plt.plot(dms, sig_max, linestyle='-', color='black', lw=0.5)
        # plt.colorbar()
        plt.ylabel('$\sigma$')
        plt.xlabel('distance modulus')
        plt.xlim(22,27)
        plt.ylim(2,6.5)
        # plt.clim(0,1.5)
        
        plt.savefig('significance_{:s}_{:3.1f}.pdf'.format(title_string,fwhm))
    
# def association(ra, dec, a, b, phi):
#     from scipy.stats import multivariate_normal
#     mu = [0, 0]
#     
#     rot = [[np.cosd(phi), -1*np.sind(phi)],
#            [np.sind(phi), np.cosd(phi)]]
#     e1 = rot*[a, 0]
#     e2 = rot*[0, b]
#     
#     cov = a * ((e1 * e1)/())
#     
#     # cov = [[2.0, 0.3], [0.3, 0.5]]
#     
#     
#     x, y = np.mgrid[-1:1:.01, -1:1:.01]
#     pos = np.empty(x.shape + (2,))
#     pos[:, :, 0] = x; pos[:, :, 1] = y
#     rv = multivariate_normal(mu, cov)
#     plt.contourf(x, y, rv.pdf(pos))
#     
#     
#     
#     return pct
    
def magfilter(fwhm, fwhm_string, dm, dm_string, filter_file, filter_string, dm2=0.0):
    # print "Getting fits files..."
    # Load the FITS header using astropy.io.fits
    for file_ in os.listdir("./"):
        if file_.endswith("i.fits"):
            fits_file_i = file_
 
    for file_ in os.listdir("./"):
        if file_.endswith("g.fits"):
            fits_file_g = file_
    
    # downloadSDSSgal(fits_file_g, fits_file_i)
              
    fits_i = fits.open(fits_file_i)
    # fits_g = fits.open(fits_file_g)
    # print "Opened fits files:",fits_file_g,"&",fits_file_i
    
    # objid = fits_i[0].header['OBJECT']
    title_string = fits_file_i.split('_')[0]       # get the name part of the filename.

    # set up some filenames
    # mag_file = 'calibrated_mags.dat'
    mag_file = "AGC249525_daophot.dat.cut"

    # read in magnitudes, colors, and positions(x,y)
    # gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir,fwhm_sr= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8,11),unpack=True)
    # gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
    idr,rar,decr,ixr,iyr,am_g,g_ir,g_ierrr,am_i,i_ir,i_ierrr,g_magr,i_magr,gmir,chi,sharp,ebv,gfwhmr,fwhm_sr = np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),unpack=True)
    gxr, gyr = ixr, iyr
    # print len(gxr), "total stars"
    # fwhm_sr = np.ones_like(gxr)
    # filter out the things with crappy color errors
    mag_error_cut = 0.99
    color_error_cut = np.sqrt(2.0)*mag_error_cut
    
    
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
    cutleft = [i for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))]
    with open('cutleft.txt', 'w+') as spud:
        for i,item in enumerate(cutleft):
            print(item, file=spud)
    
    i_ierrAVG, bedges, binid = ss.binned_statistic(i_mag,i_ierr,statistic='median',bins=10,range=[15,25])
    gmi_errAVG, bedges, binid = ss.binned_statistic(i_mag,gmi_err,statistic='median',bins=10,range=[15,25])
    
    bcenters = (bedges[:-1] + bedges[1:]) / 2
    bxvals = [3.75,3.75,3.75,3.75,3.75,3.75,3.75,3.75,3.75,3.75]
    # print bcenters
    # print i_ierrAVG
    # print gmi_errAVG
    # print len(gx), "after color+mag error cut"
    # nid = np.loadtxt(mag_file,usecols=(0,),dtype=int,unpack=True)
    pixcrd = list(zip(ix,iy))
    
    
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
    # fits_g.close()
    
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
            print('{0:12.4f} {1:12.4f} {2:10.5f} {3:9.5f} {4:8.2f} {5:8.2f} {6:8.2f}'.format(ixBright[i],iyBright[i], i_radBright[i], i_decdBright[i], g_magBright[i], i_magBright[i], g_magBright[i]-i_magBright[i]), file=f1)
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
            print('{0:12.4f} {1:12.4f} {2:10.5f} {3:9.5f} {4:8.2f} {5:8.2f} {6:8.2f}'.format(ixRed[i],iyRed[i], i_radRed[i], i_decdRed[i], g_magRed[i], i_magRed[i], g_magRed[i]-i_magRed[i]), file=f1)
        f1.close()
    
    if dm2 > 0.0 and filter_string != 'none':
        dms = np.arange(dm,dm2,0.01)
        search = open('search_{:3.1f}.txt'.format(fwhm),'w+')
    else:
        dms = [dm]
        search = open('spud.txt'.format(fwhm),'w+')
    
    sig_bins = []
    sig_cens = []
    sig_max = []
    
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
        
        xy_points = list(zip(i_ra,i_dec))
        
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
        
        # xedgesg, x_centg, yedgesg, y_centg, Sg, x_cent_Sg, y_cent_Sg, pltsigg, tblg = galaxyMap(fits_file_i, fwhm, dm, filter_file)
        
        xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl = grid_smooth(i_ra_f, i_dec_f, fwhm, width, height)
        # corr = signal.correlate2d(S, Sg, boundary='fill', mode='full')
        # print corr
        
        pct, d_bins, d_cens = distfit(n_in_filter,S[x_cent_S][y_cent_S],title_string,width,height,fwhm,dm)
        pct_hi = 0.0 #getHIcoincidence(x_cent_S, y_cent_S, title_string, ra_corner, dec_corner, width, height, dm)
        
        sig_bins.append(d_bins)
        sig_cens.append(d_cens)
        sig_max.append(S[x_cent_S][y_cent_S])
        
        if pct > 100 :
            pct, bj,cj = distfit(n_in_filter,S[x_cent_S][y_cent_S],title_string,width,height,fwhm,dm, samples=25000)
        
        # make a circle to highlight a certain region
        cosd = lambda x : np.cos(np.deg2rad(x))
        sind = lambda x : np.sin(np.deg2rad(x))
        x_circ = [yedges[y_cent] + 3.0*cosd(t) for t in range(0,359,1)]
        y_circ = [xedges[x_cent] + 3.0*sind(t) for t in range(0,359,1)]
        
        verts_circ = list(zip(x_circ,y_circ))
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
                print('{:8.4f} {:8.4f}'.format(rCentx, rCenty), file=rc)
        else :
            rCentx, rCenty = np.loadtxt('refCircle.center', usecols=(0,1), unpack=True)
                
        x_circr = [rCentx + 3.0*cosd(t) for t in range(0,359,1)]
        y_circr = [rCenty + 3.0*sind(t) for t in range(0,359,1)]
        
        verts_circr = list(zip(x_circr,y_circr))
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
        index_fc = [i for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
        with open('index_fc.txt', 'w+') as spud1:
            for i,item in enumerate(index_fc):
                print(item, file=spud1)
        
        # print len(i_mag_fc), 'filter stars in circle'
        
        # print 'max i mag in circle = ', min(i_mag_fc)
        
        rs = np.array([51, 77, 90, 180])
        for r in rs:
            x_circ = [yedges[y_cent] + r/60.*cosd(t) for t in range(0,359,1)]
            y_circ = [xedges[x_cent] + r/60.*sind(t) for t in range(0,359,1)]

            verts_circ = list(zip(x_circ,y_circ))
            circ_filter = Path(verts_circ)

            stars_circ = circ_filter.contains_points(xy_points)
            i_x_fc = [ix[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]
            i_y_fc = [iy[i] for i in range(len(i_mag)) if (stars_circ[i] and stars_f[i])]

            fcirc_file = 'circle'+repr(r)+'.txt'
            with open(fcirc_file,'w+') as f3:
                for i,x in enumerate(i_x_fc):
                    print(i_x_fc[i], i_y_fc[i], file=f3)
        
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
        
        hi_x_circ, hi_y_circ = getHIellipse(title_string, ra_corner, dec_corner)
        hi_c_ra, hi_c_dec = getHIellipse(title_string, ra_corner, dec_corner, centroid=True)        
        hi_pix_x,hi_pix_y = w.wcs_world2pix(hi_c_ra,hi_c_dec,1)
        
        sep, sep3d = dist2HIcentroid(ra_c_d, dec_c_d, hi_c_ra, hi_c_dec, mpc)
        # print hi_pix_x, hi_pix_y
        # print ra_cr, dec_cr
        # print ra_c, dec_c
        # print hi_c_ra, hi_c_dec
        
        print("m-M = {:5.2f} | d = {:4.2f} Mpc | α = {:s}, δ = {:s}, Δʜɪ = {:5.1f}' | N = {:4d} | σ = {:6.3f} | ξ = {:6.3f}% | η = {:6.3f}%".format(dm, mpc, ra_c_d, dec_c_d, sep/60., n_in_filter, S[x_cent_S][y_cent_S], pct, pct_hi*100.))
        print('{:5.2f} {:4.2f} {:s} {:s} {:5.1f} {:4d} {:6.3f} {:6.3f} {:6.3f}'.format(dm, mpc, ra_c_d, dec_c_d, sep/60., n_in_filter, S[x_cent_S][y_cent_S], pct, pct_hi*100.), file=search)        

        #iraf.imutil.hedit(images=fits_g, fields='PV*', delete='yes', verify='no')
        #iraf.imutil.hedit(images=fits_i, fields='PV*', delete='yes', verify='no') 
        if pct > 97.:
            with open(ds9_file,'w+') as ds9:
                print("fk5;circle({:f},{:f},2') # color=yellow width=2 label=ref".format(ra_cr, dec_cr), file=ds9)
                print("fk5;circle({:f},{:f},2') # color=magenta width=2 label=detection".format(ra_c, dec_c), file=ds9)
                print("fk5;circle({:f},{:f},0.5') # color=black width=2 label=HI".format(hi_c_ra, hi_c_dec), file=ds9)
            
            f1 = open(filter_reg, 'w+')
            for i in range(len(i_x_f)) :
                print('image;circle({0:8.2f},{1:8.2f},10) # color=green width=2 edit=0 move=0 delete=0 highlite=0'.format(i_x_f[i],i_y_f[i]), file=f1)
            f1.close()
                
            f1 = open(mark_file, 'w+')
            for i in range(len(i_x_f)) :
                print('{0:8.2f} {1:8.2f} {2:12.8f} {3:12.8f} {4:8.2f} {5:8.2f} {6:8.2f} {7:7.3f}'.format(i_x_f[i],i_y_f[i],i_rad_f[i],i_decd_f[i],i_mag_f[i],g_mag_f[i],gmi_f[i],fwhm_sf[i]), file=f1)
            f1.close()
            
            f2 = open(circ_file, 'w+')
            for i in range(len(i_x_c)) :
                print('{0:8.2f} {1:8.2f} {2:12.8f} {3:12.8f} {4:8.2f} {5:8.2f} {6:8.2f} {7:7.3f}'.format(i_x_c[i],i_y_c[i],i_rad_c[i],i_decd_c[i],i_mag_c[i],i_mag_c[i]+gmi_c[i],gmi_c[i],fwhm_sc[i]), file=f2)
            f2.close()
            
            with open(fcirc_file,'w+') as f3:
                for i,x in enumerate(i_x_fc):
                    print(i_x_fc[i], i_y_fc[i], fwhm_sfc[i], file=f3)
                    
            with open(circles_file,'w+') as f4:
                print(circ_pix_x, circ_pix_y, file=f4)
            
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
            
            # if disp_flag :
            #     from pyraf import iraf
            #     iraf.tv(_doprint=0)
            #     iraf.unlearn(iraf.tv.tvmark)
            #     iraf.tv.tvmark.setParam('label',"no")
            #     iraf.tv.tvmark.setParam('pointsize',7)
            #     iraf.tv.tvmark.setParam('mark',"circle")
            #     # iraf.tv.display(image=fits_file_g, frame=1)
            #     iraf.tv.display(image=fits_file_i, frame=1)
            #     iraf.tv.tvmark(frame=1, coords=mark_file, radii="14,15,16", color=207)
            #     iraf.tv.tvmark(frame=1, coords=circ_file, radii="20,21,22", color=209)
            #     iraf.tv.tvmark(frame=1, coords='brightStars22.reg', radii="26,27,28", color=205)
            #     iraf.tv.tvmark(frame=1, coords='redStars175.reg', radii="32,33,34", color=204)
            #     iraf.tv.tvmark(frame=1, coords=center_file, txsize=4, mark="plus", color=208, label="yes")
            #     iraf.tv.tvmark(frame=1, coords=center_file, radii=mark_radii, color=208)
            #     iraf.tv.tvmark(frame=1, coords=center_file, radii=anni_radii, color=208)
            #     iraf.tv.tvmark(frame=1, coords=center_file, radii=anno_radii, color=208)
            #     # iraf.tv.tvmark(frame=2, coords=mark_file, radii="14,15,16", color=207, label="yes")
            #     # iraf.tv.tvmark(frame=2, coords=circ_file, radii="20,21,22", color=209)
            #     # iraf.tv.tvmark(frame=2, coords=center_file, txsize=4, mark="plus", color=208, label="yes")
            #     # iraf.tv.tvmark(frame=2, coords=center_file, radii=mark_radii, color=208)
        
            # setup the pdf output
            pp = PdfPages('f_'+ out_file)
            plt.clf()
            plt.figure(figsize=(8,6))
            # plot
            # print "Plotting for m-M = ",dm
            ax0 = plt.subplot(2,2,1)
            plt.scatter(i_ra, i_dec,  color='black', marker='o', s=1, edgecolors='none')
            plt.plot(x_circ,y_circ,linestyle='-', color='magenta')
            plt.plot(x_circr,y_circr,linestyle='-', color='gold')
            plt.plot(hi_x_circ,hi_y_circ,linestyle='-', color='limegreen')
            # plt.scatter(i_ra_c, i_dec_c,  color='red', marker='o', s=3, edgecolors='none')
            plt.scatter(i_ra_f, i_dec_f,  c='red', marker='o', s=10, edgecolors='none')
            # plt.clim(0,2)
            # plt.colorbar()
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
            plt.scatter(gmi_f, i_mag_f,  color='red', marker='o', s=15, edgecolors='none')

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
            cbar_S.set_label('$\sigma$ from local mean')
            # cbar_S.tick_params(labelsize=10)
            plt.plot(x_circ,y_circ,linestyle='-', color='magenta')
            plt.plot(x_circr,y_circr,linestyle='-', color='gold')
            plt.plot(hi_x_circ,hi_y_circ,linestyle='-', color='limegreen')
            # X, Y = np.meshgrid(xedges,yedges)
            # ax3.pcolormesh(X,Y,grid_gaus)
            plt.xlabel('RA (arcmin)')
            plt.ylabel('Dec (arcmin)')
            plt.title('smoothed stellar density')
            # plt.ylabel('Dec (arcmin)')
            plt.xlim(0,max(i_ra))
            plt.ylim(0,max(i_dec))
            ax2.set_aspect('equal')
        
            # ax3 = plt.subplot(2,2,4)
            ax3 = plt.subplot2grid((2,4), (1,2))
            plt.scatter(gmi_c, i_mag_c,  color='black', marker='o', s=3, edgecolors='none')
            plt.scatter(gmi_fc, i_mag_fc,  color='red', marker='o', s=15, edgecolors='none')    
            plt.tick_params(axis='y',left='on',right='on',labelleft='off',labelright='off')
            ax0.yaxis.set_label_position('left')
            plt.title('detection')
            plt.xlabel('$(g-i)_0$')
            plt.ylabel('$i_0$')
            plt.ylim(25,15)
            plt.xlim(-1,4)
            # ax3.set_aspect(0.5)    
        
            ax4 = plt.subplot2grid((2,4), (1,3), sharey=ax3)
            plt.scatter(gmi_cr, i_mag_cr,  color='black', marker='o', s=3, edgecolors='none')
            plt.scatter(gmi_fcr, i_mag_fcr,  color='red', marker='o', s=15, edgecolors='none')    
            plt.tick_params(axis='y',left='on',right='on',labelleft='off',labelright='on')
            plt.title('reference')
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
    
    # make the overall significance plot
    if len(dms) > 1 :        
        dm_sigplot(dms, sig_bins, sig_max, fwhm, title_string)
    
    return dm, mpc, ra_c_d, dec_c_d, sep, n_in_filter, S[x_cent_S][y_cent_S], pct, pct_hi*100. 
    # if imexam_flag :
    #     from pyraf import iraf
    #     iraf.unlearn(iraf.tv.imexamine, iraf.rimexam)
    #     iraf.tv.rimexam.setParam('radius',int(fwhm_i))
    #     iraf.tv.rimexam.setParam('rplot',12.)
    #     iraf.tv.rimexam.setParam('fittype','gaussian')
    #     iraf.tv.rimexam.setParam('center','yes')
    #     iraf.tv.imexamine(input=fits_file_i, frame=2)

def main(argv):
    home_root = os.environ['HOME']
    # set defaults for command line flags
    imexam_flag = False
    disp_flag = False
    # fwhm = 3.0
    # fwhm_string = repr(fwhm)
    filter_file = os.path.dirname(os.path.abspath(__file__))+'/filter.txt'
    filter_string = 'old'
    dm2 = 0.0

    try:
        opts, args = getopt.getopt(argv,"h",["fwhm=","dm=","dm2="])
    except getopt.GetoptError:
        print('magfilter.py --fwhm=<fwhm in arcmin> --dm=<DM in mag> --dm=<DM in mag>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('magfilter.py --fwhm=<fwhm in arcmin> --dm=<DM in mag> --dm=<DM in mag>')
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
    magfilter(fwhm, fwhm_string, dm, dm_string, filter_file, filter_string, dm2=dm2)

if __name__ == "__main__":
    main(sys.argv[1:])    
