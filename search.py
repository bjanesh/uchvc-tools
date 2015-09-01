#! /usr/local/bin/python
import os, sys, getopt, warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.path import Path
from matplotlib import cm
from astropy import wcs
from astropy.io import fits
try :
    from scipy import ndimage
except ImportError :
    print 'bad import'

def main(argv):
    # bin the filtered stars into a grid with pixel size XXX
    # print "Binning for m-M =",dm
    bins = 165
    width = 22 
    home_root = os.environ['HOME']
    imexam_flag = False
    fwhm = 2.0        # in arcmin (7.5 pixels = 1 arcmin)
    disp_flag = False

    filter_file = home_root+'/uchvc/filter.txt'
    filter_string = 'old'

    

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
    title_string = fits_file_i[0:9]        # get the first 9 characters of the filename.

    # set up some filenames
    out_file = title_string + '_search.txt'
    mag_file = 'calibrated_mags.dat'
    # mark_file = 'f_list_' + filter_string + '_' + fwhm_string + '_' + dm_string + '_' + title_string + '.dat'
    # circ_file = 'c_list_' + filter_string + '_' + fwhm_string + '_' + dm_string + '_' + title_string + '.dat'


    # read in magnitudes, colors, and positions(x,y)
    gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
    # print len(gxr), "total stars"

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
    # print 'Corner RA:',ra_c_d,':: Corner Dec:',dec_c_d

    fwhm_i = fits_i[0].header['FWHMPSF']
    fwhm_g = fits_g[0].header['FWHMPSF']

    # print 'Image FWHM :: g = {0:5.3f} : i = {1:5.3f}'.format(fwhm_g,fwhm_i)

    fits_i.close()
    fits_g.close()

    # split the ra and dec out into individual arrays and transform to arcmin from the corner
    i_ra = [abs((world[i,0]-ra_c)*60) for i in range(len(world[:,0]))]
    i_dec = [abs((world[i,1]-dec_c)*60) for i in range(len(world[:,1]))]
    # also preserve the decimal degrees for reference
    i_rad = [world[i,0] for i in range(len(world[:,0]))]
    i_decd = [world[i,1] for i in range(len(world[:,1]))]

    # the color-magnitude filter we're going to use in abs. mag. 
    # based on Girardi ugriz isochrones, see Walsh et al. and Betsey's thesis chapter
    if filter_file.endswith("iso_filter.txt") :
        gi_iso,i_iso = np.loadtxt(filter_file, usecols=(0,1),unpack=True)
    else :
        g_iso,i_iso = np.loadtxt(filter_file, usecols=(8,10),unpack=True)
        gi_iso = g_iso - i_iso
    
    
    with open(out_file, 'w+') as f1:
        print >> f1, '# dm mpc x_cent y_cent ra_cent dec_cent S(x,y) N %'
        for dm in np.arange(22.80,23.40,0.01):
            mpc = pow(10,((dm + 5.)/5.))/1000000.
            # print '================================================================================'
            # print 'DM = {0:6.3f}, dist (Mpc) = {1:6.3f}'.format(dm, mpc)
            i_m_iso = i_iso + dm     # scale the filter to the DM: to find the apparent mag. just add the DM
        
            verts = zip(gi_iso,i_m_iso)        # set up the Path necessary for testing membership
            cm_filter = Path(verts)
        
            # calculate the color error by adding instrumental mag errors in quadrature (not correct?...)
            stars_f = list(gmi)
        
            # pad_dir = np.loadtxt(filter_file, usecols=(2,),dtype=str,unpack=True)
            # make color-magnitude ordered pairs (necessary for membership testing)
            filter_sig = 1.0
            for i in range(len(gmi)) :
                nsteps_color = int(abs((float(filter_sig)*gmi_err[i])//0.001))
                nsteps_mag = int(abs((float(filter_sig)*i_ierr[i])//0.001))
        
                if nsteps_color == 0 :
                    nsteps_color = 1
                if nsteps_mag == 0 :
                    nsteps_mag = 1
        
                cm_points_l = [(gmi[i]-0.01*j*gmi_err[i],i_mag[i]) for j in range(nsteps_color)]
                cm_points_r = [(gmi[i]+0.01*j*gmi_err[i],i_mag[i]) for j in range(nsteps_color)]
                cm_points_u = [(gmi[i],i_mag[i]-0.01*j*i_ierr[i]) for j in range(nsteps_mag)]
                cm_points_d = [(gmi[i],i_mag[i]+0.01*j*i_ierr[i]) for j in range(nsteps_mag)]
        
                stars_f_l = cm_filter.contains_points(cm_points_l)        
                stars_f_r = cm_filter.contains_points(cm_points_r)        
                stars_f_u = cm_filter.contains_points(cm_points_u)        
                stars_f_d = cm_filter.contains_points(cm_points_d)        
        
                stars_f[i] = any(stars_f_l) | any(stars_f_r) | any(stars_f_u) | any(stars_f_d)  
        
            check = [stars_f[i] for i in range(len(stars_f)) if (stars_f[i])]
            # print filter_sig,len(check)
        
            xyd_points = zip(i_rad,i_decd)
            xy_points = zip(i_ra,i_dec)
        
            # make new vectors containing only the filtered points
        
            i_mag_f = [i_mag[i] for i in range(len(i_mag)) if (stars_f[i])]
            gmi_f = [gmi[i] for i in range(len(i_mag)) if (stars_f[i])]
            i_ra_f = [i_ra[i] for i in range(len(i_mag)) if (stars_f[i])]
            i_dec_f = [i_dec[i] for i in range(len(i_mag)) if (stars_f[i])]
            i_rad_f = [i_rad[i] for i in range(len(i_mag)) if (stars_f[i])]
            i_decd_f = [i_decd[i] for i in range(len(i_mag)) if (stars_f[i])]
            i_x_f = [ix[i] for i in range(len(i_mag)) if (stars_f[i])]
            i_y_f = [iy[i] for i in range(len(i_mag)) if (stars_f[i])]
            n_in_filter = len(i_mag_f)
        
    
        
            grid, xedges, yedges = np.histogram2d(i_dec_f, i_ra_f, bins=[bins,bins], range=[[0,width],[0,width]])
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
            #     for j in range(143) :
            #         ie = i + 22
            #         je = j + 22
            #         grid_mean = np.mean(grid_gaus[i:ie,j:je])
            #         grid_sigma = np.std(grid_gaus[i:ie,j:je])
            #         S[(11+i),(11+j)] = (grid_gaus[(11+i),(11+j)]-grid_mean)/grid_sigma
            #         # print i,j,grid_mean,grid_sigma,S[(22+i),(22+j)]
        
            # find the maximum point in the grid and center the circle there
            x_cent, y_cent = np.unravel_index(grid_gaus.argmax(),grid_gaus.shape)
            # print 'Max of gaussian convolved grid located at:','('+'{0:6.3f}'.format(yedges[y_cent])+','+'{0:6.3f}'.format(xedges[x_cent])+')'
            # # print grid_gaus[x_cent][y_cent]
            # print 'Value of S at above:','{0:6.3f}'.format(S[x_cent][y_cent])
        
            x_cent_S, y_cent_S = np.unravel_index(S.argmax(),S.shape)
            # print 'Max of S located at:','('+'{0:6.3f}'.format(yedges[y_cent_S])+','+'{0:6.3f}'.format(xedges[x_cent_S])+')'
            # print grid_gaus[x_cent][y_cent]
            # print 'Value of S at above:','{0:6.3f}'.format(S[x_cent_S][y_cent_S])
        
            # print 'Number of bins above S_th: {0:4d}'.format(len(above_th))
        
            sig_values_r = list()
            for i in range(100) :
                random_ra = max(i_ra)*np.random.random_sample((n_in_filter,))
                random_dec = max(i_dec)*np.random.random_sample((n_in_filter,))
                random_xy = zip(random_ra,random_dec)
                grid_r, xedges_r, yedges_r = np.histogram2d(random_dec, random_ra, bins=[bins,bins], range=[[0,width],[0,width]])
                hist_points_r = zip(xedges_r,yedges_r)
                grid_gaus_r = ndimage.filters.gaussian_filter(grid_r, sig, mode='constant', cval=0)
                S_r = np.array(grid_gaus_r*0)
        
                grid_mean_r = np.mean(grid_gaus_r)
                grid_sigma_r = np.std(grid_gaus_r)
                S_r = (grid_gaus_r-grid_mean_r)/grid_sigma_r
        
                above_th_r = [(int(i),int(j)) for i in range(len(S_r)) for j in range(len(S_r[i])) if (S_r[i][j] >= S_th)]
                x_cent_r, y_cent_r = np.unravel_index(grid_gaus_r.argmax(),grid_gaus_r.shape)
                sig_values_r.append(S_r[x_cent_r][y_cent_r])
        
            pct_calc = [sig_values_r[i] for i in range(len(sig_values_r)) if (sig_values_r[i] < S[x_cent][y_cent])]
            percentile = (float(len(pct_calc))/100.0)*100.0
        
            # print '================================================================================'
            # print 'Significance values:'
            # print 'Number of stars in filter:','{0:4d}'.format(n_in_filter)
            # print 'Mean value of peak S in random samples:','{0:6.3f}'.format(np.mean(sig_values_r))
            # print 'Sigma of peak S in random samples:','{0:6.3f}'.format(np.std(sig_values_r))
            # 
            # print 'Real value is{0:6.3f} sigma above mean'.format((S[x_cent][y_cent]-np.mean(sig_values_r))/np.std(sig_values_r))
            # print 'Percent samples with peak S below real value:','{0:6.3f}'.format(percentile)
            # print '================================================================================'
            
            circ_c_x = (yedges[y_cent]/60.)+ra_c
            circ_c_y = (xedges[x_cent]/60.)+dec_c
            circ_pix_x, circ_pix_y = w.wcs_world2pix(circ_c_x,circ_c_y,1)
            
            print >> f1, '{0:9.3f} {1:9.3f} {2:6.3f} {3:6.3f} {4:6.3f} {5:6.3f} {6:6.3f} {7:6.3f} {8:6.3f} {9:4d} {10:3d}'.format(100.0, 100.0, dm, mpc, yedges[y_cent_S], xedges[x_cent_S], circ_c_x, circ_c_y, S[x_cent_S][y_cent_S], n_in_filter, int(percentile))
    
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
