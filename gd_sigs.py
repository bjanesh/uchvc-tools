# -*- coding: utf-8 -*-
import os
import sys
import warnings
import numpy as np
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import scipy.stats as ss
from collections import OrderedDict
from magfilter import deg2HMS, grid_smooth, getHIellipse, dist2HIcentroid, make_filter, filter_sources, distfit, dm_sigplot

def main():
    objects = OrderedDict([('AGC198511',25.30), ('AGC215417',22.72), ('HI1151+20',24.70), ('AGC238626',26.01), ('AGC249000',25.30), ('AGC249320',25.34), ('AGC249525',26.07), ('AGC258237',24.02), ('AGC268069',24.24), ('AGC198606',22.89), ('HI0959+19',23.08)])
    filter_file = os.path.dirname(os.path.abspath(__file__))+'/filter.txt'
    
    fwhm_sm = 2.0
    fwhm_sm_string = '2.0'
    
    # plt.clf()
    fig = plt.figure(figsize=(10,7))
    outer = gridspec.GridSpec(4,3, wspace=0.1, hspace=0.1)
    for i, obj in enumerate(objects.keys()):
        print obj
        inner = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[i], wspace=0.1, hspace=0.1)
        # set up some filenames
        folder = '/Volumes/galileo/uchvc/targets/'+obj.lower()+'/'
        mag_file = folder+'calibrated_mags.dat'
        fits_file_i = folder+obj+'_i.fits'                  
        fits_i = fits.open(fits_file_i)
        
        dm = objects[obj]
        
        # read in magnitudes, colors, and positions(x,y)
        # gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir,fwhm_sr= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8,11),unpack=True)
        gxr,gyr,g_magr,g_ierrr,ixr,iyr,i_magr,i_ierrr,gmir= np.loadtxt(mag_file,usecols=(0,1,2,3,4,5,6,7,8),unpack=True)
        # print len(gxr), "total stars"
        fwhm_sr = np.ones_like(gxr)
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
        
        fits_i.close()
        
        # split the ra and dec out into individual arrays and transform to arcmin from the corner
        i_ra = [abs((world[i,0]-ra_corner)*60) for i in range(len(world[:,0]))]
        i_dec = [abs((world[i,1]-dec_corner)*60) for i in range(len(world[:,1]))]
        # also preserve the decimal degrees for reference
        i_rad = [world[i,0] for i in range(len(world[:,0]))]
        i_decd = [world[i,1] for i in range(len(world[:,1]))]
        
        search = open(obj+'_search'+fwhm_sm_string+'.txt','w+')
        
        sig_bins = []
        sig_cens = []
        sig_max = []
        
        dms = np.arange(22.0,26.99,0.1)
        
        for dm in dms:
            mpc = pow(10,((dm + 5.)/5.))/1000000.
            dm_string = '{:5.2f}'.format(dm).replace('.','_')
        
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
            
            xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl = grid_smooth(i_ra_f, i_dec_f, fwhm_sm, width, height)
            pct, d_bins, d_cens = distfit(n_in_filter,S[x_cent_S][y_cent_S],obj,width,height,fwhm_sm,dm)
            
            sig_bins.append(d_bins)
            sig_cens.append(d_cens)
            sig_max.append(S[x_cent_S][y_cent_S])
            
            circ_c_x = ra_corner-(yedges[y_cent]/60.)
            circ_c_y = (xedges[x_cent]/60.)+dec_corner
            circ_pix_x, circ_pix_y = w.wcs_world2pix(circ_c_x,circ_c_y,1)
            ra_c, dec_c = w.all_pix2world(circ_pix_x, circ_pix_y,1)
            ra_c_d,dec_c_d = deg2HMS(ra=ra_c, dec=dec_c, round=False)
            
            hi_c_ra, hi_c_dec = getHIellipse(obj, ra_corner, dec_corner, centroid=True)
            sep = dist2HIcentroid(ra_c_d, dec_c_d, hi_c_ra, hi_c_dec)
            print 'm-M = {:5.2f} | d = {:4.2f} Mpc | α = {:s}, δ = {:s}, Δʜɪ = {:5.1f}" | N = {:4d} | σ = {:6.3f} | ξ = {:6.3f}%'.format(dm, mpc, ra_c_d, dec_c_d, sep, n_in_filter, S[x_cent_S][y_cent_S], pct)
            print >> search, '{:5.2f} {:4.2f} {:s} {:s} {:5.1f} {:4d} {:6.3f} {:6.3f}'.format(dm, mpc, ra_c_d, dec_c_d, sep, n_in_filter, S[x_cent_S][y_cent_S], pct)
        
        ax1 = plt.Subplot(fig, inner[0])
        ax1.imshow(np.transpose(sig_bins), cmap=plt.cm.Reds, extent=(22, 27, 22, 2))#, origin=origin)
        ax1.plot(dms, sig_max, linestyle='-', color='black', lw=0.5)
        ax1.set_ylabel('$\sigma$')
        ax1.set_xlabel('distance modulus')
        ax1.set_xlim(22,27)
        ax1.set_ylim(2,6.5)
        if obj.startswith('HI1151'):
            ax1.set_title('AGC219656', size='small')
        else:
            ax1.set_title(obj, size='small')
        ax1.set_aspect(0.4)
        fig.add_subplot(ax1)
                    
    outer.tight_layout(fig)
    plt.savefig('detections_significance.pdf')
        
    pass

if __name__ == '__main__':
    main()

