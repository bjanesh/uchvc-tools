from __future__ import print_function
import os
import sys
import warnings
import numpy as np
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.path import Path
import scipy.stats as ss
from collections import OrderedDict
from magfilter import deg2HMS, grid_smooth, getHIellipse, dist2HIcentroid, make_filter, make_youngpop, filter_sources

def main():
    objects = OrderedDict([('AGC249320',25.28), ('AGC258242',25.05), ('AGC268074',22.10), ('HI0959+19',23.08)])
    smooths = OrderedDict([('AGC249320',2.0),   ('AGC258242',3.0),   ('AGC268074',2.0),   ('HI0959+19',2.0)])
    filter_file = os.path.dirname(os.path.abspath(__file__))+'/filter.txt'
    young_file = os.path.dirname(os.path.abspath(__file__))+'/filter_young.txt'

    # plt.clf()
    fig = plt.figure(figsize=(8.75,8.75))
    outer = gridspec.GridSpec(4, 1, wspace=0.05, hspace=0.01)
    for i, obj in enumerate(objects.keys()):
        inner = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[i], wspace=0.05, hspace=0.05)
        # set up some filenames
        # if os.path.isfile('/data/uchvc/targets/'+obj.lower()+'/'):
        folder = '/data/uchvc/targets/'+obj.lower()+'/'
        # else:
            # folder = '/Volumes/discovery1/uchvc/targets/'+obj.lower()+'/'
        mag_file = folder+'calibrated_mags.dat'
        fits_file_i = folder+obj+'_i.fits'                  
        fits_i = fits.open(fits_file_i)
        
        dm = objects[obj]
        mpc = pow(10,((dm + 5.)/5.))/1000000.
        
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
        i_mag = np.array([i_magr[i] for i in range(len(gxr)) if (abs(gmi_errr[i] < color_error_cut and i_ierrr[i] < mag_error_cut))])
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
        
        fits_i.close()
        
        # split the ra and dec out into individual arrays and transform to arcmin from the corner
        i_ra = [abs((world[i,0]-ra_corner)*60) for i in range(len(world[:,0]))]
        i_dec = [abs((world[i,1]-dec_corner)*60) for i in range(len(world[:,1]))]
        # also preserve the decimal degrees for reference
        i_rad = [world[i,0] for i in range(len(world[:,0]))]
        i_decd = [world[i,1] for i in range(len(world[:,1]))]
        
        cm_filter, gi_iso, i_m_iso = make_filter(dm, filter_file)
        gi_young, i_m_young = make_youngpop(dm, young_file)
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
        
        fwhm = smooths[obj]
        xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl = grid_smooth(i_ra_f, i_dec_f, fwhm, width, height)
        
        # xedges, x_cent, yedges, y_cent, S, x_cent_S, y_cent_S, pltsig, tbl = grid_smooth(i_ra, i_dec, 2.0, width, height)
        
        cosd = lambda x : np.cos(np.deg2rad(x))
        sind = lambda x : np.sin(np.deg2rad(x))
        
        # hi_c_x, hi_c_y = abs((hi_c_ra-ra_corner)*60), abs((hi_c_dec-dec_corner)*60)
        # 
        # t = np.array(range(0,359,1))
        # ell = np.array([a*cosd(t) , b*sind(t)])
        # rot = np.array([[cosd(pa) , -sind(pa)],[sind(pa) , cosd(pa)]])
        # ell_rot = np.zeros((2,ell.shape[1]))
        # for i in range(ell.shape[1]):
        #     ell_rot[:,i] = np.dot(rot,ell[:,i])
        # hi_x_circ, hi_y_circ = hi_c_x+ell_rot[0,:], hi_c_y+ell_rot[1,:]
        # hi_x_circ = [hi_c_x + a*cosd(t) for t in range(0,359,1)]
        # hi_y_circ = [hi_c_y + b*sind(t) for t in range(0,359,1)]
        
        x_circ = [yedges[y_cent] + 3.0*cosd(t) for t in range(0,359,1)]
        y_circ = [xedges[x_cent] + 3.0*sind(t) for t in range(0,359,1)]
        circ_c_x = ra_corner-(yedges[y_cent]/60.)
        circ_c_y = (xedges[x_cent]/60.)+dec_corner
        circ_pix_x, circ_pix_y = w.wcs_world2pix(circ_c_x,circ_c_y,1)
        ra_c, dec_c = w.all_pix2world(circ_pix_x, circ_pix_y,1)
        ra_c_d,dec_c_d = deg2HMS(ra=ra_c, dec=dec_c, round=False)

        hi_x_circ, hi_y_circ = getHIellipse(obj, ra_corner, dec_corner)
        hi_c_ra, hi_c_dec = getHIellipse(obj, ra_corner, dec_corner, centroid=True)
        sep, sep3d = dist2HIcentroid(ra_c_d, dec_c_d, hi_c_ra, hi_c_dec, mpc)
        print(obj, sep, sep3d)
        
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
        if not os.path.isfile(folder+'refCircle.center'):
            rCentx = 16.0*np.random.random()+2.0
            rCenty = 16.0*np.random.random()+2.0
            with open('refCircle.center','w+') as rc:
                print('{:8.4f} {:8.4f}'.format(rCentx, rCenty), file=rc)
        else :
            rCentx, rCenty = np.loadtxt(folder+'refCircle.center', usecols=(0,1), unpack=True)
                
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
        
        # fig = plt.figure(figsize=(8.5,8.5))
        ax1 = plt.Subplot(fig, inner[0])
        if os.path.isfile(folder+'i_gmi_compl.gr.out'):
            gmiCompl, iCompl = np.loadtxt(folder+'i_gmi_compl.gr.out',usecols=(0,1),unpack=True)
            ax1.plot(gmiCompl,iCompl, linestyle='--', color='green')
        # if os.path.isfile(folder+'i_gmi_compl.gr2.out'):
        #     gmiCompl, iCompl = np.loadtxt(folder+'i_gmi_compl.gr2.out',usecols=(0,1),unpack=True)
        #     ax1.plot(gmiCompl,iCompl, linestyle='--', color='red')
        # if os.path.isfile(folder+'i_gmi_compl2.out'):
        #     iCompl,gmiCompl = np.loadtxt(folder+'i_gmi_compl2.out',usecols=(0,1),unpack=True)
        #     ax1.plot(gmiCompl,iCompl, linestyle='--', color='blue')
        
        ax1.plot(gi_iso,i_m_iso,linestyle='-', lw=0.5, color='blue')
        ax1.plot(gi_young,i_m_young,linestyle='--', lw=0.5, color='blue')
        ax1.scatter(gmi, i_mag,  color='black', marker='o', s=1, edgecolors='none')
        ax1.scatter(gmi_f, i_mag_f,  color='red', marker='o', s=5, edgecolors='none')
        ax1.errorbar(bxvals, bcenters, xerr=i_ierrAVG, yerr=gmi_errAVG, linestyle='None', color='black', capsize=0, ms=0)
        ax1.tick_params(labelsize='small')
        ax1.yaxis.set_label_position('left')
        ax1.set_xticks([-1,0,1,2,3,4])
        ax1.set_yticks([15, 17, 19, 21, 23, 25])
        ax1.set_ylabel('$i_0$', size='small')
        ax1.set_xlabel('$(g-i)_0$', size='small')
        ax1.set_ylim(25,15)
        ax1.set_xlim(-1,4)
        ax1.set_aspect(0.5)
        ax1.set_title('$m-M = ${:5.2f} | $d = ${:4.2f} Mpc'.format(dm, mpc), size='small')
        if obj.startswith('HI1151'):
            ax1.text(-0.4, 0.75, 'AGC219656', size='small', weight='bold', rotation=90, transform=ax1.transAxes)
        else:
            ax1.text(-0.4, 0.75, obj, size='small', weight='bold', rotation=90, transform=ax1.transAxes)
        fig.add_subplot(ax1)
        
        ax2 = plt.Subplot(fig, inner[1])
        extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
        im = ax2.imshow(S, extent=extent, interpolation='nearest',cmap=cm.gray)
        # cbar_S = plt.colorbar()
        # cbar_S.set_label('$\sigma$ from local mean')
        ax2.plot(hi_x_circ,hi_y_circ,linestyle='-', color='limegreen')
        ax2.plot(x_circ,y_circ,linestyle='-', color='magenta')
        ax2.plot(x_circr,y_circr,linestyle='-', color='gold')
        ax2.tick_params(labelsize='small')        
        ax2.yaxis.set_label_position('left')
        ax2.set_xticks([0,5,10,15,20])
        ax2.set_yticks([0,5,10,15,20])
        ax2.set_xlabel('RA (arcmin)', size='small')
        ax2.set_ylabel('Dec (arcmin)', size='small')
        ax2.set_title("stellar density map", weight='bold', size='small')
        ax2.set_xlim(0,max(i_ra))
        ax2.set_ylim(0,max(i_dec))
        # im.set_clim(-3, 3)
        ax2.set_aspect('equal')
        fig.add_subplot(ax2)

        ax3 = plt.Subplot(fig, inner[2])
        if os.path.isfile(folder+'i_gmi_compl.gr.out'):
            gmiCompl, iCompl = np.loadtxt(folder+'i_gmi_compl.gr.out',usecols=(0,1),unpack=True)
            ax3.plot(gmiCompl,iCompl, linestyle='--', color='green')
        ax3.plot(gi_iso,i_m_iso,linestyle='-', lw=0.5, color='blue')
        ax3.plot(gi_young,i_m_young,linestyle='--', lw=0.5, color='blue')
        ax3.scatter(gmi_c, i_mag_c,  color='black', marker='o', s=1, edgecolors='none')
        ax3.scatter(gmi_fc, i_mag_fc,  color='red', marker='o', s=5, edgecolors='none')
        ax3.errorbar(bxvals, bcenters, xerr=i_ierrAVG, yerr=gmi_errAVG, linestyle='None', color='black', capsize=0, ms=0)
        ax3.tick_params(labelsize='small')
        ax3.yaxis.set_label_position('left')
        ax3.set_xticks([-1,0,1,2,3,4])
        ax3.set_yticks([15, 17, 19, 21, 23, 25])
        ax3.set_ylabel('$i_0$', size='small')
        ax3.set_xlabel('$(g-i)_0$', size='small')
        ax3.set_ylim(25,15)
        ax3.set_xlim(-1,4)
        ax3.set_aspect(0.5)
        ax3.set_title('object circle', color='magenta', weight='bold', size='small')
        fig.add_subplot(ax3)

        ax4 = plt.Subplot(fig, inner[3])
        if os.path.isfile(folder+'i_gmi_compl.gr.out'):
            gmiCompl, iCompl = np.loadtxt(folder+'i_gmi_compl.gr.out',usecols=(0,1),unpack=True)
            ax4.plot(gmiCompl,iCompl, linestyle='--', color='green')
        ax4.plot(gi_iso,i_m_iso,linestyle='-', lw=0.5, color='blue')
        ax4.plot(gi_young,i_m_young,linestyle='--', lw=0.5, color='blue')
        ax4.scatter(gmi_cr, i_mag_cr,  color='black', marker='o', s=1, edgecolors='none')
        ax4.scatter(gmi_fcr, i_mag_fcr,  color='red', marker='o', s=5, edgecolors='none')
        ax4.errorbar(bxvals, bcenters, xerr=i_ierrAVG, yerr=gmi_errAVG, linestyle='None', color='black', capsize=0, ms=0)
        ax4.tick_params(labelsize='small')
        ax4.yaxis.set_label_position('left')
        ax4.set_xticks([-1,0,1,2,3,4])
        ax4.set_yticks([15, 17, 19, 21, 23, 25])
        ax4.set_ylabel('$i_0$', size='small')
        ax4.set_xlabel('$(g-i)_0$', size='small')
        ax4.set_ylim(25,15)
        ax4.set_xlim(-1,4)
        ax4.set_aspect(0.5)
        ax4.set_title('reference circle', color='gold', weight='bold', size='small')
        fig.add_subplot(ax4)
        fig.tight_layout()
    

    outer.tight_layout(fig)
    plt.savefig('marginals.pdf')
    # plt.savefig('{:s}.pdf'.format(obj))
        
    pass

if __name__ == '__main__':
    main()
