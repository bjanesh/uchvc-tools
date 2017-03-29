#! /usr/local/bin/python
import os
import numpy as np
from pyraf import iraf
from uchvc_cal import js_calibrate

iraf.images(_doprint=0)
iraf.tv(_doprint=0)
iraf.ptools(_doprint=0)
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.photcal(_doprint=0)
iraf.apphot(_doprint=0)  
iraf.imutil(_doprint=0)

epadu = 1.268899
title_string = 'AGC249525'
coords_file = 'region_coords.dat'

if not os.path.isfile(title_string+'_i_sh_masked.fits'):

#   ix,iy,imag = np.loadtxt('calibrated_mags.dat',usecols=(4,5,6),unpack=True)
#     with open('bright_stars.dat','w+') as f1:
#         for i in range(len(ix)):
#             if imag[i] < 18.0 :
#                 print >> f1, ix[i], iy[i], imag[i]

    # iraf.tv.display(image=title_string+'_i_sh.fits', frame=1)
    #
    # iraf.unlearn(iraf.tv.tvmark)
    # iraf.tv.tvmark.setParam('label',"no")
    # iraf.tv.tvmark.setParam('pointsize',7)
    # iraf.tv.tvmark.setParam('mark',"circle")
    # iraf.tv.tvmark(frame=1, coords='bright_stars.dat', radii="98,99,100,101,102,103", color=208)

    while not os.path.isfile('regions.txt') :
        print 'Mask out bright stars indicated and other obvious things and save as regions.txt'
        raw_input("Press Enter when finished:")

    iraf.images.imcopy(title_string+'_i_sh.fits',title_string+'_i_sh_masked.fits',verbose="yes")

    m3,m4,m5,m6 = np.loadtxt('regions.txt',usecols=(2,3,4,5),unpack=True)
    bgmean_i = 623.670898438
    print "i sky value is:",np.mean(bgmean_i)

    bgmean_g = 191.777572632
    print "g sky value is:",np.mean(bgmean_g)

    for i in range(len(m3)) :
        x1 = m3[i] - (m5[i]/2.)
        x2 = m3[i] + (m5[i]/2.)
        y1 = m4[i] - (m6[i]/2.)
        y2 = m4[i] + (m6[i]/2.)

        if (x1 < 0):
            x1 = 1
        if (y1 < 0):
            y1 = 1
        if (x2 > 11000):
            x2 = 11000
        if (y2 > 11000):
            y2 = 11000
        iraf.unlearn(iraf.imreplace)
        iraf.imutil.imreplace(images=title_string+"_i_sh_masked.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=0.0)

if not os.path.isfile(title_string+'_g_sh_masked.fits'):

    ix,iy,imag = np.loadtxt('calibrated_mags.dat',usecols=(4,5,6),unpack=True)
    with open('bright_stars.dat','w+') as f1:
        for i in range(len(ix)):
            if imag[i] < 18.0 :
                print >> f1, ix[i], iy[i], imag[i]

    iraf.images.imcopy(title_string+'_g_sh.fits',title_string+'_g_sh_masked.fits',verbose="yes")

    m3,m4,m5,m6 = np.loadtxt('regions.txt',usecols=(2,3,4,5),unpack=True)

    print "g sky value is:",np.mean(bgmean_g)

    for i in range(len(m3)) :
        x1 = m3[i] - (m5[i]/2.)
        x2 = m3[i] + (m5[i]/2.)
        y1 = m4[i] - (m6[i]/2.)
        y2 = m4[i] + (m6[i]/2.)

        if (x1 < 0):
            x1 = 1
        if (y1 < 0):
            y1 = 1
        if (x2 > 11000):
            x2 = 11000
        if (y2 > 11000):
            y2 = 11000
        iraf.unlearn(iraf.imreplace)
        iraf.imutil.imreplace(images=title_string+"_g_sh_masked.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=0.0)

if not os.path.isfile('ones_mask.fits'):
    iraf.images.imarith(title_string+'_g_sh.fits', '*', 0.0, 'zeros.fits',verbose="yes")
    iraf.images.imarith('zeros.fits', '+', 1.0, 'ones_mask.fits',verbose="yes")

    m3,m4,m5,m6 = np.loadtxt('regions.txt',usecols=(2,3,4,5),unpack=True)

    for i in range(len(m3)) :
        x1 = m3[i] - (m5[i]/2.)
        x2 = m3[i] + (m5[i]/2.)
        y1 = m4[i] - (m6[i]/2.)
        y2 = m4[i] + (m6[i]/2.)

        if (x1 < 0):
            x1 = 1
        if (y1 < 0):
            y1 = 1
        if (x2 > 11000):
            x2 = 11000
        if (y2 > 11000):
            y2 = 11000
        iraf.unlearn(iraf.imreplace)
        iraf.imutil.imreplace(images="ones_mask.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=0.0)
# iraf.tv.display(image=title_string+'_i_sh_masked.fits', frame=1)
#
# iraf.unlearn(iraf.tv.tvmark)
# iraf.tv.tvmark.setParam('label',"no")
# iraf.tv.tvmark.setParam('pointsize',7)
# iraf.tv.tvmark.setParam('mark',"circle")
# iraf.tv.tvmark(frame=1, coords=coords_file, radii="1633,1634,1635,1636,1637,1638,1639,1697,1698,1699,1700,1701,1702,1703,1797,1798,1799,1800,1801,1802,1803", color=208)
# 
iraf.unlearn(iraf.apphot.phot, iraf.datapars, iraf.photpars, iraf.centerpars, iraf.fitskypars)
iraf.apphot.phot.setParam('interactive',"no")
iraf.apphot.phot.setParam('verify',"no")
iraf.datapars.setParam('datamin',0.)
iraf.datapars.setParam('datamax',50000.)
iraf.datapars.setParam('gain',"gain")
iraf.datapars.setParam('ccdread',"rdnoise")
iraf.datapars.setParam('exposure',"exptime")
iraf.datapars.setParam('airmass',"airmass")
iraf.datapars.setParam('filter',"filter")
iraf.datapars.setParam('obstime',"time-obs")
iraf.datapars.setParam('sigma',"INDEF")
iraf.photpars.setParam('zmag',0.)
iraf.centerpars.setParam('calgorithm',"none")
iraf.centerpars.setParam('cbox',9.)
iraf.centerpars.setParam('maxshift',3.)
iraf.fitskypars.setParam('salgorithm',"median")
iraf.fitskypars.setParam('dannulus',100.)
iraf.datapars.setParam('fwhmpsf',6.197) 
iraf.photpars.setParam('apertures','409,500,591,682,773,818') 
iraf.fitskypars.setParam('annulus',450.)
#
iraf.apphot.phot(image=title_string+"_i_sh_masked.fits", coords=coords_file, output="mag_est_i.dat")

txdump_out = open('phot_region_i.txdump','w+')
iraf.ptools.txdump(textfiles='mag_est_i.dat', fields="id,sum,msky,stdev,nsky", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

iraf.datapars.setParam('fwhmpsf',6.197)

iraf.apphot.phot(image=title_string+"_g_sh_masked.fits", coords=coords_file, output="mag_est_g.dat")

txdump_out = open('phot_region_g.txdump','w+')
iraf.ptools.txdump(textfiles='mag_est_g.dat', fields="id,sum,msky,stdev,nsky", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

iraf.apphot.phot(image="ones_mask.fits", coords=coords_file, output="mag_area.dat")
txdump_out = open('area.txdump','w+')
iraf.ptools.txdump(textfiles='mag_area.dat', fields="id,sum", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

# calculate the magnitude cf. apphot.phot, but *manually* subtract the area determined by the ones mask
areas = np.loadtxt('area.txdump', usecols=(1,2,3,4,5,6), unpack=True)
sum_i = np.loadtxt('phot_region_i.txdump',usecols=(1,2,3,4,5,6),unpack=True)
sky_i, stdev_i, nsky_i = np.loadtxt('phot_region_i.txdump',usecols=(7,8,9),unpack=True)
sum_g = np.loadtxt('phot_region_g.txdump',usecols=(1,2,3,4,5,6),unpack=True)
sky_g, stdev_g, nsky_g = np.loadtxt('phot_region_g.txdump',usecols=(7,8,9),unpack=True)

fl_i = sum_i - areas * sky_i
mag_i = -2.5*np.log10(fl_i) + 2.5*np.log10(300.)
error_i = np.sqrt(fl_i/ epadu + areas * stdev_i**2 + areas**2 * stdev_i**2 / nsky_i)
merr_i = 1.0857 * error_i / fl_i

fl_g = sum_g - areas * sky_g
mag_g = -2.5*np.log10(fl_g) + 2.5*np.log10(300.)
error_g = np.sqrt(fl_g/ epadu + areas * stdev_g**2 + areas**2 * stdev_g**2 / nsky_g)
merr_g = 1.0857 * error_g / fl_g

print fl_i,mag_i,merr_i
print fl_g,mag_g,merr_g

flux_g, flux_i, merrs_g, merrs_i = [], [], [], []
rs = np.array([45, 55, 65, 75, 85, 90])
for r in rs:
    fcirc_file = 'circle'+repr(r)+'.txt'

    iraf.fitskypars.setParam('dannulus',10.)
    iraf.datapars.setParam('fwhmpsf',6.197) 
    iraf.photpars.setParam('apertures',7.) 
    iraf.fitskypars.setParam('annulus',10.)
    #
    iraf.apphot.phot(image=title_string+"_g_sh.fits", coords=fcirc_file, output="mag_min_g.dat")
    txdump_out = open('phot_indiv_g.txdump','w+')
    iraf.ptools.txdump(textfiles='mag_min_g.dat', fields="id,mag,merr,flux,area,stdev,nsky", expr='yes', headers='no', Stdout=txdump_out)
    txdump_out.close()

    iraf.fitskypars.setParam('dannulus',10.)
    iraf.datapars.setParam('fwhmpsf',6.197) 
    iraf.photpars.setParam('apertures',7.) 
    iraf.fitskypars.setParam('annulus',10.)
    #
    iraf.apphot.phot(image=title_string+"_i_sh.fits", coords=fcirc_file, output="mag_min_i.dat")
    txdump_out = open('phot_indiv_i.txdump','w+')
    iraf.ptools.txdump(textfiles='mag_min_i.dat', fields="id,mag,merr,flux,area,stdev,nsky", expr='yes', headers='no', Stdout=txdump_out)
    txdump_out.close()

    fluxes_i, areas_i, stdevs_i, nskys_i = np.loadtxt('phot_indiv_i.txdump',usecols=(3,4,5,6),unpack=True)
    fluxes_g, areas_g, stdevs_g, nskys_g = np.loadtxt('phot_indiv_g.txdump',usecols=(3,4,5,6),unpack=True)
    flux_i.append(np.sum(fluxes_i))
    flux_g.append(np.sum(fluxes_g))

    area_i = np.sum(areas_i)
    stdev_i = np.median(stdevs_i)
    nsky_i = np.sum(nskys_i)

    error_i = np.sqrt(np.sum(fluxes_i) / epadu + area_i * stdev_i**2 + area_i**2 * stdev_i**2 / nsky_i)
    merrs_i.append(1.0857 * error_i / np.sum(fluxes_i))

    area_g = np.sum(areas_g)
    stdev_g = np.median(stdevs_g)
    nsky_g = np.sum(nskys_g)

    error_g = np.sqrt(np.sum(fluxes_g) / epadu + area_g * stdev_g**2 + area_g**2 * stdev_g**2 / nsky_g)
    merrs_g.append(1.0857 * error_g / np.sum(fluxes_g))

# print fl_i, flux_i, merr_i
# print fl_g, flux_g, merr_g

mags_i = np.hstack((-2.5*np.log10(fl_i)+2.5*np.log10(300.0), -2.5*np.log10(np.array(flux_i))+2.5*np.log10(300.0)))
mags_g = np.hstack((-2.5*np.log10(fl_g)+2.5*np.log10(300.0), -2.5*np.log10(np.array(flux_g))+2.5*np.log10(300.0)))
me_i = np.hstack((merr_i, merrs_i))
me_g = np.hstack((merr_g, merrs_g))
# mags_i = -2.5*np.log10(flux_i)+2.5*np.log10(300.0)
# mags_g = -2.5*np.log10(flux_g)+2.5*np.log10(300.0)

# print mags_i, mags_g

eps_g, std_eps_g, zp_g, std_zp_g, eps_i, std_eps_i, zp_i, std_zp_i = js_calibrate(img1 = title_string+"_g_sh.fits", img2 = title_string+"_i_sh.fits", verbose=False)

# values determined by ralf/daniel @ wiyn
kg = 0.20
kr = 0.12
ki = 0.058

cal_A_g =  0.0568
cal_A_i =  0.0292

tolerance = 0.0001
g_0 = mags_g - kg*1.043537637
i_0 = mags_i - ki*1.037019711

dm = 26.07
i_sun = 4.58
m_hi = 3.2E6
rs = np.array([45, 55, 65, 75, 85, 90, 45, 55, 65, 75, 85, 90])
for i,r in enumerate(rs):
    color_guess = 0.0
    color_diff = 1.0
    while abs(color_diff) > tolerance:
    	g_cal = g_0[i] + eps_g*color_guess + zp_g
    	i_cal = i_0[i] + eps_i*color_guess + zp_i

    	color_new = g_cal - i_cal
    	color_diff = color_guess-color_new
    	color_guess = color_new
    	# print g_0[i], g_cal, i_0[i], i_cal, color_new

    g_mag = g_cal - cal_A_g
    i_mag = i_cal - cal_A_i
    gmi = g_mag - i_mag
    e_gmi = np.sqrt(me_g[i]**2 + me_i[i]**2)
    g_abs = g_mag-dm
    i_abs = i_mag-dm
    

    mtol = np.power(10,0.518*gmi-0.152)
    l_star = np.power(10,(i_sun-i_abs)/2.5)
    m_star = l_star*mtol
    hitostar = m_hi/m_star
    print '{:3d} {:5.2f} {:4.2f} {:5.2f} {:4.2f} {:4.2f} {:4.2f} {:5.2f} {:5.2f} {:4.2f} {:3.1e} {:3.1e} {:5.1f}'.format(r,g_mag,me_g[i],i_mag,me_i[i],g_mag-i_mag,e_gmi,g_abs,i_abs,mtol,l_star,m_star,hitostar)
