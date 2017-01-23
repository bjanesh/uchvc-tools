#!/usr/bin/env python
import os
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from astropy.stats import sigma_clipped_stats
from photutils import source_properties, properties_table
import matplotlib.pyplot as plt
from uchvc_cal import js_calibrate
from pyraf import iraf

iraf.images(_doprint=0)
iraf.tv(_doprint=0)
iraf.ptools(_doprint=0)
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.photcal(_doprint=0)
iraf.apphot(_doprint=0)  
iraf.imutil(_doprint=0)

title_string = 'AGC249525'
fits_g = 'AGC249525_g_sh.fits'
fits_i = 'AGC249525_i_sh.fits'
coords_file = 'fc_list_old_3.0_26.07_AGC249525.dat'

# iraf phot first
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
iraf.fitskypars.setParam('dannulus',4.)
iraf.datapars.setParam('fwhmpsf',6.197) 
iraf.photpars.setParam('apertures',8.) 
iraf.fitskypars.setParam('annulus',10.)
iraf.apphot.phot(image=title_string+"_i_sh.fits", coords=coords_file, output="mag_test_i.dat")

txdump_out = open('phot_test_i.txdump','w+')
iraf.ptools.txdump(textfiles='mag_test_i.dat', fields="id,mag,merr,sum,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

iraf.datapars.setParam('fwhmpsf',6.197)

iraf.apphot.phot(image=title_string+"_g_sh.fits", coords=coords_file, output="mag_test_g.dat")

txdump_out = open('phot_test_g.txdump','w+')
iraf.ptools.txdump(textfiles='mag_test_g.dat', fields="id,mag,merr,sum,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

g_iraf, ge_iraf, gf_iraf, gsky_iraf = np.loadtxt('phot_test_g.txdump', usecols=(1,2,3,4), unpack=True)
i_iraf, ie_iraf, if_iraf, isky_iraf = np.loadtxt('phot_test_i.txdump', usecols=(1,2,3,4), unpack=True)

# now try python
x, y = np.loadtxt(coords_file, usecols=(0,1), unpack=True)
positions = np.array(zip(x,y))

hdu_g = fits.open(fits_g)
hdu_i = fits.open(fits_i)

apertures = CircularAperture(positions, r=8.)
annulus_apertures = CircularAnnulus(positions, r_in=10., r_out=14.)
print apertures.area()
ap_mask = apertures.to_mask(method='subpixel', subpixels=7)
dummy = np.ones_like(hdu_g[0].data)
ann_mask = annulus_apertures.to_mask(method='center')
ap_g = [m.apply(hdu_g[0].data) for i,m in enumerate(ap_mask)]
ap_i = [m.apply(hdu_i[0].data) for i,m in enumerate(ap_mask)]
area_g = [np.sum(m.apply(dummy)) for i,m in enumerate(ap_mask)]
area_i = [np.sum(m.apply(dummy)) for i,m in enumerate(ap_mask)]

print area_g, area_i
# plt.imshow(ap_g[0], interpolation='nearest')
# plt.show()
ann_g = [m.apply(hdu_g[0].data, fill_value=-999.) for i,m in enumerate(ann_mask)]
ann_i = [m.apply(hdu_i[0].data, fill_value=-999.) for i,m in enumerate(ann_mask)]

flux_g = np.array([np.sum(a) for j,a in enumerate(ap_g)])
flux_i = np.array([np.sum(a) for j,a in enumerate(ap_i)])

bkg_med_g = np.array([sigma_clipped_stats(a, iters=0, mask_value=0.)[1] for j,a in enumerate(ann_g)])
bkg_med_i = np.array([sigma_clipped_stats(a, iters=0, mask_value=0.)[1] for j,a in enumerate(ann_i)])

# bkg_med_g = np.array([sigma_clipped_stats(a, iters=10, mask_value=0.)[1] for j,a in enumerate(ann_g)])
# bkg_med_i = np.array([sigma_clipped_stats(a, iters=10, mask_value=0.)[1] for j,a in enumerate(ann_i)])

# for j,a in enumerate(ann_g):
# 	# ann_keep = np.where(ann_g != 0)
# 	# bkg_mean_g = np.median(ann_g[ann_keep])
# 	# bkg_mean_i = np.median(ann_i[ann_keep])
# 	mean, med_g, std = sigma_clipped_stats(ann_g[j], mask_value=-999.)
# 	mean, med_i, std = sigma_clipped_stats(ann_g[j], mask_value=-999.)

# print bkg_med_g, bkg_med_i

# phot_tbl_g = aperture_photometry(hdu_g[0].data, apertures)
# phot_tbl_g_sky = aperture_photometry(hdu_g[0].data, annulus_apertures)

# phot_tbl_i = aperture_photometry(hdu_i[0].data, apertures)
# phot_tbl_i_sky = aperture_photometry(hdu_i[0].data, annulus_apertures)

# bkg_mean_g = phot_tbl_g_sky['aperture_sum'] / annulus_apertures.area()
# bkg_mean_i = phot_tbl_i_sky['aperture_sum'] / annulus_apertures.area()
print apertures.area()
bkg_sum_g = bkg_med_g * apertures.area()
final_sum_g = flux_g - bkg_sum_g
# phot_tbl_g['residual_aperture_sum'] = final_sum_g

bkg_sum_i = bkg_med_i * apertures.area()
final_sum_i = flux_i - bkg_sum_i
# phot_tbl_i['residual_aperture_sum'] = final_sum_i

g_py = -2.5*np.log10(final_sum_g) + 2.5*np.log10(300.0)
i_py = -2.5*np.log10(final_sum_i) + 2.5*np.log10(300.0)

for i,m in enumerate(g_iraf):
	print gf_iraf[i], flux_g[i], gsky_iraf[i], bkg_med_g[i], if_iraf[i], flux_i[i],  isky_iraf[i], bkg_med_i[i]
	# print gf_iraf[i], flux_g[i], g_iraf[i], g_py[i], gsky_iraf[i]-bkg_med_g[i], if_iraf[i], flux_i[i], i_iraf[i], i_py[i], isky_iraf[i]-bkg_med_i[i]
