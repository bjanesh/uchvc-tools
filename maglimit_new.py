#! /usr/local/bin/python
import os
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from astropy.stats import sigma_clipped_stats
from photutils import detect_sources
from photutils.utils import random_cmap
from scipy.ndimage import binary_dilation
from photutils import source_properties, properties_table
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval, LinearStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from uchvc_cal import js_calibrate
import pyregion

title_string = 'AGC249525'
fits_g = 'AGC249525_g_sh.fits'
fits_i = 'AGC249525_i_sh.fits'

hdu_g = fits.open(fits_g)
hdu_i = fits.open(fits_i)
# 
# # first guess backgrounds
# bgs_g = 5.017
# bgs_i = 11.907
# bg_g = 193.252
# bg_i = 630.374
# 
# # first get a good value for the background in each image by using an object mask
# threshold_g = bg_g + (3. * bgs_g)
# threshold_i = bg_i + (3. * bgs_i)
# 
# segm_g = detect_sources(hdu_g[0].data, threshold_g, npixels=25)
# source_mask1 =  segm_g.data.astype(np.bool)
# selem = np.ones((30,30))    # dilate using a 25x25 box
# source_mask_g = binary_dilation(source_mask1, selem)
# # source_mask_g = source_mask1
# 
# labels = [1, 5, 20, 50, 75, 80]
# props = source_properties(hdu_g[0].data, segm_g, labels=labels)
# tbl = properties_table(props)
# print tbl
# 
# # rand_cmap = random_cmap(segm_g.max + 1, random_state=12345)
# # interval=ZScaleInterval()
# # vmin, vmax = interval.get_limits(hdu_g[0].data)
# # norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
# # plt.imshow(hdu_g[0].data, origin='lower', cmap='Greys_r', norm=norm)
# # plt.imshow(segm_g, origin='lower', cmap=rand_cmap, alpha=0.6)
# # plt.show()
# 
# hot_pix_mask_g = (hdu_g[0].data > 58000.0).astype(int)
# dead_pix_mask_g = (hdu_g[0].data < 1.0).astype(int)
# full_mask_g = source_mask_g + hot_pix_mask_g + dead_pix_mask_g
# 
# segm_i = detect_sources(hdu_i[0].data, threshold_i, npixels=25)
# source_mask1 =  segm_i.data.astype(np.bool)
# selem = np.ones((15,15))    # dilate using a 25x25 box
# source_mask_i = binary_dilation(source_mask1, selem)
# 
# hot_pix_mask_i = (hdu_i[0].data > 58000.0).astype(int)
# dead_pix_mask_i = (hdu_i[0].data < 1.0).astype(int)
# full_mask_i = source_mask_i + hot_pix_mask_i + dead_pix_mask_i
# 
# mean_g, median_g, std_g = sigma_clipped_stats(hdu_g[0].data, mask=full_mask_g, sigma=3.0, iters=1)
# mean_i, median_i, std_i = sigma_clipped_stats(hdu_i[0].data, mask=full_mask_i, sigma=3.0, iters=1)
# # 
# print mean_g, median_g, std_g
# print mean_i, median_i, std_i

r2 = pyregion.open('regions_phys.reg')
full_mask = r2.get_mask(hdu=hdu_g[0])
# full_mask_i = full_mask_g #r2.get_mask(hdu=hdu_i[0])

# plt.imshow(full_mask_g)
# plt.show()

median_g = 191.9227
median_i = 623.8584

positions = [(5819.6140, 5535.8705)]
radii = [91., 136., 182., 227., 273., 318., 364., 409., 455., 500., 545., 591., 636., 682., 727., 772., 818.]
apertures = [CircularAperture(positions, r=r) for r in radii]
annulus = CircularAnnulus(positions, r_in=900., r_out=1000.)
ann_mask = annulus.to_mask(method='exact')
ann_g = ann_mask[0].apply(hdu_g[0].data) * ann_mask[0].apply(np.abs(full_mask-1))
ann_i = ann_mask[0].apply(hdu_i[0].data) * ann_mask[0].apply(np.abs(full_mask-1))
ann_keep = np.where(ann_g != 0)
bkg_mean_g = np.median(ann_g[ann_keep])
bkg_mean_i = np.median(ann_i[ann_keep])

print bkg_mean_g, bkg_mean_i

aper_mask = apertures[0].to_mask(method='exact')
dbl_mask_aper = aper_mask[0].apply(hdu_g[0].data-bkg_mean_g, fill_value=-999.) * aper_mask[0].apply(np.abs(full_mask-1), fill_value=-999.)
plt.imshow(dbl_mask_aper)
plt.show()

mean, median, std = sigma_clipped_stats(dbl_mask_aper, mask_value=-999.)
print mean, median, std

# aperture_area = np.sum(dbl_mask_aper)

# ann_table_g = aperture_photometry(hdu_g[0].data, annulus, mask=full_mask)
# ann_table_i = aperture_photometry(hdu_i[0].data, annulus, mask=full_mask)
#
# bkg_mean_g = ann_table_g['aperture_sum'][0] / annulus.area()
# bkg_mean_i = ann_table_i['aperture_sum'][0] / annulus.area()
#
# print bkg_mean_g, bkg_mean_i

phot_table_g = aperture_photometry(hdu_g[0].data-bkg_mean_g, apertures, mask=full_mask)
phot_table_i = aperture_photometry(hdu_i[0].data-bkg_mean_i, apertures, mask=full_mask)

# bkg_sum_g = bkg_mean_g * aperture_area
# final_sum_g = phot_table_g['aperture_sum'] - bkg_sum_g
# phot_table_g['bgsub_aperture_sum'] = final_sum_g
#
# bkg_sum_i = bkg_mean_i * aperture_area
# final_sum_i = phot_table_i['aperture_sum'] - bkg_sum_i
# phot_table_i['bgsub_aperture_sum'] = final_sum_i

mags_g, mags_i = np.zeros_like(np.array(radii)), np.zeros_like(np.array(radii))
for i,rad in enumerate(radii):
	mags_g[i] = -2.5*np.log10(phot_table_g['aperture_sum_'+repr(i)][0]) + 2.5*np.log10(300.0)
	mags_i[i] = -2.5*np.log10(phot_table_i['aperture_sum_'+repr(i)][0]) + 2.5*np.log10(300.0)

print mags_g, mags_i

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

g_mag, i_mag, gmi, g_abs, i_abs = np.zeros_like(np.array(radii)), np.zeros_like(np.array(radii)), np.zeros_like(np.array(radii)), np.zeros_like(np.array(radii)), np.zeros_like(np.array(radii))

for i,mg in enumerate(g_0):
	color_guess = 0.0
	color_diff = 1.0
	while abs(color_diff) > tolerance:
	    g_cal = g_0[i] + eps_g*color_guess + zp_g
	    i_cal = i_0[i] + eps_i*color_guess + zp_i

	    color_new = g_cal - i_cal
	    color_diff = color_guess-color_new
	    color_guess = color_new
	    # print g_0[i], g_cal, i_0[i], i_cal, color_new

	g_mag[i] = g_cal - cal_A_g + -2.5*np.log10(2.0)
	i_mag[i] = i_cal - cal_A_i + -2.5*np.log10(2.0)
	gmi[i] = g_mag[i] - i_mag[i]
	g_abs[i] = g_mag[i]-dm
	i_abs[i] = i_mag[i]-dm

	print '{0:5.1f} {1:6.3f} {2:6.3f} {3:6.3f} {4:6.3f} {5:6.3f}'.format(radii[i]*0.11, g_mag[i],i_mag[i],gmi[i],g_abs[i],i_abs[i])

plt.clf()
plt.plot(np.array(radii)*0.11,g_mag, 'b-', label='g')
plt.plot(np.array(radii)*0.11,i_mag, 'r-', label='i')
plt.xlim(0,100)
plt.ylim(24,16)
plt.xlabel('radius (arcsec)')
plt.ylabel('apparent magnitude')
# plt.plot(np.array(radii),g_mag, 'b-')
# plt.plot(np.array(radii),g_mag, 'b-')
# plt.plot(np.array(radii),g_mag, 'b-')

plt.savefig('magnitudes.pdf')