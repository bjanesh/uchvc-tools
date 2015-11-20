from astropy.io import fits
import astropy.stats as aps
import scipy.stats as sps
import numpy as np
import os
import aplpy
import matplotlib.pyplot as plt
from pyraf import iraf

# def cmask(index,radius,array):
#   a,b = index
#   nx,ny = array.shape
#   y,x = np.ogrid[-a:nx-a,-b:ny-b]
#   mask = x*x + y*y <= radius*radius
# 
#   return(sum(array[mask]))

# define the object name and filenames
title_string = os.getcwd()[-9:].upper()
fitsMask_g = title_string + '_g_shPhotCutout.fits'
fitsMask_i = title_string + '_i_shPhotCutout.fits'
coords_file = 'region_coords.dat.CHECK'

dm = 22.89
dist = pow(10,((dm + 5.)/5.)) # convert the dm to pc

# figure out what the aperture and sky region have to be
radius = ((((1.4*60)/206265)*420000)/dist)*206265 # in arcsec

rPixels = int(radius/0.11) # pODI pixel scale is 0.11"/px
print 'Aperture radius:', radius/60., rPixels, rPixels+36, rPixels+36+100

radiiAperture = repr(rPixels-3)+','+repr(rPixels-2)+','+repr(rPixels-1)+','+repr(rPixels)+','+repr(rPixels+1)+','+repr(rPixels+2)+','+repr(rPixels+3)

radiiSky = repr(rPixels+36-3)+','+repr(rPixels+36-2)+','+repr(rPixels+36-1)+','+repr(rPixels+36)+','+repr(rPixels+36+1)+','+repr(rPixels+36+2)+','+repr(rPixels+36+3)+','+repr(rPixels+36+100-3)+','+repr(rPixels+36+100-2)+','+repr(rPixels+36+100-1)+','+repr(rPixels+36+100)+','+repr(rPixels+36+100+1)+','+repr(rPixels+36+100+2)+','+repr(rPixels+36+100+3)

# display the masked image
iraf.tv.display(image=fitsMask_g, frame=1)
iraf.tv.display(image=fitsMask_i, frame=2)
# mark the phot and sky regions
iraf.unlearn(iraf.tv.tvmark)
iraf.tv.tvmark.setParam('label',"no")
iraf.tv.tvmark.setParam('pointsize',7)
iraf.tv.tvmark.setParam('mark',"circle")
iraf.tv.tvmark(frame=1, coords=coords_file, radii=radiiAperture, color=207)
iraf.tv.tvmark(frame=1, coords=coords_file, radii=radiiSky, color=205)
iraf.tv.tvmark(frame=2, coords=coords_file, radii=radiiAperture, color=207)
iraf.tv.tvmark(frame=2, coords=coords_file, radii=radiiSky, color=205)

# set up phot parameters
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
iraf.fitskypars.setParam('snreject',100)
iraf.fitskypars.setParam('dannulus',100.)
# iraf.datapars.setParam('fwhmpsf',6.557) # don't need to set this 
iraf.photpars.setParam('apertures',float(rPixels)) 
iraf.fitskypars.setParam('annulus',float(rPixels+36))

# do the actual phot
iraf.apphot.phot(image=fitsMask_g, coords=coords_file, output="mag_est_g.iraf.dat") 
iraf.apphot.phot(image=fitsMask_i, coords=coords_file, output="mag_est_i.iraf.dat") 

# open the fits files via the astropy FITSIO library
hdulist_g = fits.open(fitsMask_g)
hdulist_i = fits.open(fitsMask_i)

gain = hdulist_i[0].header['gain']
exptime = hdulist_i[0].header['exptime']

scidata_g = hdulist_g[0].data
scidata_i = hdulist_i[0].data
hdulist_g.close()
hdulist_i.close()

# print len(scidata_g[0])
# print scidata_g[1001][1001]

b, a = 1000,1000
n = 2001

# make the aperture mask
rAperture = 847
y,x = np.ogrid[-a:n-a, -b:n-b]
# print y,x
mask1 = np.sqrt(x*x + y*y) <= rAperture - 0.5
mask0 = np.sqrt(x*x + y*y) >= rAperture + 0.5
maskfunc = np.sqrt(x*x + y*y) > rAperture - 0.5 #and (np.sqrt(x*x + y*y) < rAperture + 0.5)
# apMask = np.zeros((n, n))
apMask = rAperture - np.sqrt(x*x + y*y) + 0.5
apMask[mask1] = 1.
apMask[mask0] = 0.

# plt.imshow(apMask, cmap='gray')
# plt.xlim(0,2000)
# plt.ylim(0,2000)
# plt.colorbar()
# plt.show()

# make the outer annulus mask
# rOuter = 983
# y,x = np.ogrid[-a:n-a, -b:n-b]
# mask = np.sqrt(x*x + y*y) <= rOuter #- 0.5
# annoMask = np.zeros((n, n))
# annoMask[mask] = 1.
# 
# # make the inner annulus mask
# rInner = 883
# y,x = np.ogrid[-a:n-a, -b:n-b]
# mask = np.sqrt(x*x + y*y) <= rInner #- 0.5
# anniMask = np.ones((n, n))
# anniMask[mask] = 0.
# 
# # combine the inner and outer annulus masks to get a single annulus
# annMask1 = annoMask * anniMask

# plt.imshow(annMask, cmap='gray')
# plt.xlim(0,2000)
# plt.ylim(0,2000)
# plt.colorbar()
# plt.show()

dataMask = (scidata_g > 0.0).astype(int) * apMask
# annMask = (scidata_g > 0.0).astype(int) * annMask1

# annMaskInvert = -1.*annMask + 1.
apMaskInvert = -1.*dataMask + 1.
nPixels = np.sum(dataMask)
# plt.imshow(dataMask, cmap='gray')
# plt.xlim(0,2000)
# plt.ylim(0,2000)
# plt.colorbar()
# plt.show()
# check the number of pixels
# print np.sum(dataMask), np.sum(annMask)


# get the counts in the aperture and annulus for the g image
ap_g = apMask*scidata_g
# ann_g = annMask*scidata_g
# apCounts_g = np.sum(ap_g)
apCounts_g = np.genfromtxt('mag_est_g.iraf.dat', usecols=(1,), unpack=True, skip_header=79)
# # print np.median(ann_g)
# data = np.ma.MaskedArray(ann_g, dataMask)
# data = np.ma.masked_values(data, 0.0)
# data_clip = aps.sigma_clip(data, sig=3.0, iters=None)
# goodvals = data_clip.data[~data_clip.mask]
# nSky = len(np.nonzero(goodvals)[0])

# IRAF sky values
annMedian_g, annStd_g, nSky = np.genfromtxt('mag_est_g.iraf.dat', usecols=(0, 1, 3), unpack=True, skip_header=77, skip_footer=2)

# print '# nPixelsAperture nPixelsSky ========================'
# print nPixels, nSky
print '# g ================================================='
# ann_gClip, low, high = sps.sigmaclip(data[np.nonzero(data)], low=3.01, high=3.01)
# print np.median(ann_gClip), low, high
# annMean_g, annMedian_g, annStd_g = aps.sigma_clipped_stats(ann_g, mask=annMaskInvert, mask_val=0.0, iters=None)
# apMean_g, apMedian_g, apStd_g = aps.sigma_clipped_stats(ap_g, mask=apMaskInvert, mask_val=0.0, iters=None)
print nPixels, nSky, annMedian_g, annStd_g
flux_g = (apCounts_g - annMedian_g*nPixels)

g_i = -2.5*np.log10(flux_g/exptime)

error_g = np.sqrt(flux_g / gain + nPixels * annStd_g**2 + nPixels**2 * annStd_g**2 / nSky)
merr_g = 1.0857 * error_g / flux_g


print apCounts_g, flux_g, g_i, merr_g

print '# i ================================================='
# get the counts in the aperture and annulus for the i image
ap_i = apMask*scidata_i
# ann_i = annMask*scidata_i
# data = np.ma.MaskedArray(ann_i, dataMask)
# data = np.ma.masked_values(data, 0.0)
# data_clip = aps.sigma_clip(data, sig=3.0, iters=None)
# goodvals = data_clip.data[~data_clip.mask]
# nSky = len(np.nonzero(goodvals)[0])

# apCounts_i = np.sum(ap_i)
apCounts_i = np.genfromtxt('mag_est_i.iraf.dat', usecols=1, unpack=True, skip_header=79)
# annMean_i, annMedian_i, annStd_i = aps.sigma_clipped_stats(ann_i, mask=annMaskInvert, mask_val=0.0, iters=None)
# apMean_i, apMedian_i, apStd_i = aps.sigma_clipped_stats(ap_i, mask=apMaskInvert, mask_val=0.0, iters=None)
annMedian_i, annStd_i, nSky = np.genfromtxt('mag_est_i.iraf.dat', usecols=(0, 1, 3), skip_header=77, skip_footer=2, unpack=True )
print nPixels, nSky, annMedian_i, annStd_i
flux_i = (apCounts_i - annMedian_i*nPixels)

i_i = -2.5*np.log10(flux_i/exptime)

error_i = np.sqrt(flux_i / gain + nPixels * annStd_i**2 + nPixels**2 * annStd_i**2 / nSky)
merr_i = 1.0857 * error_i / flux_i

print apCounts_i, flux_i, i_i, merr_i

cal_A_i = 0.0526
cal_A_g = 0.1023
kg = 0.200
ki = 0.058

amg = 1.082642493
ami = 1.208449087
mu_gi = 1.054904668
zp_gi = 0.5713562363
eps_gi = 0.006664449944
zp_i = 25.89392439

dm = 22.89

# g-i = mu_gi * (g0 - i0) + ZP_gi
# i = eps_gi * (g-i) + ZP_i
g0 = g_i - (kg*amg)
i0 = i_i - (ki*ami)
gmi = mu_gi*(g0-i0) + zp_gi
i_mag = i0 + eps_gi*gmi + zp_i #- cal_A_i 
g_mag = gmi + i_mag - cal_A_g - 0.7525749893
i_mag = i_mag - cal_A_i - 0.7525749893
gmi = g_mag - i_mag

# print g_mag, i_mag, gmi

# estimate an absolute magnitude based on the input distance
absMag = i_mag - dm # in i-band magnitudes

iAbsSolar = 4.57	# sloan i solar absolute magnitude from sparke & gallagher 2e
gAbsSolar =  5.12 # sloan g solar absolute magnitude from sparke & gallagher 2e
# estimate a luminosity from the absolute magnitude (M = M_sun,i - 2.5log10(L/Lsun))
lum_i = np.power(10,(absMag - iAbsSolar)/-2.5) # in solar luminosities

# get a stellar mass estimate using g-i color and the formula from Bell+2003 (ApJSS, 149, 289)
# note that we are using the i magnitude and g-i color, +  the corresponding coefficients a_i and b_i from Table 7 of Bell+2003
logml = -0.152 + 0.518*(gmi) # first log(M/L)
mtol = np.power(10,logml)				  # get the actual number (10^x)
stellarMass = mtol * lum_i				  # mass = M/L * L

# for i in range(len(i_mag)):
print '# i mag, g mag, (g-i)'
print i_mag,g_mag,gmi
print '# absolute i mag, L, M/L, M'
print absMag,lum_i,mtol,stellarMass

# circs = np.array([847,883,983])*0.11/3600.0
# 
# fig = aplpy.FITSFigure(fitsMask_i)
# fig.show_grayscale()
# xWorld, yWorld = fig.pixel2world(1001,1001)
# print xWorld, yWorld, circs
# fig.show_circles(np.array([xWorld,xWorld,xWorld]), np.array([yWorld,yWorld,yWorld]), circs, edgecolor=['yellow','limegreen','limegreen'], linewidth=1)
# # fig.show_circles(xWorld, yWorld, circs[1])
# # fig.show_circles(xWorld, yWorld, circs[2])
# # fig.set_theme('publication')
# fig.save('photCutout_g.pdf')