import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pyraf import iraf
import sewpy
import os

# findavgfwhm = sewpy.SEW(
#     params = ["X_IMAGE", "Y_IMAGE", "FWHM_IMAGE", "FLAGS"],
#     config = {"DETECT_THRESH":200.0},
#     sexpath = "sex"
# )
# 
# out = findavgfwhm("AGC174540_g_sh.fits")["table"]
# 
# fwhms = out['FWHM_IMAGE'] # This is an astropy table.
# flags = out['FLAGS']
# ap1x = np.median(fwhms[flags == 0])
# ap2x = 2.0*ap1x
# 
# fullfwhm = sewpy.SEW(
#     params = ["X_IMAGE","Y_IMAGE","FLUX_ISO","FLUXERR_ISO","MAG_APER(2)","MAGERR_APER(2)","FWHM_IMAGE","ELLIPTICITY","THETA_SKY","FLAGS","NUMBER"],
#     config = {"DETECT_THRESH":4.0, "PHOT_APERTURES":repr(ap1x)+', '+repr(ap2x)},
#     sexpath = "sex"
# )
# 
# out2 = fullfwhm("AGC174540_g_sh.fits")["table"]
# 
# mag1x = out2['MAG_APER']
# mag2x = out2['MAG_APER_1']
# flags = out2['FLAGS']
# idno = out2['NUMBER']


x, y, mag1x, mag2x, fwhm_sex, ellip = np.loadtxt('testfwhm.cat',usecols=(0,1,4,5,8,9), unpack=True)
flags, idno = np.loadtxt('testfwhm.cat', usecols=(11,12), dtype=int, unpack=True)
these = [ i for i,id in enumerate(idno) if (flags[i] == 0)] # and ellip[i] < 0.075

diff = mag2x - mag1x

# 
bin_meds, bin_edges, binnumber = stats.binned_statistic(mag2x[these], diff[these], statistic='median', bins=40, range=(-18,-8))
bin_stds, bin_edges, binnumber = stats.binned_statistic(mag2x[these], diff[these], statistic=np.std, bins=40, range=(-18,-8))
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2
print np.median(diff[these]), np.std(diff[these])
fwhm_iraf = np.loadtxt('getfwhm.log', usecols=(10,), unpack=True)
fwhm_med = np.median(fwhm_iraf)
# print len(these)
# # print np.mean(fwhm[these]), np.std(fwhm[these])
# with open('escut_g.pos', 'w+') as f1:
#     for j,id in enumerate(these):
#         print >> f1, x[id], y[id], fwhm[id]
# 
# iraf.images(_doprint=0)
# iraf.tv(_doprint=0)
# iraf.ptools(_doprint=0)
# iraf.noao(_doprint=0)
# iraf.digiphot(_doprint=0)
# iraf.photcal(_doprint=0)
iraf.apphot(_doprint=0)  
# iraf.imutil(_doprint=0)
# 
# iraf.rimexam.setParam('radius',4.0)
# iraf.rimexam.setParam('buffer',7.0)
# iraf.rimexam.setParam('width',5.0)
# iraf.rimexam.setParam('rplot',15.0)
# iraf.rimexam.setParam('center','no')
# iraf.rimexam.setParam('fittype',"gaussian")
# iraf.rimexam.setParam('iterati',1)
# iraf.tv.imexamine(image, logfile = 'getfwhm.log', keeplog = 'yes', defkey = "a", imagecur = 'escut_g.pos', wcs = "logical", use_display='no', frame='2')

iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
iraf.apphot.phot.setParam('interactive',"no")
iraf.apphot.phot.setParam('verify',"no")
iraf.datapars.setParam('datamax',50000.)
iraf.datapars.setParam('gain',"gain")
iraf.datapars.setParam('ccdread',"rdnoise")
iraf.datapars.setParam('exposure',"exptime")
iraf.datapars.setParam('airmass',"airmass")
iraf.datapars.setParam('filter',"filter")
iraf.datapars.setParam('obstime',"time-obs")
iraf.datapars.setParam('sigma',"INDEF")
iraf.photpars.setParam('zmag',0.)
iraf.centerpars.setParam('cbox',9.)
iraf.centerpars.setParam('maxshift',3.)
iraf.fitskypars.setParam('salgorithm',"median")
iraf.fitskypars.setParam('dannulus',10.)


if not os.path.isfile('sex.phot'):
    iraf.datapars.setParam('fwhmpsf',ap1x)
    iraf.photpars.setParam('apertures',repr(ap1x)+', '+repr(2.*ap1x))
    iraf.fitskypars.setParam('annulus',4.*ap1x)
    iraf.apphot.phot(image="AGC174540_g_sh.fits", coords='testfwhm.cat', output='testfwhm.phot.1')
    with open('sex.phot','w+') as txdump_out :
        iraf.ptools.txdump(textfiles='testfwhm.phot.1', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)


if not os.path.isfile('iraf.phot'):
    iraf.datapars.setParam('fwhmpsf',fwhm_med)
    iraf.photpars.setParam('apertures',repr(fwhm_med)+', '+repr(2.*fwhm_med))
    iraf.fitskypars.setParam('annulus',4.*fwhm_med)
    iraf.apphot.phot(image="AGC174540_g_sh.fits", coords='testfwhm.cat', output='testfwhm.phot.2')
    with open('iraf.phot','w+') as txdump_out :
        iraf.ptools.txdump(textfiles='testfwhm.phot.2', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)

m1xs, m2xs = np.loadtxt('sex.phot', usecols=(1,2), unpack=True)
m1xi, m2xi = np.loadtxt('iraf.phot', usecols=(1,2), unpack=True)

# p, pcov = np.polyfit(fwhm_sex[these], fwhm_iraf, 1, cov=True)
# perr = np.sqrt(np.diag(pcov))
# mu_gi, zp_gi, std_mu_gi, std_zp_gi = p[0], p[1], perr[0], perr[1]
# do a sigma clip based on the rms of the data from the first fit
# xplt = fwhm_sex[these]
# yplt = xplt
# 
# dy = yplt - fwhm_iraf
# print dy.mean(), np.median(fwhm_sex[these]), np.median(fwhm_iraf)
# 
# plt.scatter(fwhm_sex[these], fwhm_iraf, edgecolors='none', facecolors='black', s=4)
# plt.plot(xplt, yplt, 'r-', lw=3, alpha=0.7, label='fit')
# plt.fill_between(xplt, yplt+dy.std(), yplt-dy.std(), facecolor='blue', edgecolor='none', alpha=0.2, label='2x RMS sigma clipping region')

plt.clf()
plt.scatter(m2xs[these]-m1xs[these], m2xs[these], edgecolor='none', facecolor='black', s=4)
# plt.ylim(-1.5,-0.5)
# plt.xlabel('$m_{2xfwhm} - m_{1xfwhm}$')
# plt.ylabel('$m_{2xfwhm}$')
# plt.xlim(-1.5,-0.5)
plt.savefig('testmagsex.pdf')

plt.clf()
plt.scatter(m2xi[these]-m1xi[these], m2xi[these], edgecolor='none', facecolor='black', s=4)
plt.ylim(-6,-12)
# plt.xlabel('$m_{2xfwhm} - m_{1xfwhm}$')
# plt.ylabel('$m_{2xfwhm}$')
plt.xlim(-2,1)
plt.savefig('testmagiraf.pdf')

plt.clf()
plt.fill_betweenx(bin_centers, bin_meds+bin_stds, bin_meds-bin_stds, facecolor='red', edgecolor='none', alpha=0.2, label='2x RMS sigma clipping region')
plt.plot(bin_meds, bin_centers, 'r--', lw=2)
plt.scatter(mag2x[these]-mag1x[these], mag2x[these], edgecolor='none', facecolor='black', s=4)
plt.ylim(-8,-17)
plt.xlabel('$m_{2xfwhm} - m_{1xfwhm}$')
plt.ylabel('$m_{2xfwhm}$')
plt.xlim(-1.5,-0.5)
plt.savefig('test.pdf')

# plt.scatter(mag2x[those]-mag1x[those], mag2x[those], facecolors='red', edgecolors='none', s=4)
# # plt.scatter(fwhm[these],mag[these],edgecolors='none')
# # plt.scatter(fwhm[those],mag[those],edgecolors='red', facecolors='none')

# 
# plt.clf()
# plt.scatter(mag2x[these]-mag1x[these], ellip[these], edgecolors='none', s=2)
# plt.scatter(mag2x[those]-mag1x[those], ellip[those], facecolors='red', edgecolors='none', s=4)
# plt.xlim(-2,0)
# plt.show()