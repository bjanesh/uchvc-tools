import numpy as np
import matplotlib.pyplot as plt
from pyraf import iraf

ap_avg_g = 6.870
ap_avg_i = 6.422

# ap_avg_g = 7.0
# ap_avg_i = 7.0

# cal_A_i = 0.0526
# cal_A_g = 0.1023
# 
# apcor_g = -0.2853
# apcor_i = -0.4480
# 
kg = 0.200
ki = 0.058

amg = 1.063136
ami = 1.061535

# nid,gx,gy,g_i,g_ierr,ix,iy,i_i,i_ierr = np.loadtxt('calibration.dat',usecols=(0,1,2,4,5,11,12,14,15),unpack=True)

iraf.unlearn(iraf.phot,iraf.datapars,iraf.photpars,iraf.centerpars,iraf.fitskypars)
#
iraf.apphot.phot.setParam('interactive',"no")
iraf.apphot.phot.setParam('verify',"no")
iraf.datapars.setParam('datamin',"INDEF")
iraf.datapars.setParam('datamax',50000.)
iraf.datapars.setParam('gain',"gain")
iraf.datapars.setParam('ccdread',"rdnoise")
iraf.datapars.setParam('exposure',"exptime")
iraf.datapars.setParam('airmass',"airmass")
iraf.datapars.setParam('filter',"filter")
iraf.datapars.setParam('obstime',"time-obs")
iraf.datapars.setParam('sigma',"INDEF")
iraf.photpars.setParam('zmag',0.)
iraf.centerpars.setParam('calgorithm',"centroid")
iraf.centerpars.setParam('cbox',9.)
iraf.centerpars.setParam('maxshift',3.)
iraf.fitskypars.setParam('salgorithm',"median")
iraf.fitskypars.setParam('dannulus',3.)
#
# Use an aperture that is 1 x <fwhm>, because an aperture correction
# will be applied in the calc_calib_mags step
# Using a sky annulus thatbegins at 6 x <fwhm> should be fine
# g-band
# print 'Phot-ing g band point sources, this could take a while.'
iraf.datapars.setParam('fwhmpsf',ap_avg_g)
iraf.photpars.setParam('apertures',5*ap_avg_g)
iraf.fitskypars.setParam('annulus',5*ap_avg_g)
iraf.apphot.phot(image="AGC226067_g_sh.fits", coords='brightStars18.txt', output='magcheck5x_g.mag.1')

# print 'Phot-ing i band point sources, this could take a while.'
iraf.datapars.setParam('fwhmpsf',ap_avg_i)
iraf.photpars.setParam('apertures',5*ap_avg_i)
iraf.fitskypars.setParam('annulus',5*ap_avg_i)
iraf.apphot.phot(image="AGC226067_i_sh.fits", coords='brightStars18.txt', output='magcheck5x_i.mag.1')

txdump_out = open('magcheck5x_g.txdump','w+')
iraf.ptools.txdump(textfiles='magcheck5x_g.mag.1', fields="id,mag,merr,msky,stdev,flux,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

txdump_out = open('magcheck5x_i.txdump','w+')
iraf.ptools.txdump(textfiles='magcheck5x_i.mag.1', fields="id,mag,merr,msky,stdev,flux,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

g_i5x = np.loadtxt('magcheck5x_g.txdump',usecols=(1,),unpack=True)
i_i5x = np.loadtxt('magcheck5x_i.txdump',usecols=(1,),unpack=True)
# g_i5x2 = []
# i_i5x2 = []

# g_iSall = np.loadtxt('AGC198606_g.sdssphot',usecols=(16,),unpack=True)
# i_iSall = np.loadtxt('AGC198606_i_match.sdssphot',usecols=(16,),unpack=True)
# g_iS = []
# i_iS = []
# 
# raS, decS, gS = np.loadtxt('AGC198606_g.sdss',usecols=(0,1,4),unpack=True)
# 
# ixB, iyB, raB, decB, g_magB, i_magB = np.loadtxt('brightStars18.txt', usecols=(0,1,2,3,4,5), unpack=True)

# g_magB = []
# i_magB = []
# gmiB = []

ra,dec,psfMag_u,psfMag_g,psfMag_r,psfMag_i,psfMag_z = np.loadtxt('brightSDSS18.txt',usecols=(0,1,7,8,9,10,11),unpack=True)

# for i in range(len(raB)):
#     for j in range(len(ra)):
#         if abs(raB[i]-ra[j]) < 0.000000000001 and abs(decB[i]-dec[j]) < 0.000000000001 :
#             g_iS.append(g_iSall[i])
#             i_iS.append(i_iSall[i])
#             g_i5x2.append(g_i5x[j])
#             i_i5x2.append(i_i5x[j])
#             g_magB.append(g_magBr[j])
#             i_magB.append(i_magBr[j])
#             gmiB.append(gmiBr[j])
#             # print raS[i],ra[j],decS[i],dec[j],gS[i],psfMag_g[j],g_i5x[j],g_iSall[i]
# 
# print len(g_iS)
# AGC198606 specific values
mu_gi = 0.9982918697
zp_gi = 0.7184354406
eps_gi = 0.01738599227
zp_i = 25.77313787

# g-i = mu_gi * (g0 - i0) + ZP_gi
# i = eps_gi * (g-i) + ZP_i
g0A = g_i5x - (kg*amg)
i0A = i_i5x - (ki*ami)
gmiA = mu_gi*(g0A-i0A) + zp_gi
i_magA = i0A + eps_gi*gmiA + zp_i #- cal_A_i 
g_magA = gmiA + i_magA #- cal_A_g # get rid of the galactic extinction if you're  
# i_magA = i_magA - cal_A_i        # comparing to SDSS magnitudes
# gmiA = g_magA - i_magA
# g_magB = g_magB + cal_A_g 
# i_magB = i_magB + cal_A_i
# gmiB = g_magB - i_magB

# f_zp_g = open('AGC198606_g_phot.zp')
# data_g = f_zp_g.read()
# fl_g = data_g.split('\n', 1)[0]
# zp_vals_g = fl_g.split()
# cal_zp_g = float(zp_vals_g[2])
# cal_color_g = float(zp_vals_g[0])
# cal_zp_ge = float(zp_vals_g[3])
# cal_color_ge = float(zp_vals_g[1])
# 
# f_zp_i = open('AGC198606_i_phot.zp')
# data = f_zp_i.read()
# fl_i = data.split('\n', 1)[1]
# zp_vals_i = fl_i.split()
# cal_zp_i = float(zp_vals_i[0])
# cal_zp_ie = float(zp_vals_i[1])
# 
# g0B = g_i + apcor_g 
# i0B = i_i + apcor_i
# g_magB = g0B + cal_zp_g + cal_color_g*(g0B-i0B) - cal_A_g
# i_magB = i0B + cal_zp_i - cal_A_i
# gmiB = g_magB - i_magB

x = np.arange(-50,50,0.1)
y = x
# plt.scatter(gmiA, g-i,  color='blue', marker='o', s=10, edgecolors='none', label='model mag')
plt.scatter(gmiA, psfMag_g-psfMag_i,  color='red', marker='o', s=10, edgecolors='none', label='psf mag')
# plt.scatter(gmiA, petroMag_g-petroMag_i,  color='red', marker='o', s=10, edgecolors='none', label='petro mag')
plt.plot(x,y,linestyle='--', color='black')
plt.xlabel('$(g-i)_{pODI}$')
plt.ylabel('$(g-i)_{SDSS}$')
plt.title('SDSS vs. 5x aperture pODI colors')
# plt.legend(loc=2)
plt.xlim(0,3)
plt.ylim(0,3)

# print 'Zero point :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_g,cal_zp_i)
# print 'Zero point err :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_ge,cal_zp_ie)
# 
# print 'gi color term :: eps = {0:7.4f}'.format(cal_color_g)
# print 'gi color term err :: eps = {0:7.4f}'.format(cal_color_ge)

plt.savefig('gmicomp.pdf')

plt.clf()
# plt.scatter(i_magA, i,  color='blue', marker='o', s=10, edgecolors='none', label='model mag')
plt.scatter(i_magA, psfMag_i,  color='red', marker='o', s=10, edgecolors='none', label='psf mag')
# plt.scatter(i_magA, petroMag_i,  color='red', marker='o', s=10, edgecolors='none', label='petro mag')
plt.plot(x,y,linestyle='--', color='black')
plt.xlabel('$i_{pODI}$')
plt.ylabel('$i_{sdss}$')
plt.title('SDSS vs. 5x aperture pODI magnitudes')
# plt.legend(loc=2)
plt.xlim(15,21)
plt.ylim(15,21)
# plt.xlim(-11,-5)
# plt.ylim(-11,-5)

# print 'Zero point :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_g,cal_zp_i)
# print 'Zero point err :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_ge,cal_zp_ie)
# 
# print 'gi color term :: eps = {0:7.4f}'.format(cal_color_g)
# print 'gi color term err :: eps = {0:7.4f}'.format(cal_color_ge)

plt.savefig('icomp.pdf')

plt.clf()
# plt.scatter(g_magA, g,  color='blue', marker='o', s=10, edgecolors='none', label='model mag')
plt.scatter(g_magA, psfMag_g,  color='red', marker='o', s=10, edgecolors='none', label='psf mag')
# plt.scatter(g_magA, petroMag_g,  color='red', marker='o', s=10, edgecolors='none', label='petro mag')
plt.plot(x,y,linestyle='--', color='black')
plt.xlabel('$g_{pODI}$')
plt.ylabel('$g_{sdss}$')
plt.title('SDSS vs. 5x aperture pODI magnitudes')
# plt.legend(loc=2)
plt.xlim(15,21)
plt.ylim(15,21)
# plt.xlim(-11,-5)
# plt.ylim(-11,-5)

# print 'Zero point :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_g,cal_zp_i)
# print 'Zero point err :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_ge,cal_zp_ie)
# 
# print 'gi color term :: eps = {0:7.4f}'.format(cal_color_g)
# print 'gi color term err :: eps = {0:7.4f}'.format(cal_color_ge)

plt.savefig('gcomp.pdf')

# g_magB = g_magB - cal_A_g 
# i_magB = i_magB - cal_A_i
# gmiB = g_magB - i_magB
# 
# plt.clf()
# # plt.scatter(i_magA, i,  color='blue', marker='o', s=10, edgecolors='none', label='model mag')
# plt.scatter(i_magA, i_magB,  color='red', marker='o', s=10, edgecolors='none')
# # plt.scatter(i_magA, petroMag_i,  color='red', marker='o', s=10, edgecolors='none', label='petro mag')
# plt.plot(x,y,linestyle='--', color='black')
# plt.xlabel('$i_{5x}$')
# plt.ylabel('$i_{apcor}$')
# plt.title('5x vs. aperture corrected pODI magnitudes')
# # plt.legend(loc=2)
# plt.xlim(15,21)
# plt.ylim(15,21)
# # plt.xlim(-11,-5)
# # plt.ylim(-11,-5)
# 
# # print 'Zero point :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_g,cal_zp_i)
# # print 'Zero point err :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_ge,cal_zp_ie)
# # 
# # print 'gi color term :: eps = {0:7.4f}'.format(cal_color_g)
# # print 'gi color term err :: eps = {0:7.4f}'.format(cal_color_ge)
# 
# plt.savefig('icomp2.pdf')
# 
# plt.clf()
# # plt.scatter(g_magA, g,  color='blue', marker='o', s=10, edgecolors='none', label='model mag')
# plt.scatter(g_magA, g_magB,  color='red', marker='o', s=10, edgecolors='none')
# # plt.scatter(g_magA, petroMag_g,  color='red', marker='o', s=10, edgecolors='none', label='petro mag')
# plt.plot(x,y,linestyle='--', color='black')
# plt.xlabel('$g_{5x}$')
# plt.ylabel('$g_{apcor}$')
# plt.title('5x vs. aperture corrected pODI magnitudes')
# # plt.legend(loc=2)
# plt.xlim(15,21)
# plt.ylim(15,21)
# # plt.xlim(-11,-5)
# # plt.ylim(-11,-5)
# 
# # print 'Zero point :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_g,cal_zp_i)
# # print 'Zero point err :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_zp_ge,cal_zp_ie)
# # 
# # print 'gi color term :: eps = {0:7.4f}'.format(cal_color_g)
# # print 'gi color term err :: eps = {0:7.4f}'.format(cal_color_ge)
# 
# plt.savefig('gcomp2.pdf')