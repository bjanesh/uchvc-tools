#!/usr/bin/env python

import os
import numpy as np
from pyraf import iraf
import glob
from scipy import interpolate
import matplotlib.pyplot as plt
from uchvc_cal import download_sdss, js_calibrate

fits_g = 'AGC249525_g_sh.fits'
fits_i = 'AGC249525_i_sh.fits'

# Define constants
# extinction coefficients
kg = 0.2
ki = 0.058

# Bootstrap offsets - use if needed
g_offset = 0.0
i_offset = 0.0
#
download_sdss(fits_g, fits_i, gmaglim = 22.0)
eps_g, std_eps_g, zp_g, std_zp_g, eps_i, std_eps_i, zp_i, std_zp_i, gairmass, iairmass = js_calibrate(img1 = fits_g, img2 = fits_i, verbose=False)
gairmass = 1.082642
iairmass = 1.208449

# median g-i color for initial i magnitude calculation
medgmi = 1.236

# Define sampling of total convolved i magnitude completeness
# and location of center of first i fin
binsize = 0.01
first_bin_center = 22.0
numibins = 375

compl_binsize = 0.2

# fit the completeness data with a cubic spline so we can find the 50% values
# use iraf.curfit for legacy reasons
with open('g_results.out','w+') as f:
    iraf.curfit('ctable_g.out', function='spline3', order=5, interactive='no', listdata='yes', Stdout=f)
    
with open('i_results.out','w+') as f:
    iraf.curfit('ctable_i.out', function='spline3', order=5, interactive='no', listdata='yes', Stdout=f)

# from the completeness tables (ctable.out), read in inst. mags and completenesses
g_i, complg = np.loadtxt('g_results.out', usecols=(0,1), unpack=True)
i_i, compli = np.loadtxt('i_results.out', usecols=(0,1), unpack=True)

g_ift, complgft = np.loadtxt('g_completeness.fit_results', usecols=(0,1), unpack=True)
i_ift, complift = np.loadtxt('i_completeness.fit_results', usecols=(0,1), unpack=True)

fg = interpolate.interp1d(g_i, complg, kind=3)
fi = interpolate.interp1d(i_i, compli, kind=3)
fgft = interpolate.interp1d(g_ift, complgft, kind=3)
fift = interpolate.interp1d(i_ift, complift, kind=3)

xnew = np.arange(-3.5,-0.6, 0.01)
plt.clf()
plt.scatter(g_i, complg, c='blue')
plt.scatter(i_i, compli, c='red')
plt.scatter(g_ift, complgft, c='cyan')
plt.scatter(i_ift, complift, c='magenta')
plt.plot(xnew, fg(xnew), 'b-')
plt.plot(xnew, fi(xnew), 'r-')
plt.plot(xnew, fgft(xnew), 'c-')
plt.plot(xnew, fift(xnew), 'm-')
plt.xlabel('inst. mag')
plt.ylabel('%')
plt.show()
compg = fg(xnew)
compi = fi(xnew)

# convert inst. mags to calibrated
tolerance = 0.0001
g_mag = []
i_mag = []
g_0 = xnew - kg*gairmass
i_0 = xnew - ki*iairmass
for j,mag in enumerate(g_0):
    color_guess = 0.0
    color_diff = 1.0
    while abs(color_diff) > tolerance:
        g_cal = g_0[j] + eps_g*color_guess + zp_g
        i_cal = i_0[j] + eps_i*color_guess + zp_i
        
        color_new = g_cal - i_cal
        color_diff = color_guess-color_new
        color_guess = color_new
        # print j, g_cal, i_cal, color_new
    g_mag.append(g_cal)
    i_mag.append(i_cal)
    
g_mag = np.array(g_mag)
i_mag = np.array(i_mag)
g0 = g_i - kg*gairmass 
i0 = i_i - ki*iairmass

# for j in range(len(g_mag)):
# print '{:6.3f} {:5.3f} {:6.3f} {:5.3f}'.format(g_mag[j], compg[j], i_mag[j], compi[j])
    
# create bins to loop over in i
# ibin = np.linspace(first_bin_center, first_bin_center+binsize*numibins, num=numibins+1, endpoint=True)
# num = np.zeros_like(ibin)
# compi_interp = np.zeros_like(ibin)
# compg_interp = np.zeros_like(ibin)
#
gmi_a = np.arange(-1.2,3.1,0.1)
c50_i = np.zeros_like(gmi_a)
c50 = np.zeros_like(gmi_a)
for n,gmi in enumerate(gmi_a):
    for k in range(len(i_mag)):
        i0c = i_mag[k] - eps_i*gmi - zp_i
        g0c = gmi + i0c
        compi_interp = compi[k]
        # print gmi, i0c, g0c
#     for k in range(len(ibin)):
#         i0c = ibin[k] - eps_i*gmi - zp_i
#         g0c = gmi + i0c
        for j in range(1,len(compg)-1):
#             if ibin[k] < i_mag[0]:
#                 compi_interp[k] = 1.0
#             elif (ibin[k] > i_mag[j-1] and ibin[k] <= i_mag[j]):
#                 temp1 = (ibin[k] - i_mag[j-1]) / compl_binsize
#                 temp2 = (compi[j-1]-compi[j]) * temp1
#                 compi_interp[k] = compi[j-1] - temp2
#
            if (g0c < g_0[0]):
                compg_interp = 1.0
            elif ((g0c > g_0[j-1]) and (g0c <= g_0[j])):
                temp1 = (g0c - g_0[j-1]) / binsize
                temp2 = (compg[j-1]-compg[j]) * temp1
                compg_interp = compg[j-1] - temp2
            else:
                compg_interp = 0.0
            comp = compi_interp*compg_interp
            # print comp, compi_interp, compg_interp
            if (gmi < 2.0 and comp > 0.493 and comp < .501):
                c50[n]=comp
                c50_i[n] = i_mag[k]
                print '{:5.2f} {:5.2f} {:6.4f} {:6.4f} {:6.4f}'.format(c50_i[n], gmi, compi_interp, compg_interp, comp)
            elif (gmi >= 2.0 and comp > 0.3 and comp < .501):
                c50[n]=comp
                c50_i[n] = i_mag[k]
                print '{:5.2f} {:5.2f} {:6.4f} {:6.4f} {:6.4f}'.format(c50_i[n], gmi, compi_interp, compg_interp, comp)
 

# plt.clf()
# plt.scatter(i_mag, compi)
# plt.plot(ibin, compi_interp)
# plt.scatter(g_mag, compg)
# plt.plot(ibin, compg_interp)
# plt.show()