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

epsgi = -0.0102148
zpi = 25.9233501
mugi = 1.0821236 
zpgi = 0.5868422
download_sdss(fits_g, fits_i, gmaglim = 22.0)
eps_g, std_eps_g, zp_g, std_zp_g, eps_i, std_eps_i, zp_i, std_zp_i, gairmass, iairmass = js_calibrate(img1 = fits_g, img2 = fits_i, verbose=False)
# gairmass = 1.0435380
# iairmass = 1.0370200

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

xnew = np.arange(-4.0,-0.25, 0.01)
compg = fg(xnew)
compi = fi(xnew)

# convert inst. mags to calibrated
tolerance = 0.0001
g_magjs = []
i_magjs = []
g_0 = xnew - kg*gairmass
i_0 = xnew - ki*iairmass
gmi_0 = mugi*(g_0-i_0) + zpgi
i_cal = i_0 + epsgi*(gmi_0) + zpi
g_cal = gmi_0 + i_cal

g0 = g_i - kg*gairmass
i0 = i_i - ki*iairmass
gmi0 = mugi*(g0-i0) + zpgi
i_mag = i0 + epsgi*(gmi0) + zpi
g_mag = gmi0 + i_mag

for j,mag in enumerate(g_0):
    color_guess = 0.0
    color_diff = 1.0
    while abs(color_diff) > tolerance:
        g_caljs = g_0[j] + eps_g*color_guess + zp_g
        i_caljs = i_0[j] + eps_i*color_guess + zp_i
        color_new = g_caljs - i_caljs
        color_diff = color_guess-color_new
        color_guess = color_new
        # print j, g_cal, i_cal, color_new
    g_magjs.append(g_caljs)
    i_magjs.append(i_caljs)
g_magjs = np.array(g_magjs)
i_magjs = np.array(i_magjs)

gmagjs = []
imagjs = []
# g0 = g_i - kg*gairmass 
# i0 = i_i - ki*iairmass 
for j,mag in enumerate(g0):
    color_guess = 0.0
    color_diff = 1.0
    while abs(color_diff) > tolerance:
        gcaljs = g0[j] + eps_g*color_guess + zp_g
        icaljs = i0[j] + eps_i*color_guess + zp_i
        color_new = gcaljs - icaljs
        color_diff = color_guess-color_new
        color_guess = color_new
        print j, gcaljs, icaljs, color_new
    gmagjs.append(gcaljs)
    imagjs.append(icaljs)

gmagjs = np.array(gmagjs)
imagjs = np.array(imagjs)

# print g_magjs-g_cal, i_magjs-i_cal

plt.clf()
plt.scatter(g_mag, complg, c='blue')
plt.scatter(i_mag, compli, c='red')
plt.scatter(gmagjs, complg, c='cyan')
plt.scatter(imagjs, compli, c='magenta')
plt.plot(g_cal, fg(xnew), 'b-')
plt.plot(i_cal, fi(xnew), 'r-')
plt.plot(g_magjs, fg(xnew), 'c-')
plt.plot(i_magjs, fi(xnew), 'm-')
plt.xlabel('mag')
plt.ylabel('%')
plt.savefig('compl_curves.pdf')

# for j in range(len(g_mag)):
# print '{:6.3f} {:5.3f} {:6.3f} {:5.3f}'.format(g_mag[j], compg[j], i_mag[j], compi[j])
    
# create bins to loop over in i
# ibin = np.linspace(first_bin_center, first_bin_center+binsize*numibins, num=numibins+1, endpoint=True)
# num = np.zeros_like(ibin)
# compi_interp = np.zeros_like(ibin)
# compg_interp = np.zeros_like(ibin)
#
gmi_a = np.arange(-1.2,4.1,0.1)
c50_i = np.zeros_like(gmi_a)
c50 = np.zeros_like(gmi_a)
for n,gmi in enumerate(gmi_a):
    for k in range(len(i_cal)):
        i0cjs = i_magjs[k] - eps_i*gmi - zp_i
        g0cjs = gmi + i0cjs
        i0c = i_cal[k] - epsgi*gmi - zpi
        g0c = ((gmi - zpgi)/mugi) + i0c
        compi_interp = compi[k]
        if k==0:
            print '{:4.1f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f}'.format(gmi, i0c, g0c, g0c-i0c, i0cjs, g0cjs, g0cjs-i0cjs)# i0c-i0cjs, g0c-g0cjs)
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
            if (comp > 0.495 and comp <= .50):
                c50[n]=comp
                c50_i[n] = i_cal[k]
                # print '{:5.2f} {:5.2f} {:6.4f} {:6.4f} {:6.4f}'.format(c50_i[n], gmi, compi_interp, compg_interp, comp)
 