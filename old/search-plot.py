#! /usr/local/bin/python
import aplpy
import numpy as np
import matplotlib.pyplot as plt

x, y, s = np.loadtxt('AGC268074_search.txt', usecols=(0,1,10), unpack=True)

f = aplpy.FITSFigure('AGC268074_i_sh.fits')
f.set_theme('publication')
f.show_grayscale()
f.add_grid()
ccx,ccy = f.pixel2world(x, y)

s95 = np.where(s>=95)
s90 = np.where((s>=90) & (s<95))
s85 = np.where((s>=85) & (s<90))
s80 = np.where((s>=80) & (s<85))
s75 = np.where((s>=75) & (s<80))
slo = np.where(s<75)

f.show_circles(ccx[s95], ccy[s95], radius=0.003, linestyle='-', edgecolor='red')
f.show_circles(ccx[s90], ccy[s90], radius=0.0028, linestyle='-', edgecolor='orange')
f.show_circles(ccx[s85], ccy[s85], radius=0.0026, linestyle='-', edgecolor='gold')
# f.show_circles(ccx[s80], ccy[s80], radius=0.0024, linestyle='-', edgecolor='green')
# f.show_circles(ccx[s75], ccy[s75], radius=0.0022, linestyle='-', edgecolor='blue')
# f.show_circles(ccx[slo], ccy[slo], radius=0.002, linestyle='-', edgecolor='purple')

# plt.scatter(x,y,c=s, edgecolor='none')
# plt.xlim(0,11000)
# plt.ylim(0,11000)
plt.savefig('AGC268074_search.pdf')