from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import os
def getVabs(g_mag, i_mag, dm):
    # read in the SDSS BVR data from Lupton (2005) via Stetson
    # see http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php (bottom option)
    # the os.path functions are looking for the source file in the same folder as the python script
    # this is so you can run things from different folders
    v_magr, g_magr, i_magr = np.loadtxt(os.path.dirname(os.path.abspath(__file__))+'/sdssBVR.dat', usecols=(4, 12, 16), unpack=True)
    good_g, good_i = np.loadtxt(os.path.dirname(os.path.abspath(__file__))+'/sdssBVR.dat', usecols=(21,23), dtype=bool, unpack=True)
    
    # select only the stars in a certain color range and with good photometry
    good = np.where((v_magr-g_magr < 0.5) & (v_magr-g_magr > -1.5) & good_g & good_i)
    # cut down the arrays keeping only the good stars
    v, g, i = v_magr[good], g_magr[good], i_magr[good]
    
    # i only really care what the V mag is, so i'm going to fit V-g vs. g-i and then solve later for V based on the color
    # p is the fit parameters
    p = np.polyfit(g-i, v-g, 3)
    
    # this is the fit function generated from the params
    fit = np.poly1d(p)
    
    # sort the g-i so a line will look nice
    xplt = np.sort(g-i)
    yplt = fit(xplt)
    
    # this is all just fitting diagnostic stuff, or if you need an error
    # residuals
    # find the predicted v-g values
    x_fit = g-i
    y_fit = fit(x_fit)
    y_data = v-g
    res = y_data - y_fit
    # rms
    rms = np.sqrt(np.sum(res*res)/(res.size-4))
    # variance
    var = np.sum((y_data-np.mean(y_data))**2)/(y_fit.size-1)
    # chi square
    chi_sq = np.sum(res*res/var)/(res.size-4)
    # print rms, chi_sq
    
    # diagnostic plot
    plt.scatter(g-i,v-g, edgecolors='none')
    plt.plot(xplt, yplt, c='red')
    plt.xlabel('g-i')
    plt.ylabel('v-g')
    plt.savefig('gi_to_v.pdf')
    
    # plug in the input values of g and i for a single object to get the V mag and return the absolute (shifted by the distance modulus)
    v_mag = g_mag + fit(g_mag - i_mag)
    return v_mag-dm
    
if __name__ == '__main__':
    # test values (g, i, dm)
    vabs = getVabs(16.5, 15.5, 23.5)
    print(vabs)