#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
def getVabs(g_mag, i_mag, dm):
    v_magr, g_magr, i_magr = np.loadtxt(os.path.dirname(os.path.abspath(__file__))+'/sdssBVR.dat', usecols=(4, 12, 16), unpack=True)
    good_g, good_i = np.loadtxt(os.path.dirname(os.path.abspath(__file__))+'/sdssBVR.dat', usecols=(21,23), dtype=bool, unpack=True)
    
    good = np.where((v_magr-g_magr < 0.5) & (v_magr-g_magr > -1.5) & good_g & good_i)
    # print good[0].size, v_magr.size
    v, g, i = v_magr[good], g_magr[good], i_magr[good]
    
    p = np.polyfit(g-i, v-g, 3)
    # print p
    
    fit = np.poly1d(p)
    
    xplt = np.sort(g-i)
    x_fit = g-i
    y_data = v-g
    
    yplt = fit(xplt)
    y_fit = fit(x_fit)
    
    res = y_data - y_fit
    rms = np.sqrt(np.sum(res*res)/(res.size-4))
    var = np.sum((y_data-np.mean(y_data))**2)/(y_fit.size-1)
    chi_sq = np.sum(res*res/var)/(res.size-4)
    # print rms, chi_sq
    
    # plt.scatter(g-i,v-g, edgecolors='none')
    # plt.plot(xplt, yplt, c='red')
    # plt.xlabel('g-i')
    # plt.ylabel('v-g')
    # plt.savefig('gi_to_v.pdf')
    
    # g_mag = np.array([g_mag])
    # i_mag = np.array([18.20, 20.93])
    v_mag = g_mag + fit(g_mag - i_mag)
    return v_mag-dm
    
