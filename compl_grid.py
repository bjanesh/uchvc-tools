import os
import sys
import numpy as np
from pyraf import iraf
import glob
from matplotlib import cm
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.special import erfc
import matplotlib.pyplot as plt
from odi_calibrate import download_sdss, js_calibrate, get_calibration

def cubic(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d

def erfc_p(x, a, b, m):
    return m*erfc((x-a)/(b*np.sqrt(2)))

def main():
    path = os.getcwd()
    steps = path.split('/')
    objname = steps[-1].upper()
    fits_g = objname+'_g.fits'
    fits_i = objname+'_i.fits'
    
    # Define constants
    # extinction coefficients
    kg = 0.2
    ki = 0.058
    # fit the completeness data with a cubic spline so we can find the 50% values
    # use iraf.curfit for legacy reasons
    with open('g_results.out','w+') as f:
        iraf.curfit('ctable_g.out', function='spline3', order=5, interactive='no', listdata='yes', Stdout=f)
    
    with open('i_results.out','w+') as f:
        iraf.curfit('ctable_i.out', function='spline3', order=5, interactive='no', listdata='yes', Stdout=f)
    
    # from the completeness tables (ctable.out), read in inst. mags and completenesses
    g_i, complg = np.loadtxt('g_results.out', usecols=(0,1), unpack=True)
    i_i, compli = np.loadtxt('i_results.out', usecols=(0,1), unpack=True)
    
    gi, complgu = np.loadtxt('ctable_g.out', usecols=(0,1), unpack=True)
    ii, compliu = np.loadtxt('ctable_i.out', usecols=(0,1), unpack=True)
    
    # print(g_i,i_i,complg,compli)

    fg = interpolate.interp1d(g_i, complg, kind=3)
    fi = interpolate.interp1d(i_i, compli, kind=3)
    
    xnew = np.arange(-6.0,0.89, 0.01)
    
    # download_sdss(fits_g, fits_i, gmaglim = 22.0)
    eps_g, std_eps_g, zp_g, std_zp_g, eps_i, std_eps_i, zp_i, std_zp_i = get_calibration()
    
    # get some auxiliary info from the phot output
    gXAIRMASS = np.loadtxt(fits_g[0:-5]+'_cal_js.sdssphot', usecols=(9,), dtype=str, unpack=True)
    iXAIRMASS = np.loadtxt(fits_i[0:-5]+'_cal_js.sdssphot', usecols=(9,), dtype=str, unpack=True)
    # keep the airmasses and aperture radii as single values
    gXAIRMASS, iXAIRMASS = gXAIRMASS.astype(float)[0], iXAIRMASS.astype(float)[0]
    
    # convert inst. mags to calibrated
    tolerance = 0.0001
    g_magjs = []
    i_magjs = []
    g_0 = xnew - kg*gXAIRMASS
    i_0 = xnew - ki*iXAIRMASS
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
    
    g_js = []
    i_js = []
    g0 = gi - kg*gXAIRMASS
    i0 = ii - ki*iXAIRMASS
    for j,mag in enumerate(g0):
        color_guess = 0.0
        color_diff = 1.0
        while abs(color_diff) > tolerance:
            g_cjs = g0[j] + eps_g*color_guess + zp_g
            i_cjs = i0[j] + eps_i*color_guess + zp_i
            color_new = g_cjs - i_cjs
            color_diff = color_guess-color_new
            color_guess = color_new
            # print j, g_cal, i_cal, color_new
        g_js.append(g_cjs)
        i_js.append(i_cjs)
    g_js = np.array(g_js)
    i_js = np.array(i_js)
    
    # print(g_js, i_js)
    
    p_g, pcov_g = curve_fit(erfc_p, g_js, complgu, p0=[25., 1.0, 1.0])
    p_i, pcov_i = curve_fit(erfc_p, i_js, compliu, p0=[25., 1.0, 1.0])
    
    # print(p_g)
    # print(p_i)
    
    magplot = np.arange(21.0, 27.0, 0.01)
    
    plt.clf()
    plt.scatter(g_js, complgu, c='blue')
    plt.scatter(i_js, compliu, c='red')
    plt.plot(g_magjs, fg(xnew), 'b-', label='odi_g')
    plt.plot(i_magjs, fi(xnew), 'r-', label='odi_i')
    plt.plot(magplot, erfc_p(magplot, p_g[0], p_g[1], p_g[2]), c='cyan')
    plt.plot(magplot, erfc_p(magplot, p_i[0], p_i[1], p_i[2]), c='magenta')
    plt.xlim(21,27)
    plt.ylim(-0.05, 1.05)
    plt.xlabel('calibrated magnitude')
    plt.ylabel('completeness %')
    plt.legend()
    plt.savefig('compl_curves.pdf')
    compg = fg(xnew)
    compi = fi(xnew)
    
    gmi_a = np.arange(-1.5,4.0,0.01)
    compi_gmi = np.zeros((i_magjs.size, gmi_a.size))
    for i,gmi in enumerate(gmi_a):
        for j,im in enumerate(i_magjs):
            gm = gmi+im
            minusg = np.absolute(g_magjs-gm)
            closestg = np.argmin(minusg)
            if gm > g_magjs[-1]:
                cg = 0.0
            elif gm < g_magjs[0]:
                cg = 1.0
            else:
                cg = compg[closestg]
            
            # print '{:5.2f} {:6.3f} {:6.3f} {:6.4f} {:6.4f} {:6.4f} {:6.4f}'.format(gmi, gm, g_magjs[closestg], im, compi[j], cg, compi[j]*cg)
            compi_gmi[j][i] = compi[j]*cg
    
    minuscomp = np.absolute(compi_gmi-0.5)
    index50 = np.argmin(minuscomp, axis=0)
    print(index50.size, i_magjs.size, gmi_a.size)
    line_gmi = []
    line_i = []
    with open('i_gmi_compl.gr.out','w+') as f:
        for k in range(gmi_a.size):
            if index50[k] > 0:
                line_gmi.append(gmi_a[k])
                line_i.append(i_magjs[index50[k]])
                print(gmi_a[k], i_magjs[index50[k]], compi_gmi[index50[k]][k], file=f)
    line_gmi = np.array(line_gmi)
    line_i = np.array(line_i)
    
    # print compi_gmi
    extent = [gmi_a[0], gmi_a[-1], i_magjs[-1], i_magjs[0]]
    
    plt.clf()
    plt.imshow(minuscomp, extent=extent, interpolation='nearest', cmap=cm.jet)
    plt.plot(line_gmi, line_i, color='white', linestyle='--')
    plt.xlabel('$(g-i)_0$')
    plt.ylabel('$i_0$')
    plt.colorbar()
    plt.savefig('compl_grid.pdf')
    
if __name__ == '__main__':
    main()

