from __future__ import print_function
import os, sys
import numpy as np
from pyraf import iraf
from odi_calibrate import download_sdss, js_calibrate
from sdss_fit import getVabs
from collections import OrderedDict

iraf.images(_doprint=0)
iraf.tv(_doprint=0)
iraf.ptools(_doprint=0)
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.photcal(_doprint=0)
iraf.apphot(_doprint=0)
iraf.imutil(_doprint=0)

def getHImass(object, dm):
    # print object, mpc
    uchvcdb = os.path.dirname(os.path.abspath(__file__))+'/predblist.sort.csv'
    name, mass = np.loadtxt(uchvcdb, usecols=(1,6), dtype=str, delimiter=',', unpack=True)
    # find the right row
    coord = [i for i,this in enumerate(mass) if object.upper() in name[i]][0]
    # print 'the HI mass of', name[coord], 'is', mass[coord], 'at 1 Mpc'
    
    # mpc = mpc/1000.
    mpc = pow(10,((dm + 5.)/5.))/1000000.
    logm = float(mass[coord])
    mass = mpc*mpc*10**logm  # make sure to scale by the distance in Mpc^2
    
    print('{:3.1f}'.format(np.log10(mass)))
    return mass

def main():
    epadu = 1.268899
    path = os.getcwd()
    steps = path.split('/')
    title_string = steps[-1].upper()        # which should always exist in the directory
    coords_file = 'region_coords.dat'
    # dm = 26.07
    print("computing magnitude estimates for", title_string)
    # dm = float(raw_input("Enter the distance modulus: "))
    dm = float(sys.argv[1])
    print("at a distance modulus of", dm)

    if not os.path.isfile(title_string+'_i_masked.fits'):

    #   ix,iy,imag = np.loadtxt('calibrated_mags.dat',usecols=(4,5,6),unpack=True)
    #     with open('bright_stars.dat','w+') as f1:
    #         for i in range(len(ix)):
    #             if imag[i] < 18.0 :
    #                 print >> f1, ix[i], iy[i], imag[i]

        # iraf.tv.display(image=title_string+'_i.fits', frame=1)
        #
        # iraf.unlearn(iraf.tv.tvmark)
        # iraf.tv.tvmark.setParam('label',"no")
        # iraf.tv.tvmark.setParam('pointsize',7)
        # iraf.tv.tvmark.setParam('mark',"circle")
        # iraf.tv.tvmark(frame=1, coords='bright_stars.dat', radii="98,99,100,101,102,103", color=208)

        while not os.path.isfile('regions.txt') :
            print('Mask out bright stars indicated and other obvious things and save as regions.txt')
            input("Press Enter when finished:")

        iraf.images.imcopy(title_string+'_i.fits',title_string+'_i_masked.fits',verbose="yes")

        m3,m4,m5,m6 = np.loadtxt('regions.txt',usecols=(2,3,4,5),unpack=True)
        # bgmean_i = 623.670898438
        # print "i sky value is:",np.mean(bgmean_i)
        #
        # bgmean_g = 191.777572632
        # print "g sky value is:",np.mean(bgmean_g)

        for i in range(len(m3)) :
            x1 = m3[i] - (m5[i]/2.)
            x2 = m3[i] + (m5[i]/2.)
            y1 = m4[i] - (m6[i]/2.)
            y2 = m4[i] + (m6[i]/2.)

            if (x1 < 0):
                x1 = 1
            if (y1 < 0):
                y1 = 1
            if (x2 > 11000):
                x2 = 11000
            if (y2 > 11000):
                y2 = 11000
            iraf.unlearn(iraf.imreplace)
            iraf.imutil.imreplace(images=title_string+"_i_masked.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=0.0)

    if not os.path.isfile(title_string+'_g_masked.fits'):

        ix,iy,imag = np.loadtxt('calibrated_mags.dat',usecols=(4,5,6),unpack=True)
        with open('bright_stars.dat','w+') as f1:
            for i in range(len(ix)):
                if imag[i] < 18.0 :
                    print(ix[i], iy[i], imag[i], file=f1)

        iraf.images.imcopy(title_string+'_g.fits',title_string+'_g_masked.fits',verbose="yes")

        m3,m4,m5,m6 = np.loadtxt('regions.txt',usecols=(2,3,4,5),unpack=True)

        # print "g sky value is:",np.mean(bgmean_g)

        for i in range(len(m3)) :
            x1 = m3[i] - (m5[i]/2.)
            x2 = m3[i] + (m5[i]/2.)
            y1 = m4[i] - (m6[i]/2.)
            y2 = m4[i] + (m6[i]/2.)

            if (x1 < 0):
                x1 = 1
            if (y1 < 0):
                y1 = 1
            if (x2 > 11000):
                x2 = 11000
            if (y2 > 11000):
                y2 = 11000
            iraf.unlearn(iraf.imreplace)
            iraf.imutil.imreplace(images=title_string+"_g_masked.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=0.0)

    if not os.path.isfile('ones_mask.fits'):
        iraf.images.imarith(title_string+'_g.fits', '*', 0.0, 'zeros.fits',verbose="yes")
        iraf.images.imarith('zeros.fits', '+', 1.0, 'ones_mask.fits',verbose="yes")

        m3,m4,m5,m6 = np.loadtxt('regions.txt',usecols=(2,3,4,5),unpack=True)

        for i in range(len(m3)) :
            x1 = m3[i] - (m5[i]/2.)
            x2 = m3[i] + (m5[i]/2.)
            y1 = m4[i] - (m6[i]/2.)
            y2 = m4[i] + (m6[i]/2.)

            if (x1 < 0):
                x1 = 1
            if (y1 < 0):
                y1 = 1
            if (x2 > 11000):
                x2 = 11000
            if (y2 > 11000):
                y2 = 11000
            iraf.unlearn(iraf.imreplace)
            iraf.imutil.imreplace(images="ones_mask.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=0.0)
    # iraf.tv.display(image=title_string+'_i_masked.fits', frame=1)
    #
    # iraf.unlearn(iraf.tv.tvmark)
    # iraf.tv.tvmark.setParam('label',"no")
    # iraf.tv.tvmark.setParam('pointsize',7)
    # iraf.tv.tvmark.setParam('mark',"circle")
    # iraf.tv.tvmark(frame=1, coords=coords_file, radii="1633,1634,1635,1636,1637,1638,1639,1697,1698,1699,1700,1701,1702,1703,1797,1798,1799,1800,1801,1802,1803", color=208)
    # 
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
    iraf.fitskypars.setParam('dannulus',100.)
    iraf.datapars.setParam('fwhmpsf',6.197) 
    iraf.photpars.setParam('apertures','409,500,591,682,773,818') 
    iraf.fitskypars.setParam('annulus',450.)
    #
    iraf.apphot.phot(image=title_string+"_i_masked.fits", coords=coords_file, output="mag_est_i.dat")

    txdump_out = open('phot_region_i.txdump','w+')
    iraf.ptools.txdump(textfiles='mag_est_i.dat', fields="id,sum,msky,stdev,nsky", expr='yes', headers='no', Stdout=txdump_out)
    txdump_out.close()

    iraf.datapars.setParam('fwhmpsf',6.197)

    iraf.apphot.phot(image=title_string+"_g_masked.fits", coords=coords_file, output="mag_est_g.dat")

    txdump_out = open('phot_region_g.txdump','w+')
    iraf.ptools.txdump(textfiles='mag_est_g.dat', fields="id,sum,msky,stdev,nsky", expr='yes', headers='no', Stdout=txdump_out)
    txdump_out.close()

    iraf.apphot.phot(image="ones_mask.fits", coords=coords_file, output="mag_area.dat")
    txdump_out = open('area.txdump','w+')
    iraf.ptools.txdump(textfiles='mag_area.dat', fields="id,sum", expr='yes', headers='no', Stdout=txdump_out)
    txdump_out.close()

    # calculate the magnitude cf. apphot.phot, but *manually* subtract the area determined by the ones mask
    areas = np.loadtxt('area.txdump', usecols=(1,2,3,4,5,6), unpack=True)
    sum_i = np.loadtxt('phot_region_i.txdump',usecols=(1,2,3,4,5,6),unpack=True)
    sky_i, stdev_i, nsky_i = np.loadtxt('phot_region_i.txdump',usecols=(7,8,9),unpack=True)
    sum_g = np.loadtxt('phot_region_g.txdump',usecols=(1,2,3,4,5,6),unpack=True)
    sky_g, stdev_g, nsky_g = np.loadtxt('phot_region_g.txdump',usecols=(7,8,9),unpack=True)

    fl_i = sum_i - areas * sky_i
    mag_i = -2.5*np.log10(fl_i) + 2.5*np.log10(300.)
    error_i = np.sqrt(fl_i/ epadu + areas * stdev_i**2 + areas**2 * stdev_i**2 / nsky_i)
    merr_i = 1.0857 * error_i / fl_i

    fl_g = sum_g - areas * sky_g
    mag_g = -2.5*np.log10(fl_g) + 2.5*np.log10(300.)
    error_g = np.sqrt(fl_g/ epadu + areas * stdev_g**2 + areas**2 * stdev_g**2 / nsky_g)
    merr_g = 1.0857 * error_g / fl_g

    print(fl_i,mag_i,merr_i)
    print(fl_g,mag_g,merr_g)

    flux_g, flux_i, merrs_g, merrs_i = [], [], [], []
    rs = np.array([51, 77, 90, 180]) # 2' cell, 3' cell, 3' diam, 3' radius
    for r in rs:
        fcirc_file = 'circle'+repr(r)+'.txt'

        iraf.fitskypars.setParam('dannulus',10.)
        iraf.datapars.setParam('fwhmpsf',6.197) 
        iraf.photpars.setParam('apertures',7.) 
        iraf.fitskypars.setParam('annulus',10.)
        #
        iraf.apphot.phot(image=title_string+"_g.fits", coords=fcirc_file, output="mag_min_g_{:2d}.dat".format(r))
        txdump_out = open('phot_indiv_g.txdump','w+')
        iraf.ptools.txdump(textfiles="mag_min_g_{:2d}.dat".format(r), fields="id,mag,merr,flux,area,stdev,nsky", expr='yes', headers='no', Stdout=txdump_out)
        txdump_out.close()

        iraf.fitskypars.setParam('dannulus',10.)
        iraf.datapars.setParam('fwhmpsf',6.197) 
        iraf.photpars.setParam('apertures',7.) 
        iraf.fitskypars.setParam('annulus',10.)
        #
        iraf.apphot.phot(image=title_string+"_i.fits", coords=fcirc_file, output="mag_min_i_{:2d}.dat".format(r))
        txdump_out = open('phot_indiv_i.txdump','w+')
        iraf.ptools.txdump(textfiles="mag_min_i_{:2d}.dat".format(r), fields="id,mag,merr,flux,area,stdev,nsky", expr='yes', headers='no', Stdout=txdump_out)
        txdump_out.close()

        fluxes_i, areas_i, stdevs_i, nskys_i = np.loadtxt('phot_indiv_i.txdump',usecols=(3,4,5,6),unpack=True)
        fluxes_g, areas_g, stdevs_g, nskys_g = np.loadtxt('phot_indiv_g.txdump',usecols=(3,4,5,6),unpack=True)
        flux_i.append(np.sum(fluxes_i))
        flux_g.append(np.sum(fluxes_g))

        area_i = np.sum(areas_i)
        stdev_i = np.median(stdevs_i)
        nsky_i = np.sum(nskys_i)

        error_i = np.sqrt(np.sum(fluxes_i) / epadu + area_i * stdev_i**2 + area_i**2 * stdev_i**2 / nsky_i)
        merrs_i.append(1.0857 * error_i / np.sum(fluxes_i))

        area_g = np.sum(areas_g)
        stdev_g = np.median(stdevs_g)
        nsky_g = np.sum(nskys_g)

        error_g = np.sqrt(np.sum(fluxes_g) / epadu + area_g * stdev_g**2 + area_g**2 * stdev_g**2 / nsky_g)
        merrs_g.append(1.0857 * error_g / np.sum(fluxes_g))
    
        # os.remove("mag_min_g.dat")
        # os.remove("mag_min_i.dat")
        os.remove("phot_indiv_g.txdump")
        os.remove("phot_indiv_i.txdump")

    # print fl_i, flux_i, merr_i
    # print fl_g, flux_g, merr_g

    mags_i = np.hstack((-2.5*np.log10(fl_i)+2.5*np.log10(300.0), -2.5*np.log10(np.array(flux_i))+2.5*np.log10(300.0)))
    mags_g = np.hstack((-2.5*np.log10(fl_g)+2.5*np.log10(300.0), -2.5*np.log10(np.array(flux_g))+2.5*np.log10(300.0)))
    me_i = np.hstack((merr_i, merrs_i))
    me_g = np.hstack((merr_g, merrs_g))
    # mags_i = -2.5*np.log10(flux_i)+2.5*np.log10(300.0)
    # mags_g = -2.5*np.log10(flux_g)+2.5*np.log10(300.0)

    # print mags_i, mags_g
    download_sdss(title_string+"_g.fits", title_string+"_i.fits", gmaglim = 21)
    eps_g, std_eps_g, zp_g, std_zp_g, eps_i, std_eps_i, zp_i, std_zp_i = js_calibrate(img1 = title_string+"_g.fits", img2 = title_string+"_i.fits", verbose=False)

    # values determined by ralf/daniel @ wiyn
    kg = 0.20
    kr = 0.12
    ki = 0.058

    # get the photometric calibration coefficients from Steven's help file <--
    # or from the image header/fits table/ whatever
    photcalFile = open(title_string+'_help_js.txt')
    photcal = photcalFile.read()
    photcalLines = photcal.splitlines()

    amg = float(photcalLines[25].split()[5])
    ami = float(photcalLines[26].split()[5])
    photcalFile.close()

    print(amg, ami)

    if not os.path.isfile('extinction.tbl.txt'):
        print('Fetching extinction table for',fits_h_i[0].header['RA'],fits_h_i[0].header['DEC'])
        getexttbl(fits_h_i[0].header['RA'],fits_h_i[0].header['DEC'])

    LamEff,A_over_E_B_V_SandF,A_SandF,A_over_E_B_V_SFD,A_SFD= np.genfromtxt('extinction.tbl.txt', usecols=(2,3,4,5,6),unpack=True,skip_header=27,skip_footer=12)
    A_id = np.genfromtxt('extinction.tbl.txt', usecols=(1,),dtype=str,unpack=True,skip_header=27,skip_footer=12)
    E_B_V = np.genfromtxt('extinction.tbl.txt', usecols=(2,),skip_header=1,skip_footer=42)

    for j in range(len(A_id)):
        if A_id[j] == 'g':
            cal_A_g = A_over_E_B_V_SandF[j]*0.86*E_B_V # E(B-V) is the Schlegel+ value, S&F say with their calibration
    for j in range(len(A_id)):                                  # use 0.86*E(B-V) instead. cf. S&F2011 pg 1, 2011ApJ...737..103S
        if A_id[j] == 'i':
            cal_A_i = A_over_E_B_V_SandF[j]*0.86*E_B_V

    print('Reddening correction :: g = {0:7.4f} : i = {1:7.4f}'.format(cal_A_g,cal_A_i))

    tolerance = 0.0001

    g_0 = mags_g - kg*amg
    i_0 = mags_i - ki*ami

    i_sun = 4.58
    m_hi = getHImass(title_string, dm)

    v_magr, g_magr, i_magr = np.loadtxt(os.path.dirname(os.path.abspath(__file__))+'/sdssBVR.dat', usecols=(4, 12, 16), unpack=True)
    good_g, good_i = np.loadtxt(os.path.dirname(os.path.abspath(__file__))+'sdssBVR.dat', usecols=(21,23), dtype=bool, unpack=True)

    good = np.where((v_magr-g_magr < 0.5) & (v_magr-g_magr > -1.5) & good_g & good_i)
    print(good[0].size, v_magr.size)
    v, g, i = v_magr[good], g_magr[good], i_magr[good]

    p = np.polyfit(g-i, v-g, 3)
    print(p)

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
    print(rms, chi_sq)

    # plt.scatter(g-i,v-g, edgecolors='none')
    # plt.plot(xplt, yplt, c='red')
    # plt.xlabel('g-i')
    # plt.ylabel('v-g')
    # plt.savefig('gi_to_v.pdf')
    # f = open('mag_est.txt', 'w+')
    # print '# r g     ge   i     ie   g-i  err   Mg    Mi   MV    M/L  L*      M*       HI/*'
    # print >> f, '# r g     ge   i     ie   g-i  err   Mg    Mi   MV    M/L  L*      M*       HI/*'

    rs = np.array([51, 77, 90, 180, 51, 77, 90, 180])

    with open('optical_props.txt', 'w+') as opt:
        print('# ap     g   ge     i   ie  g-i Eg-i    Mg    Mi    MV  M/L L*      MHI   M*  Hi/*')
        print('# ap     g   ge     i   ie  g-i Eg-i    Mg    Mi    MV  M/L L*      MHI   M*  Hi/*', file=opt)
        for i,r in enumerate(rs):
            color_guess = 0.0
            color_diff = 1.0
            while abs(color_diff) > tolerance:
                g_cal = g_0[i] + eps_g*color_guess + zp_g
                i_cal = i_0[i] + eps_i*color_guess + zp_i
    
                color_new = g_cal - i_cal
                color_diff = color_guess-color_new
                color_guess = color_new
                # print g_0[i], g_cal, i_0[i], i_cal, color_new
    
            g_mag = g_cal - cal_A_g
            i_mag = i_cal - cal_A_i
            gmi = g_mag - i_mag
            e_gmi = np.sqrt(me_g[i]**2 + me_i[i]**2)
            g_abs = g_mag-dm
            i_abs = i_mag-dm
        
            v_abs = getVabs(g_mag, i_mag, dm)
        
            mtol = np.power(10,0.518*gmi-0.152)
            l_star = np.power(10,(i_sun-i_abs)/2.5)
            m_star = l_star*mtol
            hitostar = m_hi/m_star
            print(' {:3d} {:5.2f} {:4.2f} {:5.2f} {:4.2f} {:4.2f} {:4.2f} {:5.2f} {:5.2f} {:5.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f} {:5.1f}'.format(r,g_mag,me_g[i],i_mag,me_i[i],g_mag-i_mag,e_gmi,g_abs,i_abs,v_abs,mtol,np.log10(l_star),np.log10(m_hi),np.log10(m_star),hitostar))
            print(' {:3d} {:5.2f} {:4.2f} {:5.2f} {:4.2f} {:4.2f} {:4.2f} {:5.2f} {:5.2f} {:5.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f} {:5.1f}'.format(r,g_mag,me_g[i],i_mag,me_i[i],g_mag-i_mag,e_gmi,g_abs,i_abs,v_abs,mtol,np.log10(l_star),np.log10(m_hi),np.log10(m_star),hitostar), file=opt)
    # 
    for file_ in ["area.txdump","mag_area.dat","mag_est_g.dat","phot_region_g.txdump","mag_est_i.dat","phot_region_i.txdump"]:
        os.remove(file_)
    

if __name__ == '__main__':
    main()
    # objects = OrderedDict([('AGC174540', 1100),
    # ('AGC198511',  380),
    # ('AGC198606',  880),
    # ('HI1037+21',  520),
    # ('AGC215417',  350),
    # ('HI1151+20',  900),
    # ('AGC226067',  880),
    # ('AGC227987',  820),
    # ('AGC229326',  250),
    # ('AGC238626',  450),
    # ('AGC238713', 1740),
    # ('AGC249000', 1900),
    # ('AGC249282',  630),
    # ('AGC249320', 1150),
    # ('AGC249323', 1900),
    # ('AGC249525', 2200),
    # ('AGC258237', 2400),
    # ('AGC258242', 1000),
    # ('AGC258459', 1340),
    # ('AGC268069',  700),
    # ('AGC268074',  260),
    # ('HI0959+19',  410),
    # ('HI1050+24',  270)])
    # for i,obj in enumerate(objects.keys()):
    #     getHImass(obj, objects[obj])