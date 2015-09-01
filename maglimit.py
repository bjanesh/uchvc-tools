import os
import numpy as np
from pyraf import iraf
from matplotlib.path import Path

iraf.images(_doprint=0)
iraf.tv(_doprint=0)
iraf.ptools(_doprint=0)
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.photcal(_doprint=0)
iraf.apphot(_doprint=0)  
iraf.imutil(_doprint=0)

# specify the target and the file with the coordinate peak
title_string = os.getcwd()[-9:].upper() # get the target name from the folder name
coords_file = 'region_coords.dat.CHECK'
# title_string = input('Enter the root of the image name e.g. AGC123456: ')

# make sure you specify the distance to use for mass estimate, usually best distance
# dist = 385000 		# in pc
# dm = input('Enter the distance modulus to the source (m-M) : ') # ask for input
dm = 22.89
dist = pow(10,((dm + 5.)/5.)) # convert the dm to pc

# figure out what the aperture and sky region have to be
radius = ((((1.4*60)/206265)*420000)/dist)*206265 # in arcsec

rPixels = int(radius/0.11) # pODI pixel scale is 0.11"/px
print 'Aperture radius:', radius/60., rPixels, rPixels+36, rPixels+36+100

radiiAperture = repr(rPixels-3)+','+repr(rPixels-2)+','+repr(rPixels-1)+','+repr(rPixels)+','+repr(rPixels+1)+','+repr(rPixels+2)+','+repr(rPixels+3)

radiiSky = repr(rPixels+36-3)+','+repr(rPixels+36-2)+','+repr(rPixels+36-1)+','+repr(rPixels+36)+','+repr(rPixels+36+1)+','+repr(rPixels+36+2)+','+repr(rPixels+36+3)+','+repr(rPixels+36+100-3)+','+repr(rPixels+36+100-2)+','+repr(rPixels+36+100-1)+','+repr(rPixels+36+100)+','+repr(rPixels+36+100+1)+','+repr(rPixels+36+100+2)+','+repr(rPixels+36+100+3)

# if there isn't an i band masked image, guide the user through creating one
if not os.path.isfile(title_string+'_i_sh_masked.fits'):
    # load in the positions and magnitudes of calibrated stars
    ix,iy,imag,gmi = np.loadtxt('calibrated_mags.dat',usecols=(4,5,6,8),unpack=True)
    # make a list of bright stars
    with open('bright_stars.dat','w+') as f1:
        for i in range(len(ix)):
            if imag[i] < 18.0 :
                print >> f1, ix[i], iy[i], imag[i]
    # make a list of all stars
    with open('all_stars.dat','w+') as f1:
        for i in range(len(ix)):
            print >> f1, ix[i], iy[i], imag[i]
            
    with open('red_stars.dat','w+') as f1:
        for i in range(len(ix)):
            if gmi[i] > 2.0 :
                print >> f1, ix[i], iy[i], imag[i]
                
    with open('blue_stars.dat','w+') as f1:
        for i in range(len(ix)):
            if gmi[i] < 0.0 :
                print >> f1, ix[i], iy[i], imag[i]
    
    # # display the i band image
    # iraf.tv.display(image=title_string+'_i_sh.fits', frame=1)
    # # mark the positions from the lists above
    # iraf.unlearn(iraf.tv.tvmark)
    # iraf.tv.tvmark.setParam('label',"no")
    # iraf.tv.tvmark.setParam('pointsize',7)
    # iraf.tv.tvmark.setParam('mark',"circle")
    # # iraf.tv.tvmark(frame=1, coords='bright_stars.dat', radii="48,49,50,51,52,53", color=208)
    # iraf.tv.tvmark(frame=1, coords='all_stars.dat', radii="58,59,60,61,62,63", color=209)
    # # iraf.tv.tvmark(frame=1, coords='red_stars.dat', radii="64,65,66,67", color=204)
    # iraf.tv.tvmark(frame=1, coords='f_list_old_2.0_'+str(dm)+'_'+title_string+'.dat', radii="68,69,70,71,72,73", color=212)
    # iraf.tv.tvmark(frame=1, coords=coords_file, radii=radiiAperture, color=207)
    # iraf.tv.tvmark(frame=1, coords=coords_file, radii=radiiSky, color=205)
    
    # wait for the user to make a regions file from ds9
    while not os.path.isfile('regions.txt') :
        print 'Mask out bright stars indicated and other obvious things and save as regions.txt in IRAF/PROS'
        raw_input("Press Enter when finished:")
    
    # make a copy of the original image
    iraf.images.imcopy(title_string+'_i_sh.fits',title_string+'_i_sh_masked.fits',verbose="yes")
    
    # load in regions and calculate the mean sky background
    m3,m4,m5,m6 = np.loadtxt('regions.txt',usecols=(2,3,4,5),unpack=True)
    bgmean_i = 705.0328 #np.loadtxt('bgvals_i.txt',usecols=(2,),unpack=True)
    print "i sky value is:",bgmean_i
    
    # change regions into proper format (see kathy's .cl scripts) and replace regions with mean sky background
    # suggestion from John: don't fill with the sky value, instead set a very negative number and exclude these
    # pixels with a datamin parameter
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
        iraf.imutil.imreplace(images=title_string+"_i_sh_masked.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=0.0)

# do the same for the g-band image
if not os.path.isfile(title_string+'_g_sh_masked.fits'):
    iraf.images.imcopy(title_string+'_g_sh.fits',title_string+'_g_sh_masked.fits',verbose="yes")

    m3,m4,m5,m6 = np.loadtxt('regions.txt',usecols=(2,3,4,5),unpack=True)
    bgmean_g = 256.2528 #np.loadtxt('bgvals_g.txt',usecols=(2,),unpack=True)
    print "g sky value is:",bgmean_g

    # suggestion from John: don't fill with the sky value, instead set a very negative number and exclude these
    # pixels with a datamin parameter
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
        iraf.imutil.imreplace(images=title_string+"_g_sh_masked.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=0.0)

# if not os.path.isfile(title_string+'_i_sh_maskR.fits'):
#     iraf.images.imcopy(title_string+'_i_sh.fits',title_string+'_i_sh_maskR.fits',verbose="yes")
# 
#     m3,m4,m5,m6 = np.loadtxt('regionsRed.txt',usecols=(2,3,4,5),unpack=True)
#     bgmean_i = np.loadtxt('bgvals_i.txt',usecols=(2,),unpack=True)
#     print "i sky value is:",np.mean(bgmean_i)
# 
#     for i in range(len(m3)) :
#         x1 = m3[i] - (m5[i]/2.)
#         x2 = m3[i] + (m5[i]/2.)
#         y1 = m4[i] - (m6[i]/2.)
#         y2 = m4[i] + (m6[i]/2.)
# 
#         if (x1 < 0): 
#             x1 = 1
#         if (y1 < 0): 
#             y1 = 1
#         if (x2 > 11000):
#             x2 = 11000
#         if (y2 > 11000):
#             y2 = 11000
#         iraf.unlearn(iraf.imreplace)
#         iraf.imutil.imreplace(images=title_string+"_i_sh_maskR.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=np.mean(bgmean_i))
#         
# if not os.path.isfile(title_string+'_g_sh_maskR.fits'):
#     iraf.images.imcopy(title_string+'_g_sh.fits',title_string+'_g_sh_maskR.fits',verbose="yes")
# 
#     m3,m4,m5,m6 = np.loadtxt('regionsRed.txt',usecols=(2,3,4,5),unpack=True)
#     bgmean_g = np.loadtxt('bgvals_g.txt',usecols=(2,),unpack=True)
#     print "g sky value is:",np.mean(bgmean_g)
# 
#     for i in range(len(m3)) :
#         x1 = m3[i] - (m5[i]/2.)
#         x2 = m3[i] + (m5[i]/2.)
#         y1 = m4[i] - (m6[i]/2.)
#         y2 = m4[i] + (m6[i]/2.)
# 
#         if (x1 < 0): 
#             x1 = 1
#         if (y1 < 0): 
#             y1 = 1
#         if (x2 > 11000):
#             x2 = 11000
#         if (y2 > 11000):
#             y2 = 11000
#         iraf.unlearn(iraf.imreplace)
#         iraf.imutil.imreplace(images=title_string+"_g_sh_maskR.fits["+repr(int(x1))+":"+repr(int(x2))+","+repr(int(y1))+":"+repr(int(y2))+"]", value=np.mean(bgmean_g))
# 
# # display the masked image
# iraf.tv.display(image=title_string+'_g_sh_masked.fits', frame=1)
# iraf.tv.display(image=title_string+'_i_sh_masked.fits', frame=2)
# # mark the phot and sky regions
# iraf.unlearn(iraf.tv.tvmark)
# iraf.tv.tvmark.setParam('label',"no")
# iraf.tv.tvmark.setParam('pointsize',7)
# iraf.tv.tvmark.setParam('mark',"circle")
# iraf.tv.tvmark(frame=1, coords=coords_file, radii=radiiAperture, color=207)
# iraf.tv.tvmark(frame=1, coords=coords_file, radii=radiiSky, color=205)
# iraf.tv.tvmark(frame=2, coords=coords_file, radii=radiiAperture, color=207)
# iraf.tv.tvmark(frame=2, coords=coords_file, radii=radiiSky, color=205)

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
iraf.apphot.phot(image=title_string+"_i_shPhotCutout.fits", coords=coords_file, output="mag_est_i.dat") 

# txdump to a file for easier data handling
txdump_out = open('phot_region_i.txdump','w+')
iraf.ptools.txdump(textfiles='mag_est_i.dat', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

# iraf.datapars.setParam('fwhmpsf',7.045) # don't need to set this

# do the actual phot
iraf.apphot.phot(image=title_string+"_g_shPhotCutout.fits", coords=coords_file, output="mag_est_g.dat") 

# txdump to a file for easier data handling
txdump_out = open('phot_region_g.txdump','w+')
iraf.ptools.txdump(textfiles='mag_est_g.dat', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
txdump_out.close()

# # do the actual phot
# iraf.apphot.phot(image=title_string+"_i_sh_maskR.fits", coords=coords_file, output="mag_est_iR.dat") 
# 
# # txdump to a file for easier data handling
# txdump_out = open('phot_region_iR.txdump','w+')
# iraf.ptools.txdump(textfiles='mag_est_iR.dat', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
# txdump_out.close()
# 
# # iraf.datapars.setParam('fwhmpsf',7.045) # don't need to set this
# 
# # do the actual phot
# iraf.apphot.phot(image=title_string+"_g_sh_maskR.fits", coords=coords_file, output="mag_est_gR.dat") 
# 
# # txdump to a file for easier data handling
# txdump_out = open('phot_region_gR.txdump','w+')
# iraf.ptools.txdump(textfiles='mag_est_gR.dat', fields="id,mag,merr,msky,stdev,rapert,xcen,ycen,ifilter,xairmass,image", expr='yes', headers='no', Stdout=txdump_out)
# txdump_out.close()

# get the zero points from steven's files
# no i-band color term, use line 2
with open(title_string+'_i_phot.zp') as f_zp_i :
	data = f_zp_i.read()
	fl_i = data.split('\n', 1)[1]
	zp_vals_i = fl_i.split()
	cal_zp_i = float(zp_vals_i[0])
	cal_zp_ie = float(zp_vals_i[1])

# g-band uses the color term on line 1
with open(title_string+'_g_phot.zp') as f_zp_g :	
	data_g = f_zp_g.read()
	fl_g = data_g.split('\n',1)[0]
	zp_vals_g = fl_g.split()
	cal_zp_g = float(zp_vals_g[2])
	cal_color_g = float(zp_vals_g[0])
	cal_zp_ge= float(zp_vals_g[3])
	cal_color_ge = float(zp_vals_g[1])

# read in phot magnitudes from the txdumps
mags_iB = np.loadtxt('phot_region_i.txdump',usecols=(1,),unpack=True)
mags_gB = np.loadtxt('phot_region_g.txdump',usecols=(1,),unpack=True)
# mags_iR = np.loadtxt('phot_region_iR.txdump',usecols=(1,),unpack=True)
# mags_gR = np.loadtxt('phot_region_gR.txdump',usecols=(1,),unpack=True)

# get the positions and magnitudes for all the filter stars 
ixr, iyr, mags_ifr, mags_gfr = np.loadtxt('f_list_old_2.0_'+str(dm)+'_'+title_string+'.dat',usecols=(0,1,4,5),unpack=True)
xy_points = zip(ixr,iyr)
# get the position of the circle center
ccxr, ccyr = np.loadtxt('region_coords.dat',usecols=(0,1),unpack=True)
# make a testing circle
cosd = lambda x : np.cos(np.deg2rad(x))
sind = lambda x : np.sin(np.deg2rad(x))
x_circr = [ccxr + rPixels*cosd(t) for t in range(0,359,1)]
y_circr = [ccyr + rPixels*sind(t) for t in range(0,359,1)]
verts_circr = zip(x_circr,y_circr)
rcirc_filter = Path(verts_circr)
# test to see if the points are in the circle
stars_circr = rcirc_filter.contains_points(xy_points)

mags_if = np.array([mags_ifr[i] for i in range(len(mags_ifr)) if (stars_circr[i])])
mags_gf = np.array([mags_gfr[i] for i in range(len(mags_gfr)) if (stars_circr[i])])

print len(mags_if),'filter stars in the aperture'
# absolute magnitudes
iabsf = mags_if - dm
gabsf = mags_gf - dm

iAbsSolar = 4.57	# sloan i solar absolute magnitude from sparke & gallagher 2e
gAbsSolar =  5.12 # sloan g solar absolute magnitude from sparke & gallagher 2e

# find the intrinsic flux for each star in the filter
flux_i = np.power(10,((iabsf-iAbsSolar)/-2.5))
flux_g = np.power(10,((gabsf-gAbsSolar)/-2.5))

# find the total luminosity
fluxTi = np.sum(flux_i)
fluxTg = np.sum(flux_g)

# and finally the absolute magnitude again
mags_iC = iAbsSolar - 2.5*np.log10(fluxTi)
mags_gC = gAbsSolar - 2.5*np.log10(fluxTg)

print mags_iC,mags_gC-mags_iC 

mags_i = np.array([mags_iB])
mags_g = np.array([mags_gB])

# calculate calibrated magnitudes
# -0.7525749893 is doubling the flux, since we said 2 arcmin is the half-light radius
# this is a BIG assumption

mu_gi = 1.055
zp_gi = 0.571
eps_gi = 0.007
zp_i = 25.894

kg = 0.200
ki = 0.058

amg = 1.082642493
ami = 1.208449087

cal_A_g = 0.1023
cal_A_i = 0.0526

# g-i = mu_gi * (g0 - i0) + ZP_gi
# i = eps_gi * (g-i) + ZP_i
g0 = mags_g - (kg*amg)
i0 = mags_i - (ki*ami)
gmi = mu_gi*(g0-i0) + zp_gi

i_mag = i0 + eps_gi*gmi + zp_i #- cal_A_i 
g_mag = gmi + i_mag - cal_A_g - 0.7525749893
i_mag = i_mag - cal_A_i - 0.7525749893

i_mag = np.append(i_mag,mags_iC+dm)
g_mag = np.append(g_mag,mags_gC+dm)
gmi = g_mag - i_mag

# estimate an absolute magnitude based on the input distance
absMag = i_mag - dm # in i-band magnitudes

# estimate a luminosity from the absolute magnitude (M = M_sun,i - 2.5log10(L/Lsun))
lum_i = np.power(10,(absMag - iAbsSolar)/-2.5) # in solar luminosities

# get a stellar mass estimate using g-i color and the formula from Bell+2003 (ApJSS, 149, 289)
# note that we are using the i magnitude and g-i color, +  the corresponding coefficients a_i and b_i from Table 7 of Bell+2003
logml = -0.152 + 0.518*(g_mag-i_mag) # first log(M/L)
mtol = np.power(10,logml)				  # get the actual number (10^x)
stellarMass = mtol * lum_i				  # mass = M/L * L

for i in range(len(i_mag)):
	print '# i mag, g mag, (g-i)'
	print i_mag[i],g_mag[i],g_mag[i]-i_mag[i]
	print '# absolute i mag, L, M/L, M'
	print absMag[i],lum_i[i],mtol[i],stellarMass[i]