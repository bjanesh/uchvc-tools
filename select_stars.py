from matplotlib.path import Path
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS

iso_filename = 'iso_t+12_z-1.7.dat'
photcalib_filename = 'phot_calib_final'
wcs_source_image_filename = 'LeoII_i_trim.fits'
output_cmd_pdf = 'testcmd.pdf'

width = 0.2 # in mags
# dm = 21.69

###################
g_m_iso, i_m_iso = np.loadtxt(iso_filename, usecols=(0,1), unpack=True)
g_ierr, ix, iy, i_ierr, g_mag,i_mag,gmi = np.loadtxt(photcalib_filename, usecols=(5,6,7,10,11,12,13), unpack=True)
gmi_err = np.sqrt(g_ierr**2 + i_ierr**2)

# if you get an error about loading the WCS, uncomment the following lines to
# delete the pipeline WCS keywords from the header
# from pyraf import iraf
# iraf.imutil.hedit(images=wcs_source_image_filename, fields='PV*', delete='yes', verify='no')
w = WCS(wcs_source_image_filename)
ra, dec = w.all_pix2world(ix, iy, 1)

# i_m_iso = i_iso + dm
gi_iso = g_m_iso - i_m_iso

colors_left = gi_iso - width/2.0
colors_right = gi_iso + width/2.0

colors = np.concatenate((colors_left, np.flipud(colors_right)))
mags = np.concatenate((i_m_iso, np.flipud(i_m_iso)))

verts = zip(colors, mags)		# set up the Path necessary for testing membership
cm_filter = Path(verts)

stars_f = np.empty_like(gmi, dtype=bool)
for i in range(len(gmi)) : 
    nsteps_color = int(abs((gmi_err[i])/0.001))
    nsteps_mag = int(abs((i_ierr[i])/0.001))
    # print i, nsteps_color, nsteps_mag
    
    if nsteps_color == 0 :
        nsteps_color = 1
    if nsteps_mag == 0 :
        nsteps_mag = 1
    
    # point smearing - use the 1sigma error bars to test membership
    
    cm_points_l = [(gmi[i]-0.01*j*gmi_err[i],i_mag[i]) for j in range(nsteps_color)]
    cm_points_r = [(gmi[i]+0.01*j*gmi_err[i],i_mag[i]) for j in range(nsteps_color)]
    cm_points_u = [(gmi[i],i_mag[i]-0.01*j*i_ierr[i]) for j in range(nsteps_mag)]
    cm_points_d = [(gmi[i],i_mag[i]+0.01*j*i_ierr[i]) for j in range(nsteps_mag)]
            
    stars_f_l = cm_filter.contains_points(cm_points_l)		
    stars_f_r = cm_filter.contains_points(cm_points_r)		
    stars_f_u = cm_filter.contains_points(cm_points_u)		
    stars_f_d = cm_filter.contains_points(cm_points_d)		
    
    stars_f[i] = any(stars_f_l) | any(stars_f_r) | any(stars_f_u) | any(stars_f_d)  

# print stars_f
check = [stars_f[i] for i in range(len(stars_f)) if (stars_f[i])]
print len(check), "stars in filter"

iplt, giplt = i_mag[np.where(stars_f)], gmi[np.where(stars_f)]

plt.scatter(gmi, i_mag, c='black', s=2, edgecolors='none')
plt.scatter(giplt, iplt, c='red', s=5, edgecolors='none')
plt.plot(colors, mags)
plt.xlim(-2,5)
plt.ylim(25,15)
plt.savefig(output_cmd_pdf)

print 'printing!'
with open('members.txt', 'w+') as f:
    for i,chk in enumerate(stars_f):
        if chk:
            # ra, dec = w.all_pix2world(ix[i], iy[i], 1)
            print >> f, i_mag[i], ra[i], dec[i]