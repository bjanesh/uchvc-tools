import numpy as np
import matplotlib.pyplot as plt
# XCEN YCEN RA DEC MAG MAGERR ELL PA FWHM ID FILTER

mag_r, magerr_r, fwhm_r = np.loadtxt('/Volumes/galileo/uchvc/targets/agc198606_copy/point_source_info', usecols=(4,5,8), unpack=True)
filterid = np.loadtxt('/Volumes/galileo/uchvc/targets/agc198606_copy/point_source_info', usecols=(10,), dtype=str, unpack=True)
keepi = [j for j,item in enumerate(filterid) if (item.endswith('i'))]

mag_i, magerr_i, fwhm_i = mag_r[keepi], magerr_r[keepi], fwhm_r[keepi]

leftcut = np.loadtxt('/Volumes/galileo/uchvc/targets/agc198606_copy/cutleft.txt', usecols=(0,), dtype=int, unpack=True)
fc = np.loadtxt('/Volumes/galileo/uchvc/targets/agc198606_copy/index_fc.txt', usecols=(0,), dtype=int, unpack=True)

mag, magerr, fwhm = mag_i[leftcut], magerr_i[leftcut], fwhm_i[leftcut]

plt.scatter(mag, fwhm, c='black')
plt.scatter(mag[fc], fwhm[fc], c='red')
plt.hlines([np.median(fwhm)], -12, 0, colors='red', linestyle='dashed')
plt.xlabel('$mag_{1x}$')
plt.ylabel('fwhm')
plt.xlim(-12,0)
plt.ylim(0,20)
plt.savefig('fwhmcheck_old.pdf')


# mag_r, fwhm_r = np.loadtxt('/Volumes/galileo/uchvc/targets/agc198606/escut_i.pos', usecols=(4,5,8), unpack=True)
# filterid = np.loadtxt('/Volumes/galileo/uchvc/targets/agc198606/escut_i.pos', usecols=(10,), dtype=str, unpack=True)
# keepi = [j for j,item in enumerate(filterid) if (item.endswith('i'))]
# 
# mag_i, magerr_i, fwhm_i = mag_r[keepi], magerr_r[keepi], fwhm_r[keepi]
# 
# leftcut = np.loadtxt('/Volumes/galileo/uchvc/targets/agc198606_copy/cutleft.txt', usecols=(0,), dtype=int, unpack=True)
# fc = np.loadtxt('/Volumes/galileo/uchvc/targets/agc198606_copy/index_fc.txt', usecols=(0,), dtype=int, unpack=True)
# 
# mag, magerr, fwhm = mag_i[leftcut], magerr_i[leftcut], fwhm_i[leftcut]
# 
# plt.scatter(mag, fwhm, c='black')
# plt.scatter(mag[fc], fwhm[fc], c='red')
# plt.hlines([np.median(fwhm)], -12, 0, colors='red', linestyle='dashed')
# plt.xlabel('$mag_{1x}$')
# plt.ylabel('fwhm')
# plt.xlim(-12,0)
# plt.ylim(0,20)
# plt.savefig('fwhmcheck_old.pdf')