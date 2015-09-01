import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from scipy.stats import lognorm

bins = 165
width = 22 
fwhm = 2.0

sig = ((bins/width)*fwhm)/2.355

alphas = []
betas = []
locs = []
ns = []
for n in range(100,450,10) :
	print n
	ns.append(n)
# n = 123
# with open('valuesLeoP.txt','w+') as f1:
	sig_values_r = []
	for i in range(25000) :
		random_ra = 20.0*np.random.random_sample((n,))
		random_dec = 20.0*np.random.random_sample((n,))
		random_xy = zip(random_ra,random_dec)
		grid_r, xedges_r, yedges_r = np.histogram2d(random_dec, random_ra, bins=[bins,bins], range=[[0,width],[0,width]])
		hist_points_r = zip(xedges_r,yedges_r)
		grid_gaus_r = ndimage.filters.gaussian_filter(grid_r, sig, mode='constant', cval=0)
		S_r = np.array(grid_gaus_r*0)
		
		grid_mean_r = np.mean(grid_gaus_r)
		grid_sigma_r = np.std(grid_gaus_r)
		S_r = (grid_gaus_r-grid_mean_r)/grid_sigma_r
		
		x_cent_r, y_cent_r = np.unravel_index(grid_gaus_r.argmax(),grid_gaus_r.shape)
		sig_values_r.append(S_r[x_cent_r][y_cent_r])
		# print >> f1, S_r[x_cent_r][y_cent_r]
	al,loc,beta=lognorm.fit(sig_values_r)
	alphas.append(al)
	betas.append(beta)
	locs.append(loc)
	# pct_calc = [sig_values_r[i] for i in range(len(sig_values_r)) if (sig_values_r[i] < S_th)]
	# percentile = (float(len(pct_calc))/1000.0)*100.0
	# print n, S_th, percentile

ax0 = plt.subplot(2,2,1)	
plt.scatter(ns,alphas,c='r', edgecolors='none')
# plt.ylim(0,1.1)
# plt.xlim(2,12)
plt.xlabel('sample size')
plt.ylabel('alpha')

ax1 = plt.subplot(2,2,2)	
plt.scatter(ns,betas,c='g', edgecolors='none')
# plt.ylim(0,1.1)
# plt.xlim(2,12)
plt.xlabel('sample size')
plt.ylabel('beta')

ax2 = plt.subplot(2,2,3)	
plt.scatter(ns,locs,c='b', edgecolors='none')
# plt.ylim(0,1.1)
# plt.xlim(2,12)
plt.xlabel('sample size')
plt.ylabel('mu')

plt.savefig('fitstuff.pdf')
