def distfit(n,dists,title):
	import numpy as np
	import matplotlib.pyplot as plt
	# from scipy.optimize import curve_fit
	from scipy.stats import lognorm
	from scipy import ndimage
	
	# n = 279
	bins = 165
	width = 22 
	fwhm = 2.0
	sig = ((bins/width)*fwhm)/2.355
	valsLP = []
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
		valsLP.append(S_r[x_cent_r][y_cent_r])
	# valsLP = np.loadtxt('valuesLeoP.txt', usecols=(0,), unpack=True)
	# vals = np.loadtxt('values.txt', usecols=(0,), unpack=True)
	
	# bins, edges = np.histogram(vals, bins=400, range=[2,22], normed=True)
	# centers = (edges[:-1] + edges[1:])/2.
	# plt.scatter(centers, bins, edgecolors='none')
	
	x = np.linspace(2, 22, 4000)
	
	# al,loc,beta=lognorm.fit(vals)
	# print al, loc, beta
	# # plt.plot(x, lognorm.pdf(x, al, loc=loc, scale=beta),'r-', lw=5, alpha=0.6, label='lognormal AGC198606')
	# print lognorm.cdf(dists, al, loc=loc, scale=beta)
	
	bins, edges = np.histogram(valsLP, bins=400, range=[2,22], normed=True)
	centers = (edges[:-1] + edges[1:])/2.
	
	
	# x = np.linspace(2, 22, 4000)
	# dists = np.array([3.958,3.685,3.897,3.317])
	al,loc,beta=lognorm.fit(valsLP)
	# print al, loc, beta
	plt.plot(x, lognorm.pdf(x, al, loc=loc, scale=beta),'r-', lw=2, alpha=0.6, label='lognormal distribution')
	print 'Significance of detection:','{0:6.3f}%'.format(100.0*lognorm.cdf(dists, al, loc=loc, scale=beta))
	
	plt.scatter(centers, bins, edgecolors='none', label='histogram of $\sigma$ from 25000 \nuniform random samples')
	# print chisqg(bins, lognorm.pdf(centers, al, loc=loc, scale=beta))
	
	
	ax = plt.subplot(111)
	plt.plot([dists,dists],[-1.0,2.0],'k--', lw=2, alpha=1.0, label='best '+title+' detection') 
	# plt.plot([4.115,4.115],[-1.0,2.0],'k--', lw=2, alpha=1.0, label='Leo P detection at 1.74 Mpc')
	# plt.plot([3.897,3.897],[-1.0,2.0],'k-', lw=5, alpha=0.6, label='d=417 kpc')
	# plt.plot([3.317,3.317],[-1.0,2.0],'k-', lw=5, alpha=0.4, label='d=427 kpc')
	plt.ylim(0,1.1)
	plt.xlim(2,12)
	plt.xlabel('$\sigma$ above local mean')
	plt.ylabel('$P(\sigma = X)$')
	plt.legend(loc='best', frameon=True)
	ax.set_aspect(3)
	# plt.show()
	plt.savefig(title+'dist.pdf')
