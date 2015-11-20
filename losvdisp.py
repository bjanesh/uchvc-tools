import scipy.stats as ss
import numpy as np
import matplotlib.pyplot as plt

vgsrall, rgcall = np.loadtxt('vgsrall.dat', usecols=(0,1), unpack=True)
vgsrnosgr, rgcnosgr = np.loadtxt('vgsrnosgr.dat', usecols=(0,1), unpack=True)

vdispall, bedges, binid = ss.binned_statistic(rgcall,vgsrall,statistic=np.std,bins=10,range=[0,60])
vdispnosgr, bedges, binid = ss.binned_statistic(rgcnosgr,vgsrnosgr,statistic=np.std,bins=10,range=[0,60])
bcenters = (bedges[:-1] + bedges[1:]) / 2

plt.scatter(bcenters,vdispall, color='red', edgecolors='none', label='with Sgr')
plt.scatter(bcenters,vdispnosgr, color='blue', edgecolors='none', label='without Sgr')

rgc = np.arange(0,60,0.1)
vlosfit = 111*np.exp(-1*rgc/354)
plt.plot(rgc,vlosfit, linestyle='--', color='black', label='Xue+ 11 fit')
plt.xlim(0,60)
plt.ylim(40,160)
plt.legend()
plt.xlabel(r'$R_{gc}$ (kpc)')
plt.ylabel(r'$\sigma_{los}$ (km s$^{-1}$)')
plt.savefig('test.pdf')