import aplpy
import numpy as np

fig = aplpy.FITSFigure("AGC198606_i.fits")
fig.show_grayscale()
fig.save("plain.pdf")
fig.close()

fig = aplpy.FITSFigure("AGC198606_i.fits")
fig.show_grayscale()

x1_bg, y1_bg, x2_bg, y2_bg = np.loadtxt("bgvals_AGC198606_i.fits.txt", usecols=(0,1,2,3), unpack=True)
xc_bg = x1_bg+(x2_bg-x1_bg)/2.0
yc_bg = y1_bg+(y2_bg-y1_bg)/2.0

rac_bg, decc_bg = fig.pixel2world(xc_bg,yc_bg)

fig.show_rectangles(rac_bg, decc_bg, 0.0003, 0.0003, color='red', lw=0.25)

fig.save("bgboxes.pdf")
fig.close()

fig = aplpy.FITSFigure("AGC198606_i.fits")
fig.show_grayscale()

x_dao, y_dao = np.loadtxt("AGC198606_i.fits.coo.1", usecols=(0,1), unpack=True)
ra_dao, dec_dao = fig.pixel2world(x_dao,y_dao)

fig.show_markers(ra_dao, dec_dao, marker='x', facecolor='red', lw=0.25, s=10)

fig.save("daofind.pdf")
fig.close()

fig = aplpy.FITSFigure("AGC198606_i.fits")
fig.show_grayscale()

x_mat, y_mat = np.loadtxt("tol7_i.pos", usecols=(0,1), unpack=True)
ra_mat, dec_mat = fig.pixel2world(x_mat,y_mat)

fig.show_markers(ra_mat, dec_mat, marker='x', facecolor='red', lw=0.25, s=10)

fig.save("match.pdf")
fig.close()

fig = aplpy.FITSFigure("AGC198606_i.fits")
fig.show_grayscale()

x_apcor, y_apcor = np.loadtxt("apcor_stars_i.txt", usecols=(0,1), unpack=True)
ra_apcor, dec_apcor = fig.pixel2world(x_apcor,y_apcor)

fig.show_markers(ra_apcor, dec_apcor, marker='o', edgecolor='red', facecolor='none', lw=1, s=30)

fig.save("apcor.pdf")
fig.close()

fig = aplpy.FITSFigure("AGC198606_i.fits")
fig.show_grayscale()

x_escut, y_escut = np.loadtxt("escut_i.pos", usecols=(0,1), unpack=True)
ra_escut, dec_escut = fig.pixel2world(x_escut,y_escut)

fig.show_markers(ra_escut, dec_escut, marker='x', facecolor='red', lw=0.25, s=10)

fig.save("escut.pdf")
fig.close()