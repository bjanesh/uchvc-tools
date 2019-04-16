import aplpy
import numpy as np

fig = aplpy.FITSFigure("ch2_fig24a_div.fits")
fig.show_grayscale(vmin=0.803899, vmax=1.13039)
fig.axis_labels.hide()
fig.tick_labels.hide()
fig.ticks.hide()
fig.save("OTA33_med.pdf")
fig.close()

fig = aplpy.FITSFigure("ch2_fig24b.fits")
fig.show_grayscale(vmin=0.803899, vmax=1.13039)
fig.axis_labels.hide()
fig.tick_labels.hide()
fig.ticks.hide()
fig.save("OTA33_med_smooth.pdf")
fig.close()