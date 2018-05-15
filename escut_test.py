import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from photutils.datasets import make_random_gaussians_table, make_gaussian_sources_image


def genStars(fwhm):
    stdev = fwhm/2.355
    n_sources = 100
    param_ranges = [('amplitude', [1000,35000]),
                    ('x_mean', [0, 500]),
                    ('y_mean', [0, 300]),
                    ('x_stddev', [stdev-0.1,stdev+0.1]),
                    ('y_stddev', [stdev-0.1,stdev+0.1]),
                    ('theta', [0, np.pi])]
    param_ranges = OrderedDict(param_ranges)
    sources = make_random_gaussians_table(n_sources, param_ranges,
                                          random_state=12345)
    return sources
    
def genSources():
    n_sources = 100
    param_ranges = [('amplitude', [1000,10000]),
                    ('x_mean', [0, 500]),
                    ('y_mean', [0, 300]),
                    ('x_stddev', [3,12]),
                    ('y_stddev', [3,12]),
                    ('theta', [0, np.pi])]
    param_ranges = OrderedDict(param_ranges)
    sources = make_random_gaussians_table(n_sources, param_ranges,
                                          random_state=12345)
    return sources

def main():
    stars = genStars(6.0)
    sources = genSources()
    shape = (300, 500)
    image1 = make_gaussian_sources_image(shape, stars)
    image2 = make_gaussian_sources_image(shape, sources)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    ax1.imshow(image1, origin='lower', interpolation='nearest')
    ax2.imshow(image2, origin='lower', interpolation='nearest')
    plt.show()

if __name__ == '__main__':
    main()

