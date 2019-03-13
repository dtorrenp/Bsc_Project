# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:39:12 2019

@author: daniel
"""

# -*- coding: utf-8 -*-
#%%
import open_planck as op
import abs_mag as ab
import pandas as pd
import warnings
import numpy as np
import timeit
import seaborn as sns
from astroquery.sdss import SDSS
from astropy.table import Table, unique
from astropy import coordinates as coords
from astropy.table import vstack
import astropy.units as u

import matplotlib.pyplot as plt
#%%
start = timeit.default_timer()

#%%
#ugriz  filters
filter_one = "petroMag_u"
filter_two = "petroMag_r"
dot_size = 0.35
bins = 100
magnitude_limit = 18
mag_lim_sdss = 23
z_lim = 0.2
number_of_compact_sources = len(t_data_PLANCK)
num_req = 4000
width = 10



t = Table.read('photometry.dat', format='ascii.daophot')
#%%
Abs_mag , colour = ab.cal_abs_mag_astro(SDSS_DATA[0],filter_one,filter_two)
Abs_mag_sdss , colour_sdss = ab.cal_abs_mag_astro(sdss_galaxies[0],filter_one,filter_two)

reduc_mask = (-4 < colour) & (colour < 1.3) & (-15.8 < Abs_mag) & (Abs_mag < -3)

reduc_mask_sdss = (-4 < colour_sdss) & (colour_sdss < 1.3) & (-15.8 < Abs_mag_sdss) & (Abs_mag_sdss < -3)

reduc_colour = colour[reduc_mask]
reduc_Abs_mag = Abs_mag[reduc_mask]

reduc_colour_sdss = colour_sdss[reduc_mask_sdss]
reduc_Abs_mag_sdss = Abs_mag_sdss[reduc_mask_sdss]

#CONVERT TO DATAFRAME
df = pd.DataFrame(columns=['colour', 'abs mag'])
df['colour'] = np.array(reduc_colour)
df['abs mag'] = np.array(reduc_Abs_mag)

df_sdss = pd.DataFrame(columns=['colour', 'abs mag'])
df_sdss['colour'] = np.array(reduc_colour_sdss)
df_sdss['abs mag'] = np.array(reduc_Abs_mag_sdss)

#%%

plt.figure()
plt.title("SDSS galaxy contours with Planck data")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
ax = sns.kdeplot(df_sdss['abs mag'], df_sdss['colour'], cmap='Reds', shade=False, shade_lowest= False, cbar=True )#.plot_joint(plt.scatter, c='black', s=dot_size, marker='.')
ax.scatter(reduc_Abs_mag,reduc_colour, s = dot_size)
plt.savefig("Plots/Colour_Magnitude_Petrosian_plot_contour_SDSS_GALAXIES_{}_{}___length_{}.png".format(filter_one,filter_two,length))

plt.figure()
plt.title("Colour Histogram SDSS vs Planck- ({} - {})".format(filter_one,filter_two))
plt.xlabel("Colour ({} - {})".format(filter_one,filter_two))
plt.ylabel("Count")
plt.hist(reduc_colour,bins = bins)
plt.hist(reduc_colour_sdss, bins = bins)
plt.savefig("Plots/Colour Histogram_SDSS_vs_Planck_{}_{}___length_{}_bins_{}.png".format(filter_one,filter_two,length,bins))