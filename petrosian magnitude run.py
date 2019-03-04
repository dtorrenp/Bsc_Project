# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:39:12 2019

@author: daniel
"""

# -*- coding: utf-8 -*-
#%%
import open_planck as op
import produce_petrosian_gal as pg
import abs_mag as ab

import warnings
import timeit

import matplotlib.pyplot as plt
#%%
start = timeit.default_timer()

#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
#ugriz  filters
filter_one = "petroMag_r"
dot_size = 0.35
bins = 100
magnitude_limit = 20
z_lim = 0.2
number_of_compact_sources = len(t_data_PLANCK)

#%%
if filter_one == "petroMag_u":
    filter_two = "petroMag_g"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "petroMag_g":
    filter_two = "petroMag_r"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "petroMag_r":
    filter_two = "petroMag_i"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "petroMag_i":
    filter_two = "petroMag_z"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
else:
    raise Exception("Invalid Filter, choose one of petroMag_u,petroMag_g,petroMag_r,petroMag_i,petroMag_z ")
#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
#ugriz  filters
with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    SDSS_DATA = pg.loop_through_planck_petro(t_data_PLANCK,filter_one,number_of_compact_sources,magnitude_limit,z_lim)
    length = len(SDSS_DATA[0])
#    
    SDSS_DATA[0].write("Data/sdss_galaxies_data_petrosian_{}___lengthGalaxies_{}".format(filter_one,length), format = 'ascii')
    
#%%

Abs_mag , colour = ab.cal_abs_mag_astro(SDSS_DATA[0],filter_one,filter_two)
    
#%%

plt.figure()
plt.title("Colour Magnitude Petrosian plot")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
plt.scatter(Abs_mag,colour,s = dot_size)
plt.savefig("Plots/Colour_Magnitude_Petrosian_plot_{}_{}___length_{}.png".format(filter_one,filter_two,length))

plt.figure()
plt.title("Colour Magnitude Petrosian Density plot")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
plt.hist2d(Abs_mag,colour, bins = 100)


#%%
end = timeit.default_timer()
print("time taken : {}s".format(round(end-start,2)))