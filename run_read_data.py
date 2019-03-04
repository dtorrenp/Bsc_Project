# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 09:45:37 2019

@author: daniel
"""
# -*- coding: utf-8 -*-
#%%

import abs_mag as ab
import produce_plots as pp
import top_ten_plate_pics as tt
import timeit

from astropy.table import Table

#%%
start = timeit.default_timer()

#%%
#ugriz  filters
filter_one = "u"
dot_size = 0.35
bins = 100
magnitude_limit = 20
length = 2000

#%%

COLOUR_MAGNITUDE_PLOT = True
RED_SHIFT_PLOT = True
GALAXY_SEPARATION_PLOT = True
MAGNITUDE_PLOT = True
MAGNITUDE_DIFFERENCE_ONE_TWO_PLOT = True
MAGNITUDE_DIFFERENCE_TWO_THREE_PLOT = True
COLOUR_PLOT = True
COLOUR_ABS_MAGNITUDE_PLOT = True
GALAXY_MAGNITUDE_PLOT = True
GALAXY_PLATE_PICS  = False

#COLOUR_MAGNITUDE_PLOT = False
#RED_SHIFT_PLOT = False
#GALAXY_SEPARATION_PLOT = False
#MAGNITUDE_PLOT = False
#MAGNITUDE_DIFFERENCE_ONE_TWO_PLOT = False
#MAGNITUDE_DIFFERENCE_TWO_THREE_PLOT = False
#COLOUR_PLOT = False
#COLOUR_ABS_MAGNITUDE_PLOT = False
#GALAXY_MAGNITUDE_PLOT = False
#GALAXY_PLATE_PICS  = True

#COLOUR_MAGNITUDE_PLOT = False
#RED_SHIFT_PLOT = False
#GALAXY_SEPARATION_PLOT = False
#MAGNITUDE_PLOT = False
#MAGNITUDE_DIFFERENCE_ONE_TWO_PLOT = False
#MAGNITUDE_DIFFERENCE_TWO_THREE_PLOT = False
#COLOUR_PLOT = False
#COLOUR_ABS_MAGNITUDE_PLOT = False
#GALAXY_MAGNITUDE_PLOT = False
#GALAXY_PLATE_PICS  = False

#%%
if filter_one == "u":
    filter_two = "g"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "g":
    filter_two = "r"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "r":
    filter_two = "i"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "i":
    filter_two = "z"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
else:
    raise Exception("Invalid Filter, choose one of u,g,r,i,z ")
#%%
    
SDSS_DATA = Table.read("Data/sdss_data_{}__length_{}".format(filter_one,length), format = 'ascii')

#%%
Abs_mag , colour = ab.cal_abs_mag_astro(SDSS_DATA[0],filter_one,filter_two)

#def reduced_data_set(DATA):
#    galaxies_copy = DATA[0]
#    galaxy_mask = Abs_mag < -10
#    reduced_data_set = galaxies_copy[galaxy_mask]
#    x,y = ab.cal_abs_mag_astro(reduced_data_set,filter_one,filter_two)
#    pp.produce_colour_magnitude_reduc(x,y,filter_one,filter_two,dot_size)

#%%
if COLOUR_MAGNITUDE_PLOT == True:
    pp.produce_colour_magnitude(Abs_mag,colour,filter_one,filter_two,dot_size)
if RED_SHIFT_PLOT == True:
    pp.produce_red_shift(SDSS_DATA[0],filter_one,bins)
if GALAXY_SEPARATION_PLOT == True:
    pp.produce_galaxy_separation(SDSS_DATA[3],filter_one,bins)
if MAGNITUDE_PLOT == True:
    pp.produce_magnitude(SDSS_DATA[4],filter_one,bins)
if MAGNITUDE_DIFFERENCE_ONE_TWO_PLOT == True:
    pp.produce_mag_difference_one_two(SDSS_DATA[5],filter_one,bins)
if MAGNITUDE_DIFFERENCE_TWO_THREE_PLOT == True:
    pp.produce_mag_difference_two_three(SDSS_DATA[6],filter_one,bins)
if COLOUR_PLOT == True:
    pp.produce_colour(colour,filter_one,filter_two,bins,lim_low_mag,lim_up_mag)
if COLOUR_ABS_MAGNITUDE_PLOT == True:
    pp.produce_abs_mag(Abs_mag,filter_one,bins,lim_low_abs_mag,lim_up_abs_mag)
if GALAXY_MAGNITUDE_PLOT == True:
    pp.produce_galaxy_magnitude(SDSS_DATA[0][filter_one],filter_one,bins)

#%%
if GALAXY_PLATE_PICS == True:
    tt.top_ten_pics(SDSS_DATA[0],filter_one)
    
#%%
end = timeit.default_timer()
print("time taken : {}s".format(round(end-start,2)))

