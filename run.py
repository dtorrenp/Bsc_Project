# -*- coding: utf-8 -*-
#%%
import open_planck as op
import produce_galaxies_stars_other as pg
import abs_mag as ab
import produce_plots as pp
import top_ten_plate_pics as tt
import warnings
import timeit
import numpy as np
#%%
start = timeit.default_timer()

#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
#ugriz  filters
filter_one = "modelMag_r"
dot_size = 0.35
bins = 100
magnitude_limit = 20
z_lim = 0.6
number_of_compact_sources = len(t_data_PLANCK)

#%%

COLOUR_MAGNITUDE_PLOT = True
COLOUR_MAGNITUDE_AP_TYPE_PLOT = True
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
#COLOUR_MAGNITUDE_AP_TYPE_PLOT = False
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
#COLOUR_MAGNITUDE_AP_TYPE_PLOT = False
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
if filter_one == "modelMag_u":
    filter_two = "modelMag_g"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "modelMag_g":
    filter_two = "modelMag_r"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "modelMag_r":
    filter_two = "modelMag_i"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
elif filter_one == "modelMag_i":
    filter_two = "modelMag_z"
    lim_low_mag = -1.5
    lim_up_mag = 0.2
    lim_low_abs_mag = -18
    lim_up_abs_mag = 10
else:
    raise Exception("Invalid Filter, choose one of modelMag_u,modelMag_g,modelMag_r,modelMag_i,modelMag_z ")
#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
#ugriz  filters
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    #loop_through_planck(table,f,length)
    SDSS_DATA = pg.loop_through_planck(t_data_PLANCK,filter_one,number_of_compact_sources,magnitude_limit,z_lim)
    #[new_table_galaxy,new_table_star,new_table_other,sep_list,mag_list,mag_list_diff_one_two,mag_list_diff_two_three]
    length = len(SDSS_DATA[0])
    
    #SDSS_DATA.write("Data/sdss_data_{}__length_{}".format(filter_one,length), format = 'ascii')
#    
    SDSS_DATA[0].write("Data/sdss_galaxies_data_{}___lengthGalaxies_{}".format(filter_one,length), format = 'ascii')
    SDSS_DATA[1].write("Data/star_{}___lengthGalaxies_{}".format(filter_one,length), format = 'ascii')
    SDSS_DATA[2].write("Data/other_{}___lengthGalaxies_{}".format(filter_one,length), format = 'ascii')
    np.save("Data/separation_{}___lengthGalaxies_{}".format(filter_one,length),SDSS_DATA[3])
    
#%%
galaxies = SDSS_DATA[0]
exp_mask = galaxies["Aperture Type"] == 1
de_v_mask = np.invert(exp_mask)

exp_galaxies = galaxies[exp_mask]
de_v_galaxies = galaxies[de_v_mask]

Abs_mag_exp , colour_exp = ab.cal_abs_mag_astro(exp_galaxies,filter_one,filter_two)
Abs_mag_de_v , colour_de_v = ab.cal_abs_mag_astro(de_v_galaxies,filter_one,filter_two)

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
if COLOUR_MAGNITUDE_AP_TYPE_PLOT == True:
    pp.produce_colour_magnitude_ap_type(Abs_mag_exp,colour_exp,Abs_mag_de_v,colour_de_v,filter_one,filter_two,dot_size)
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
