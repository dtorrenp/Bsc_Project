# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 11:47:33 2019

@author: daniel
"""
#%%
import matplotlib.pyplot as plt

#%%#TURNS OFF PICTURES PRODUCED IN CONSOLE, LOOK IN YOUR FOLDER

#%%
def produce_colour_magnitude(Abs_mag,colour,filter_one,filter_two,dot_size):
    length = len(Abs_mag)
    plt.figure()
    plt.title("Colour Magnitude plot")
    plt.xlabel("Absolute magnitude ({})".format(filter_one))
    plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
    plt.scatter(Abs_mag,colour,s = dot_size)
    plt.savefig("Plots/Colour_Magnitude_plot_{}_{}___length_{}.png".format(filter_one,filter_two,length))

#%%
def produce_magnitude(mag_list,filter_one,bins):
    length = len(mag_list)
    plt.figure()
    plt.title("Magnitude ({})".format(filter_one))
    plt.xlabel("Magnitude ({})".format(filter_one))
    plt.ylabel("Count")
    plt.hist(mag_list,bins)
    plt.savefig("Plots/Magnitude Historgram_{}___length_{}_bins_{}".format(filter_one,length,bins))
    
#%%
def produce_red_shift(galaxies,filter_one,bins):
    length = len(galaxies)
    plt.figure()
    plt.title("Red shift of galaxies")
    plt.xlabel("Red_shift")
    plt.ylabel("Count")
    plt.hist(galaxies["Z"],bins)
    plt.savefig("Plots/Redshift Historgram_{}___length_{}_bins_{}".format(filter_one,length,bins))
    
#%%
def produce_mag_difference_one_two(mag_list_diff_one_two,filter_one,bins):
    length = len(mag_list_diff_one_two)
    plt.figure()
    plt.title("Difference in Magnitude between first and second object {}".format(filter_one))
    plt.xlabel("Magnitude  Difference ({})".format(filter_one))
    plt.ylabel("Count")
    plt.hist(mag_list_diff_one_two,bins)
    plt.savefig("Plots/Difference in Magnitude between first and second object Historgram_{}___length_{}_bins_{}".format(filter_one,length,bins))
    
#%%
def produce_mag_difference_two_three(mag_list_diff_two_three,filter_one,bins):
    length = len(mag_list_diff_two_three)
    plt.figure()
    plt.title("Difference in Magnitude between second and third object {}".format(filter_one))
    plt.xlabel("Magnitude  Difference ({})".format(filter_one))
    plt.ylabel("Count")
    plt.hist(mag_list_diff_two_three,bins)
    plt.savefig("Plots/Difference in Magnitude between second and third object Historgram_{}___length_{}_bins_{}".format(filter_one,length,bins))
    
#%%
def produce_galaxy_separation(sep_list,filter_one,bins):
    length = len(sep_list)
    plt.figure()
    plt.title("Galaxy separation using Great Circle - {}".format(filter_one))
    plt.xlabel("Separation ({})".format(filter_one))
    plt.ylabel("Count")
    plt.hist(sep_list,bins)
    plt.savefig("Plots/Galaxy separation using Great Circle Historgram_{}___length_{}_bins_{}".format(filter_one,length,bins))
    
#%%
def produce_colour(colour,filter_one,filter_two,bins,lim_low,lim_up):
    length = len(colour)
    plt.figure()
    plt.title("Colour Histogram - ({} - {})".format(filter_one,filter_two))
    plt.xlabel("Colour ({} - {})".format(filter_one,filter_two))
    plt.ylabel("Count")
    plt.hist(colour,bins)
    plt.savefig("Plots/Colour Histogram_{}_{}___length_{}_bins_{}.png".format(filter_one,filter_two,length,bins))
    
    Colour_mask = (lim_low < colour) & (colour < lim_up)
    reduc_colour = colour[Colour_mask]
    reduc_length = len(reduc_colour)
    plt.figure()
    plt.title("Reduced Colour Histogram - ({} - {})".format(filter_one,filter_two))
    plt.xlabel("Colour ({} - {})".format(filter_one,filter_two))
    plt.ylabel("Count")
    plt.hist(reduc_colour,bins)
    plt.savefig("Plots/Reduced Colour Histogram_{}_{}___length_{}_bins_{}_limits_{}_{}.png".format(filter_one,filter_two,reduc_length,bins,lim_low,lim_up))
    
#%%
def produce_abs_mag(abs_mag,filter_one,bins,lim_low,lim_up):
    length = len(abs_mag)
    plt.figure()
    plt.title("Absolute Magnitude Histogram - ({})".format(filter_one))
    plt.xlabel("Absolute Magnitude ({})".format(filter_one))
    plt.ylabel("Count")
    plt.hist(abs_mag,bins)
    plt.savefig("Plots/Absolute Magnitude Histogram_{}___length_{}_bins_{}.png".format(filter_one,length,bins))
    
    abs_mag_mask = (lim_low < abs_mag) & (abs_mag < lim_up)
    reduc_abs_mag = abs_mag[abs_mag_mask]
    reduc_length = len(reduc_abs_mag)
    plt.figure()
    plt.title("Reduced Absolute Magnitude Histogram - ({})".format(filter_one))
    plt.xlabel("Absolute Magnitude ({})".format(filter_one))
    plt.ylabel("Count")
    plt.hist(reduc_abs_mag,bins)
    plt.savefig("Plots/Reduced Absolute Magnitude Histogram_{}___length_{}_bins_{}_limits_{}_{}.png".format(filter_one,reduc_length,bins,lim_low,lim_up))
    
#%%
def produce_galaxy_magnitude(mag_list,filter_one,bins):
    length = len(mag_list)
    plt.figure()
    plt.title("Magnitude Galaxies ({})".format(filter_one))
    plt.xlabel("Magnitude ({})".format(filter_one))
    plt.ylabel("Count")
    plt.hist(mag_list,bins)
    plt.savefig("Plots/Magnitude Galaxies Historgram_{}___length_{}_bins_{}".format(filter_one,length,bins))