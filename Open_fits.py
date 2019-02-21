# -*- coding: utf-8 -*-
#%%
from astroquery.sdss import SDSS

from astropy.io import fits
from astropy.table import Table, unique
from astropy import coordinates as coords
from astropy.table import vstack
import astropy.units as u

import timeit
import numpy as np
import matplotlib.pyplot as plt

import warnings

#%%
start = timeit.default_timer()

#%%
##Planck data
hdulist_PLANCK = fits.open('/Users/daniel/Documents/Bsc_Project/COM_PCCS-Catalogs_vPR2/COM_PCCS_857_R2.01.fits')
hdu_BinTable_PLANCK = hdulist_PLANCK[1]
t_data_PLANCK = Table(hdu_BinTable_PLANCK.data)
colnames_PLANCK = t_data_PLANCK.colnames

bins = 100

def loop_through_planck(table,f):
    length = 100#len(table)
    
    new_table_galaxy = Table(names=(table.colnames))
    new_table_star = Table(names=(table.colnames))
    new_table_other = Table(names=(table.colnames))
    
    if f == 'r':
        sep_list = []
    mag_list = []
    mag_list_diff_one_two = []
    mag_list_diff_two_three = []
    for i in np.arange(length):
        
        
        righta = table["RA"][i]
        declin = table["DEC"][i]
        
        if np.abs(table["GLAT"][i]) < 20:
            continue
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        

        xid = SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z','class'], photoobj_fields=['ra','dec','u','g','r','i'])
        
        if xid is None:
            continue
        xid = unique(xid,keys=['ra', 'dec'])
        xid.sort(f)
        sub_length = len(xid)        
        
        mag_list = np.concatenate((mag_list,xid[f]))
        
#        if len(xid) == 1:
#            mag_list_diff_one_two.append(xid[f][0])
        if len(xid) == 2:
            mag_list_diff_one_two.append(xid[f][1] - xid[f][0])
        elif len(xid) > 2:
            mag_list_diff_one_two.append(xid[f][1] - xid[f][0])
            mag_list_diff_two_three.append(xid[f][2] - xid[f][1])

        if sub_length == 1:
            number_objects_to_look_at = 1
        elif sub_length == 2:
            number_objects_to_look_at = 1
        else:
            number_objects_to_look_at = 1
        
        for v in np.arange(number_objects_to_look_at):
            max_mag_object = xid[v]
            
            if max_mag_object["class"] == "GALAXY":
                new_table_galaxy = vstack([new_table_galaxy, max_mag_object])
                
                if f == 'r':
                    righta_max = max_mag_object["ra"]
                    declin_max = max_mag_object["dec"]
                    sep = pos.separation(coords.SkyCoord( ra = righta_max*u.degree , dec = declin_max*u.degree))
                    sep_list.append(sep.arcsecond)
                    
            elif max_mag_object["class"] == "STAR":
                new_table_star = vstack([new_table_star, max_mag_object])
            else:
                new_table_other = vstack([new_table_other, max_mag_object])
    
    if f == 'r':
        return new_table_galaxy,new_table_star,new_table_other,sep_list,mag_list,mag_list_diff_one_two,mag_list_diff_two_three
    else:
        return new_table_galaxy,new_table_star,new_table_other
    
#%%

def cal_abs_mag_astro(table,f,f_plus):
    #ugriz  filters
    red_mask = table["Z"] > 0
    red_table = table[red_mask]
    z_val = np.array(red_table["Z"])
    colour = np.asarray(red_table[f_plus]) - np.asarray(red_table[f])
    Abs_mag = np.asarray(red_table[f]) - 5*np.log10(4.28e8*z_val) + 5
    
    # ABSOLUTE MAGNITUDE NOT FEASBLE
    #wathc out for infinity magnitude!!!!!
    return Abs_mag, colour

#%%
#ugriz  filters
with warnings.catch_warnings():
    
    warnings.simplefilter("ignore")
    
#    sdss_galaxies_u, star_u, other_u = loop_through_planck(t_data_PLANCK, "u")
#    sdss_galaxies_u.write("sdss_galaxies_u", format = 'ascii')
#    star_u.write("star_u", format = 'ascii')
#    other_u.write("other_u", format = 'ascii')
#    np.save("separation_u_data",sep_list)
#    Abs_mag_u , colour_u_g = cal_abs_mag_astro(sdss_galaxies_u,"u","g")
#
#    print("1 down")
    
#    sdss_galaxies_g, star_g, other_g = loop_through_planck(t_data_PLANCK, "g")
#    sdss_galaxies_g.write("sdss_galaxies_g", format = 'ascii')
#    star_g.write("star_g", format = 'ascii')
#    other_g.write("other_g", format = 'ascii')
#    Abs_mag_g , colour_g_r = cal_abs_mag_astro(sdss_galaxies_g,"g","r")
#
#    print("2 down")
#    
    sdss_galaxies_r, star_r, other_r, sep_list,mag_list,mag_list_diff_one_two,mag_list_diff_two_three = loop_through_planck(t_data_PLANCK, "r")
#    sdss_galaxies_r.write("sdss_galaxies_r", format = 'ascii')
#    star_r.write("star_r", format = 'ascii')
#    other_r.write("other_r", format = 'ascii')
    Abs_mag_r , colour_r_i = cal_abs_mag_astro(sdss_galaxies_r,"r","i")
#
#    print("3 down")
    
#    sdss_galaxies_i, star_i, other_i = loop_through_planck(t_data_PLANCK, "i")
#    sdss_galaxies_i.write("sdss_galaxies_i", format = 'ascii')
#    star_i.write("star_i", format = 'ascii')
#    other_i.write("other_i", format = 'ascii')
#    Abs_mag_i , colour_i_z = cal_abs_mag_astro(sdss_galaxies_r,"i","z")
#    
#    print("4 down")
    
    #%%

plt.figure(1)
plt.title("Colour Magnitude plot")
plt.xlabel("Absolute magnitude (r)")
plt.ylabel("Colour (u-g)")
plt.scatter(Abs_mag_r,colour_r_i,s = 0.45)
plt.savefig("Colour_Magnitude_plot_r__r_i.png")

#plt.figure(2)
#plt.title("Colour Magnitude plot")
#plt.xlabel("Absolute magnitude (g)")
#plt.ylabel("Colour (g-r)")
#plt.scatter(Abs_mag_g,colour_g_r,s = 0.45)
#plt.savefig("Colour_Magnitude_plot_g__g_r.png")
#    
#plt.figure(3)
#plt.title("Colour Magnitude plot")
#plt.xlabel("Absolute magnitude (r)")
#plt.ylabel("Colour (r-i)")
#plt.scatter(Abs_mag_r,colour_r_i, s = 0.45)

#plt.figure(4)
#plt.title("Colour Magnitude plot")
#plt.xlabel("Absolute magnitude (i)")
#plt.ylabel("Colour (i-z)")
#plt.scatter(Abs_mag_i,colour_i_z, s = 0.45)
#plt.savefig("Colour_Magnitude_plot_i__i_z.png")

#plt.savefig("Colour_Magnitude_plot_r__r_i.png")
#u_data_galaxy = Table.read("sdss_galaxies_u", format='ascii')
#g_data_galaxy = Table.read("sdss_galaxies_g", format='ascii')
#r_data_galaxy = Table.read("sdss_galaxies_r", format='ascii')
#i_data_galaxy = Table.read("sdss_galaxies_i", format='ascii')

#%%

plt.figure(0)
plt.title("magnitude r")
plt.xlabel("magnitude r")
plt.ylabel("Count")
plt.hist(mag_list,bins)
plt.savefig("mag_list")
    
plt.figure(12)
plt.title("mag_list_diff_one_two r")
plt.xlabel("magnitude  diff u")
plt.ylabel("Count")
plt.hist(mag_list_diff_one_two,bins)
plt.savefig("mag_list_diff_one_two")

plt.figure(13)
plt.title("mag_list_diff_two_three r")
plt.xlabel("magnitude  diff r")
plt.ylabel("Count")
plt.hist(mag_list_diff_two_three,bins)
plt.savefig("mag_list_diff_two_three")

plt.figure(5)
plt.title("Galaxy separation using Great Circle - r")
plt.xlabel("Separation (r)")
plt.ylabel("Number of galaxies")
plt.hist(sep_list,bins)
plt.savefig("Galaxy separation using Great Circle - r")

plt.figure(6)
plt.title("Colour Histogram")
plt.xlabel("Colour (r-i)")
plt.ylabel("Count")
plt.hist(colour_r_i,bins)
plt.savefig("colour_r_i")

Colour_mask = (-1.5 < colour_r_i) & (colour_r_i < 0.2)
up_colour = colour_r_i[Colour_mask]
plt.figure()
plt.title("Colour Histogram - reduced")
plt.xlabel("Colour (r-i)")
plt.ylabel("Count")
plt.hist(up_colour,bins)
plt.savefig("reduced colour_r_i")

#plt.figure(7)
#plt.title("Colour Histogram")
#plt.xlabel("Colour (g-r)")
#plt.ylabel("Count")
#plt.hist(colour_g_r)
#
#plt.figure(8)
#plt.title("Colour Histogram")
#plt.xlabel("Colour (r-i)")
#plt.ylabel("Count")
#plt.hist(colour_r_i)

plt.figure(9)
plt.title("Histrogram of absolute magnitudes - r")
plt.xlabel("Absolute magnitude (r)")
plt.ylabel("Number of galaxies")
plt.hist(Abs_mag_r,bins)
plt.savefig("Histrogram of absolute magnitudes - r")

abs_mask = (-18 < Abs_mag_r) & (Abs_mag_r < 10)
up_abs = Abs_mag_r[abs_mask]
plt.figure()
plt.title("Histrogram of absolute magnitudes - r")
plt.xlabel("Absolute magnitude (r)")
plt.ylabel("Number of galaxies")
plt.hist(up_abs,bins)
plt.savefig("Histrogram of absolute magnitudes - r")

#%%
end = timeit.default_timer()
print(end-start)