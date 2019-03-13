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
from astropy.io import fits

import matplotlib.pyplot as plt
from scipy.stats import norm
#%%
start = timeit.default_timer()

#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
#ugriz  filters
filter_one = "petroMag_i"
filter_two = "petroMag_r"
dot_size = 0.35
bins = 100
magnitude_limit = 18
mag_lim_sdss = 23
z_lim = 0.2
number_of_compact_sources = len(t_data_PLANCK)
num_req = 4000
width = 10

#Gr (-1.6 < colour) & (colour < 1) & (-18 < Abs_mag) & (Abs_mag < -4)
#Ri (-0.7 < colour) & (colour < 0.5) & (-18.1 < Abs_mag) & (Abs_mag < -4)
#gi (-2 < colour) & (colour < 0.8) & (-17 < Abs_mag) & (Abs_mag < -5)
#ug (-3 < colour) & (colour < 0.5) & (-16 < Abs_mag) & (Abs_mag < -4)
#ur (-4 < colour) & (colour < 1.3) & (-15.8 < Abs_mag) & (Abs_mag < -3)
#(-3.6 < colour_sdss) & (colour_sdss < -1) & (-15.8 < Abs_mag_sdss) & (Abs_mag_sdss < -11)

#%%
#if filter_one == "petroMag_u":
#    filter_two = "petroMag_g"
#elif filter_one == "petroMag_g":
#    filter_two = "petroMag_r"
#elif filter_one == "petroMag_r":
#    filter_two = "petroMag_i"
#elif filter_one == "petroMag_i":
#    filter_two = "petroMag_z"
#else:
#    raise Exception("Invalid Filter, choose one of petroMag_u,petroMag_g,petroMag_r,petroMag_i,petroMag_z ")
#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()
print("PLANCK OPEN")

#%%
def loop_through_planck_petro(table,f,length,mag_limit,z_lim):
    
    new_table_galaxy = Table()
    
    for i in np.arange(length):
        
        righta = table["RA"][i]
        declin = table["DEC"][i]
        
        if np.abs(table["GLAT"][i]) < 20:
            continue
        
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic

        xid = SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z','class', 'BestObjID'], photoobj_fields=['ra','dec','petroMag_u','petroMag_g','petroMag_r','petroMag_i'])

        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        
        mask_redshift_gal = (xid['Z'] < z_lim) & (xid["class"] == "GALAXY")
        xid = xid[mask_redshift_gal]
        
        if len(xid) == 0:
            continue
        
        
        data_mask = xid[f] > mag_limit
        data = xid[data_mask]
        
        if len(data) == 0:
            xid.sort(f)
            data = Table(xid[0])
            
        new_table_galaxy = vstack([new_table_galaxy, data])
        print(len(new_table_galaxy))
                    
    return [new_table_galaxy]

#%%
    
#def loop_through_sdss_gal(f, num_req, mag_limit, width):
#
#    sdss_gal = Table()
#    while len(sdss_gal) < num_req:
#        
#        righta = np.random.uniform(-180,180)
#        declin = np.random.uniform(-90,90)
#        
#        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
#        
#        galatic_pos = pos.galactic
#        
#        xid = SDSS.query_region(galatic_pos, spectro = True , radius = width*u.arcmin , specobj_fields=['Z','class', 'BestObjID'], photoobj_fields=['ra','dec','petroMag_u','petroMag_g','petroMag_r','petroMag_i'])
#        
#        if xid is None:
#            continue
#        
#        xid = unique(xid,keys=['ra', 'dec'])
#        xid.sort(f)
#        
#        xid_gal_mask = (xid["class"] == "GALAXY") & (xid[f] > mag_limit)
#        xid_reduc = xid[xid_gal_mask]
#        sdss_gal = vstack([sdss_gal, xid_reduc])
#        print(len(sdss_gal))
#        
#    sdss_gal = unique(sdss_gal,keys=['ra', 'dec'])
#    return[sdss_gal]
#    
#%%
    

#hdulist_PLANCK = fits.open('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/COM_PCCS_857_R2.01.fits')
#hdulist_SDSS = fits.open('/Users/daniel/Documents/Bsc_Project/random_galaxies_model_mag.fits')
hdulist_SDSS = fits.open('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/local code/petromag_galaxies.fits')
#hdulist_PLANCK = fits.open('C:/Users/DELL/Desktop/COM_PCCS_857_R2.01.fits')
hdu_BinTable_SDSS = hdulist_SDSS[1]
t_data_SDSS = Table(hdu_BinTable_SDSS.data)
print("SDSS OPEN")

#%%
#ugriz  filters
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    
    #SDSS_DATA = loop_through_planck_petro(t_data_PLANCK,filter_one,number_of_compact_sources,magnitude_limit,z_lim)
    #length = len(SDSS_DATA[0])
    #SDSS_DATA[0].write("Data/sdss_galaxies_data_petrosian_{}___lengthGalaxies_{}".format(filter_one,length), format = 'ascii')
    
    SDSS_DATA = [Table.read('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/local code/Data/sdss_galaxies_data_petrosian_petroMag_i___lengthGalaxies_1914', format = 'ascii')]
    length = len(SDSS_DATA[0])
    print("CROSS MATCH OPEN")
    #sdss_galaxies = loop_through_sdss_gal(filter_one, num_req,mag_lim_sdss, width)
#    print("HIEEEE")
    
    
    
    #sdss_galaxies[0].write("Data/sdss_galaxies_data_petrosian_sdssgal_{}___lengthGalaxies_{}".format(filter_one,len(sdss_galaxies)), format = 'ascii')
    
#%%
    
Planck_table = ab.cal_abs_mag_astro(SDSS_DATA[0],filter_one,filter_two)
Abs_mag = Planck_table['Absolute Magnitude {}'.format(filter_one)]
colour = Planck_table['Colour ({} - {})'.format(filter_two,filter_one)]

SDSS_table = ab.cal_abs_mag_astro_sdss(t_data_SDSS,filter_one,filter_two)
Abs_mag_sdss = SDSS_table['Absolute Magnitude {}'.format(filter_one)]
colour_sdss = SDSS_table['Colour ({} - {})'.format(filter_two,filter_one)]

print("ABSOLUTE MAGNITUDE CALCULATED")



#%%

reduc_mask = (0.15 < colour) & (colour < 0.6) & (-18 < Abs_mag) & (Abs_mag < -12)

reduc_mask_sdss = (0.15 < colour_sdss) & (colour_sdss < 0.6) & (-18 < Abs_mag_sdss) & (Abs_mag_sdss < -12)

reduc_colour = colour[reduc_mask]
reduc_Abs_mag = Abs_mag[reduc_mask]

reduc_colour_sdss = colour_sdss[reduc_mask_sdss]
reduc_Abs_mag_sdss = Abs_mag_sdss[reduc_mask_sdss]

green_valley_galaxies_mask = (0.15 < colour) & (colour < 0.6) & (-18 < Abs_mag) & (Abs_mag < -12)
green_valley_galaxies = Planck_table[green_valley_galaxies_mask]
green_valley_galaxies.write("green valley candidates_{}_filter ".format(filter_one), format = 'ascii')
#print(green_valley_galaxies["BestObjID"])

print("DATA SET REDUCED")

#%%
#CONVERT TO DATAFRAME
df_pre = pd.DataFrame(columns=['colour', 'abs mag'])
df_pre['colour'] = np.array(colour)
df_pre['abs mag'] = np.array(Abs_mag)

df_sdss = pd.DataFrame(columns=['colour', 'abs mag'])
df_sdss['colour'] = np.array(reduc_colour_sdss)
df_sdss['abs mag'] = np.array(reduc_Abs_mag_sdss)


df = pd.DataFrame(columns=['colour', 'abs mag'])
df['colour'] = np.array(reduc_colour)
df['abs mag'] = np.array(reduc_Abs_mag)

df_sdss = pd.DataFrame(columns=['colour', 'abs mag'])
df_sdss['colour'] = np.array(reduc_colour_sdss)
df_sdss['abs mag'] = np.array(reduc_Abs_mag_sdss)

print("DATA FRAME COMPLETE")


#%%
plt.figure()
plt.title("SDSS galaxy contours with Planck data")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
plt.scatter(Abs_mag,colour, s = dot_size)
plt.savefig("Plots/Colour_Magnitude_Petrosian_plot_UNEDITED_SDSS_GALAXIES_{}_{}___length_{}.png".format(filter_one,filter_two,length))

#%%
plt.figure()
plt.title("SDSS galaxy contours with Planck data")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
ax = sns.kdeplot(df_sdss['abs mag'], df_sdss['colour'], cmap='Reds', shade=False, shade_lowest= False, cbar=True, n_levels = 15 )#.plot_joint(plt.scatter, c='black', s=dot_size, marker='.')
ax.scatter(reduc_Abs_mag,reduc_colour, s = dot_size)
plt.savefig("Plots/Colour_Magnitude_Petrosian_plot_contour_SDSS_GALAXIES_{}_{}___length_{}.png".format(filter_one,filter_two,length))

plt.figure()
plt.title("Colour Histogram SDSS vs Planck- ({} - {})".format(filter_one,filter_two))
plt.xlabel("Colour ({} - {})".format(filter_one,filter_two))
plt.ylabel("Count")
plt.hist(reduc_colour_sdss, bins = bins, alpha = 0.5)
plt.hist(reduc_colour,bins = bins, alpha = 0.5)
plt.savefig("Plots/Colour Histogram_SDSS_vs_Planck_{}_{}___length_{}_bins_{}.png".format(filter_one,filter_two,length,bins))

#%%
#def gauss(x):
#    return np.exp( -0.5 * x**2 / np.sqrt( 2 * np.pi ) )
#
#popt, pcov = curve_fit(gauss, xdata, ydata)
#plt.plot(xdata, func(xdata, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

plt.figure()
plt.title("Colour Histogram Planck- ({} - {})".format(filter_one,filter_two))
plt.xlabel("Colour ({} - {})".format(filter_one,filter_two))
plt.ylabel("Count")
ax = sns.distplot(reduc_colour,fit=norm, kde=False, bins = bins)
#plt.hist(reduc_colour,bins = bins)
plt.savefig("Plots/Colour Histogram_SDSS_vs_Planck__PLANCK_{}_{}___length_{}_bins_{}.png".format(filter_one,filter_two,length,bins))

plt.figure()
plt.title("Colour Histogram SDSS ({} - {})".format(filter_one,filter_two))
plt.xlabel("Colour ({} - {})".format(filter_one,filter_two))
plt.ylabel("Count")
#ax = sns.distplot(reduc_colour_sdss,fit=norm, kde=False, bins = bins)
plt.hist(reduc_colour_sdss, bins = bins)
plt.savefig("Plots/Colour Histogram_SDSS_vs_Planck_SDSS_{}_{}___length_{}_bins_{}.png".format(filter_one,filter_two,length,bins))

print("PLOTS DONE")

#%%
end = timeit.default_timer()
print("time taken : {}s".format(round(end-start,2)))