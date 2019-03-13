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
from astroquery.simbad import Simbad
from astropy.table import Table, unique
from astropy import coordinates as coords
from astropy.table import vstack
import astropy.units as u

import matplotlib.pyplot as plt
#%%
start = timeit.default_timer()

#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
#ugriz  filters
filter_one = "petroMag_u"
filter_two = "petroMag_r"
dot_size = 0.35
bins = 100
magnitude_limit = 20
z_lim = 0.2
number_of_compact_sources = 10#len(t_data_PLANCK)

#Gr (-1.6 < colour) & (colour < 1) & (-18 < Abs_mag) & (Abs_mag < -4)
#Ri (-0.7 < colour) & (colour < 0.5) & (-18.1 < Abs_mag) & (Abs_mag < -4)
#gi (-2 < colour) & (colour < 0.8) & (-17 < Abs_mag) & (Abs_mag < -5)
#ug (-3 < colour) & (colour < 0.5) & (-16 < Abs_mag) & (Abs_mag < -4)
#ur (-4 < colour) & (colour < 1.3) & (-15.8 < Abs_mag) & (Abs_mag < -3)

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

#%%
def loop_through_planck_petro(table,f,length,mag_limit,z_lim):
    
    new_table_galaxy = Table(names=(table.colnames))
    
    for i in np.arange(length):
        
        righta = table["RA"][i]
        declin = table["DEC"][i]
        
        if np.abs(table["GLAT"][i]) < 20:
            continue
        
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        xid = Simbad.query_region(galatic_pos , radius = 2.5*u.arcmin)
        #xid = Simbad.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z','class'], photoobj_fields=['ra','dec','petroMag_u','petroMag_g','petroMag_r','petroMag_i'])
        print(xid.colnames)
        print(xid)
        #SIMBAD IS A COMPLIATION OF A VARIETY OF OTHER CATALOGS NEED TO PICK A SPECIFIC ONE
        
        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        
        mask_redshift = xid['Z'] < z_lim
        xid = xid[mask_redshift]
        
        if xid is None:
            continue
        
        xid.sort(f)
        
        for v in np.arange(len(xid)):
            
            max_mag_object = xid[v]
            max_mag_object = Table(max_mag_object)
            
            if max_mag_object[f] > mag_limit and v > 0:
                break
            
            if max_mag_object["class"] == "GALAXY":
                new_table_galaxy = vstack([new_table_galaxy, max_mag_object])
                    
    return [new_table_galaxy]

#%%
#ugriz  filters
with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    Simbad_DATA = loop_through_planck_petro(t_data_PLANCK,filter_one,number_of_compact_sources,magnitude_limit,z_lim)
    length = len(Simbad_DATA[0])
    
    Simbad_DATA[0].write("Data/Simbad_galaxies_data_petrosian_{}___lengthGalaxies_{}".format(filter_one,length), format = 'ascii')
    
#%%

Abs_mag , colour = ab.cal_abs_mag_astro(Simbad_DATA[0],filter_one,filter_two)

reduc_mask = (-4 < colour) & (colour < 1.3) & (-15.8 < Abs_mag) & (Abs_mag < -3)
reduc_colour = colour[reduc_mask]
reduc_Abs_mag = Abs_mag[reduc_mask]


#CONVERT TO DATAFRAME
df = pd.DataFrame(columns=['colour', 'abs mag'])
df['colour'] = np.array(reduc_colour)
df['abs mag'] = np.array(reduc_Abs_mag)

#%%

plt.figure()
plt.title("Colour Magnitude Petrosian plot")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
plt.scatter(Abs_mag,colour,s = dot_size)
plt.savefig("Plots/Simbad_Colour_Magnitude_Petrosian_plot_{}_{}___length_{}.png".format(filter_one,filter_two,length))

plt.figure()
plt.title("Colour Magnitude Petrosian plot reduced")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
plt.scatter(reduc_Abs_mag,reduc_colour,s = dot_size)
plt.savefig("Plots/Simbad_Colour_Magnitude_Petrosian_plot_reduced_{}_{}___length_{}.png".format(filter_one,filter_two,length))

plt.figure()
plt.title("Colour Magnitude Petrosian Density plot")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
plt.hist2d(Abs_mag,colour, bins = 100)

plt.figure()
plt.title("Colour Magnitude Petrosian Density plot reduced")
plt.xlabel("Absolute magnitude ({})".format(filter_one))
plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
plt.hist2d(reduc_Abs_mag,reduc_colour, bins = bins)

plt.figure()
ax = sns.kdeplot(df['abs mag'], df['colour'], cmap='Reds', shade=False, shade_lowest= False, cbar=True )#.plot_joint(plt.scatter, c='black', s=dot_size, marker='.')
ax.scatter(reduc_Abs_mag,reduc_colour, s = dot_size)

ax2 = sns.jointplot(x='abs mag', y='colour', data=df, kind='kde')
ax2.plot_joint(plt.scatter, c='black', s=dot_size, marker='.')

ax3 = sns.jointplot(x='abs mag', y='colour', data=df,s = dot_size).plot_joint(sns.kdeplot, zorder=1, shade=True, n_levels = 6,colour = 'black')

ax4 = sns.jointplot(x='abs mag', y='colour', data=df,s = dot_size).plot_joint(sns.kdeplot, zorder=1, shade=True, n_levels = 10,colour = 'black')

#%%
end = timeit.default_timer()
print("time taken : {}s".format(round(end-start,2)))