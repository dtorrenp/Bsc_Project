# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 11:41:37 2019

@author: daniel
"""
#%%

from astroquery.sdss import SDSS
from astropy.table import Table, unique
from astropy import coordinates as coords
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

import open_planck as op
import timeit
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

def loop_through__sep(table,f,length,mag_limit,z_lim):
    sep_list = []
    
    for i in np.arange(length):
        
        righta = table["RA"][i]
        declin = table["DEC"][i]
        
        if np.abs(table["GLAT"][i]) < 20:
            continue
        
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic

        xid = SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z','class'], photoobj_fields=['ra','dec','modelMag_u','modelMag_g','modelMag_r','modelMag_i'])

        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        
        xid.sort(f)
            
        mask_redshift_gal = (xid['Z'] < z_lim) & (xid["class"] == "GALAXY")
        xid = xid[mask_redshift_gal]
        
        if len(xid) == 0:
            continue
        
        xid.sort(f)
        data_mask = xid[f] > mag_limit
        data = xid[data_mask]
        
        if len(data) == 0:
            data = Table(xid[0])
        
        righta_max = data["ra"]
        declin_max = data["dec"]
        sep = pos.separation(coords.SkyCoord( ra = righta_max*u.degree , dec = declin_max*u.degree))
        sep_list = np.concatenate((sep_list,sep.arcsecond))

        print(len(sep_list))
                
    return sep_list

#%%
    
def rand_loop_through__sep(f,length):
    sep_list = []
    length_max = length 
    
    while len(sep_list) < length_max:
        
        righta = np.random.uniform(-180,180)
        declin = np.random.uniform(-90,90)
        
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        
        xid = SDSS.query_region(galatic_pos, spectro = False , radius = 2.5*u.arcmin , photoobj_fields=['ra','dec'])
        
        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        
        selection = np.unique(np.random.randint(0,len(xid),size=10))
        
        max_mag_object = xid[selection]
        
        righta_max = max_mag_object["ra"]
        declin_max = max_mag_object["dec"]
        sep = pos.separation(coords.SkyCoord( ra = righta_max*u.degree , dec = declin_max*u.degree))
        sep_list = np.concatenate((sep_list,sep.arcsecond))
        print(len(sep_list))
        
    return sep_list

#%%
#a = loop_through__sep(t_data_PLANCK,filter_one,number_of_compact_sources,magnitude_limit,z_lim)
#np.save("Data/separation_test_{}___lengthGalaxies_{}".format(filter_one,len(a)),a)
a = np.load("Data/separation_test_modelMag_r___lengthGalaxies_2460.npy")
print("CROSS MATCH OPEN")
#b = rand_loop_through__sep(filter_one, len(a))
#np.save("Data/separation_test_rand_{}___lengthGalaxies_{}".format(filter_one,len(b)),b)
b = np.load("Data/separation_test_rand_modelMag_r___lengthGalaxies_2462.npy")
print("RAND OPEN")

#%%
length = len(a)
plt.figure()
plt.title("separation of galaxies")
plt.xlabel("Separation")
plt.ylabel("Count")
plt.hist(a,bins, label = "Planck", alpha = 0.6)
plt.hist(np.asarray(b),bins,label = "SDSS random", alpha = 0.5)
plt.legend()
plt.savefig("Plots/separation random together Historgram_{}___length_{}".format(filter_one,length))

plt.figure()
plt.title("Separation of galaxies - Planck")
plt.xlabel("Separation")
plt.ylabel("Count")
plt.hist(a,bins)
plt.savefig("Plots/separation Planck random Historgram_{}___length_{}".format(filter_one,length))

plt.figure()
plt.title("Separation of galaxies - Random")
plt.xlabel("Separation")
plt.ylabel("Count")
plt.hist(np.asarray(b),bins)
plt.savefig("Plots/separation random Historgram_{}___length_{}".format(filter_one,length))

print("PLOTS DONE")

#%%
end = timeit.default_timer()
print("time taken : {}s".format(round(end-start,2)))
