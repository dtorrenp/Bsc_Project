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
#import produce_galaxies_stars_other as pg
#import abs_mag as ab
#import produce_plots as pp
#import top_ten_plate_pics as tt
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
            
        mask_redshift = xid['Z'] < z_lim
        xid = xid[mask_redshift]
        
        if xid is None:
            continue
        
        for v in np.arange(len(xid)):
            
            max_mag_object = xid[v]
            max_mag_object = Table(max_mag_object)
            
            if max_mag_object[f] > mag_limit and v > 0:
                break
            
            if max_mag_object["class"] == "GALAXY":
                righta_max = max_mag_object["ra"]
                declin_max = max_mag_object["dec"]
                sep = pos.separation(coords.SkyCoord( ra = righta_max*u.degree , dec = declin_max*u.degree))
                sep_list.append(sep.arcsecond[0])
                
    return sep_list

#%%
    
def rand_loop_through__sep(f,length):
    sep_list = []
    length_max = length 
    print(length_max)
    
    while len(sep_list) < length_max:
        
        righta = np.random.uniform(-180,180)
        declin = np.random.uniform(-90,90)
        
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        
        xid = SDSS.query_region(galatic_pos, spectro = False , radius = 2.5*u.arcmin , photoobj_fields=['ra','dec'])
        
        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        max_mag_object = xid[np.random.randint(0,len(xid))]
        
        righta_max = max_mag_object["ra"]
        declin_max = max_mag_object["dec"]
        sep = pos.separation(coords.SkyCoord( ra = righta_max*u.degree , dec = declin_max*u.degree))
        sep_list.append(sep.arcsecond)
        
        
    return sep_list

#%%
print("hi")
#a = loop_through__sep(t_data_PLANCK,filter_one,number_of_compact_sources,magnitude_limit,z_lim)

#np.save("Data/separation_test_{}___lengthGalaxies_{}".format(filter_one,len(a)),a)
a = np.load("Data/separation_test_modelMag_r___lengthGalaxies_3652.npy")
print("hi")
b = rand_loop_through__sep(filter_one, len(a))
print("done")

#%%
length = len(a)
plt.figure()
plt.title("separation of galaxies")
plt.xlabel("separation")
plt.ylabel("Count")
plt.hist(a,bins)
plt.hist(np.asarray(b),bins)
plt.savefig("Plots/separation random Historgram_{}___length_{}".format(filter_one,length))

#%%
end = timeit.default_timer()
print("time taken : {}s".format(round(end-start,2)))
