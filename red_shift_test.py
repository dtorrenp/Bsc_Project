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
from astropy.table import vstack
import open_planck as op
import timeit
import warnings
from astropy.io import fits

#%%
start = timeit.default_timer()

#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()
print("PLANCK OPEN")
#%%
#ugriz  filters
filter_one = "modelMag_r"
dot_size = 0.35
bins = 100
magnitude_limit = 20
z_lim = 0.6
number_of_compact_sources = len(t_data_PLANCK)
width = 10

#%%
def loop_through__redshift(table,f,length,mag_limit,z_lim):
    
    new_table_galaxy = Table()
    
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
        
        mask_gal = (xid["class"] == "GALAXY")
        xid = xid[mask_gal]
        
        if len(xid) == 0:
            continue
        
        xid.sort(f)
        data_mask = xid[f] > mag_limit
        data = xid[data_mask]
        
        if len(data) == 0:
            data = Table(xid[0])
            
        new_table_galaxy = vstack([new_table_galaxy, data])
        
        print(len(new_table_galaxy))
        
#    print(new_table_galaxy)
    return new_table_galaxy

#%%
def rand_loop_through__redshift(table,f,length,width):
    new_table_galaxy = Table()
    
    while len(new_table_galaxy) < length:
        righta = np.random.uniform(-180,180)
        declin = np.random.uniform(-90,90)
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        
        xid = SDSS.query_region(galatic_pos, spectro = True , radius = width*u.arcmin , specobj_fields=['Z'], photoobj_fields=['ra','dec'])
        
        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        new_table_galaxy = vstack([new_table_galaxy, xid])
        
        print(len(new_table_galaxy))
#    print(new_table_galaxy)
    return new_table_galaxy

#%%
hdulist_SDSS = fits.open('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/local code/fifty_thousand_redshift_unlim.fits')
#hdulist_PLANCK = fits.open('C:/Users/DELL/Desktop/COM_PCCS_857_R2.01.fits')
hdu_BinTable_SDSS = hdulist_SDSS[1]
b = Table(hdu_BinTable_SDSS.data)
print("SDSS OPEN")

#%%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    
    #a = loop_through__redshift(t_data_PLANCK,filter_one,number_of_compact_sources,magnitude_limit,z_lim)
    #a.write("Data/redshift_test_{}___lengthGalaxies_{}".format(filter_one,len(a)), format = 'ascii')
    a = Table.read('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/local code/Data/redshift_test_modelMag_r___lengthGalaxies_2737', format = 'ascii')
    print("DATA OPEN")
    
    #b = rand_loop_through__redshift(t_data_PLANCK,filter_one,len(a),width)
    #b.write("Data/redshift_test_rand_{}___lengthGalaxies_{}".format(filter_one,len(b)), format = 'ascii')
    #b = Table.read('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/local code/sdss_galaxies_data_petrosian_petroMag_r___lengthGalaxies_1716', format = 'ascii')
    #print("RAND DATA OPEN")
    
#%%
length = len(a)

plt.figure()
plt.title("redshift of galaxies")
plt.xlabel("redshift")
plt.ylabel("Count")
plt.hist(b["redshift"],bins)
plt.hist(a["Z"],bins)
plt.savefig("Plots/redshift random Historgram_{}___length_{}_bins_{}".format(filter_one,length,bins))

plt.figure()
plt.title("redshift of galaxies - planck")
plt.xlabel("redshift")
plt.ylabel("Count")
plt.hist(a["Z"],bins)
plt.savefig("Plots/redshift random Historgram_ planck_{}___length_{}_bins_{}".format(filter_one,length,bins))

plt.figure()
plt.title("redshift of galaxies - random")
plt.xlabel("redshift")
plt.ylabel("Count")
plt.hist(b["redshift"],bins)
plt.savefig("Plots/redshift random Historgram_random_{}___length_{}_bins_{}".format(filter_one,length,bins))

plt.figure()
plt.title("Normalised redshift of galaxies")
plt.xlabel("redshift")
plt.ylabel("Count")
plt.hist(b["redshift"],bins,density = True)
plt.hist(a["Z"],bins, density = True)
plt.savefig("Plots/Normalised redshift random Historgram_{}___length_{}_bins_{}".format(filter_one,length,bins))

print("PLOTS DONE")

#%%
end = timeit.default_timer()
print("time taken : {}s".format(round(end-start,2)))