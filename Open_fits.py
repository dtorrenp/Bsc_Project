#%%
from astroquery.sdss import SDSS

from astropy.io import fits
from astropy.table import Table
from astropy import coordinates as coords
from astropy.table import vstack
import astropy.units as u

import timeit
import numpy as np
import matplotlib.pyplot as plt

#%%
"""Open fits file and create "Table" object with which actual data can be accessed"""
start = timeit.default_timer()
##Planck data
hdulist_PLANCK = fits.open('/Users/daniel/Documents/Bsc_Project/COM_PCCS-Catalogs_vPR2/COM_PCCS_857_R2.01.fits')
hdu_BinTable_PLANCK = hdulist_PLANCK[1]
t_data_PLANCK = Table(hdu_BinTable_PLANCK.data)
colnames_PLANCK = t_data_PLANCK.colnames

def loop_through_planck(table,f):
    length = len(table)
    new_table = Table(names=(table.colnames))
    
    for i in np.arange(length):
        print(str(np.around(i/length,3)) + "%")
        righta = table["RA"][i]
        declin = table["DEC"][i]
        
        if np.abs(table["GLAT"][i]) < 20:
#            print("less than 20")
            continue
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        
        xid = SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z'], photoobj_fields=['u','g','r','i'])
        
        if xid is None:
#            print("NOTHING THERE")
            continue
        
        flux_list_filter = []
        
        for v in np.arange(len(xid)):
            flux_f = xid[f][v]
            flux_list_filter.append(flux_f)
            
        position = flux_list_filter.index(max(flux_list_filter))
        max_flux_object = xid[position]
        new_table = vstack([new_table, max_flux_object])
    
    return new_table
#%%

def cal_abs_mag_astro(table,f,f_plus):
    #ugriz  filters
    z_val = np.array(table["Z"])
    colour = np.asarray(table[f_plus]) - np.asarray(table[f])
    Abs_mag = np.asarray(table[f]) - 5*np.log10(4.28e8*z_val) + 5
    return Abs_mag, colour

#%%
#ugriz  filters
    
sdss_galaxies_g = loop_through_planck(t_data_PLANCK, "g")
a , b = cal_abs_mag_astro(sdss_galaxies_g,"g","r")
plt.figure(1)
plt.title("Colour Magnitude plot")
plt.xlabel("Absolute magnitude (g)")
plt.ylabel("Colour (g-r)")
plt.scatter(a,b,s = 0.45)

sdss_galaxies_u = loop_through_planck(t_data_PLANCK, "u")
c , d = cal_abs_mag_astro(sdss_galaxies_u,"u","g")
plt.figure(2)
plt.title("Colour Magnitude plot")
plt.xlabel("Absolute magnitude (u)")
plt.ylabel("Colour (u-g)")
plt.scatter(c,d,s = 0.45)

end = timeit.default_timer()
print(end-start)