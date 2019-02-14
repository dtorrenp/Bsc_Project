#%%
from astroquery.sdss import SDSS

from astropy.io import fits
from astropy.table import Table
from astropy import coordinates as coords
import astropy.units as u

import timeit
import numpy as np
import matplotlib.pyplot as plt

#%%
"""Open fits file and create "Table" object with which actual data can be accessed"""
start = timeit.default_timer()
##Planck data
hdulist_PLANCK = fits.open('/Users/daniel/Documents/Bsc_Project/COM_PCCS-Catalogs_vPR2/COM_PCCS_857-excluded_R2.01.fits')
hdu_BinTable_PLANCK = hdulist_PLANCK[1]
t_data_PLANCK = Table(hdu_BinTable_PLANCK.data)
colnames_PLANCK = t_data_PLANCK.colnames

def loop_through_planck(table,f):
    
    new_table = Table()
    
    for i in np.arange(len(table)):
        righta = table["RA"][i]
        declin = table["DEC"][i]
        pos = coords.SkyCoord(righta, declin, unit = u.degree )#NOTE SURE THIS IS CORRECT
        
        xid_test = SDSS.query_region(pos, spectro=True, radius = 2.5*u.arcmin)
        print(xid_test)
        xid = SDSS.query_region(pos, spectro=True, radius = 2.5*u.arcmin, photoobj_fields=['ra','dec','u','g','r','i','z'])
        
        flux_list_filter = []
        for v in np.arange(len(xid)):
            flux_f = xid[f][v]
            flux_list_filter.append(flux_f)
            
        max_flux_object = xid[f][flux_list_filter.index(max(flux_list_filter))]
        new_table.add_row(max_flux_object)
        
#    galactic_stars = new_table[new_table['CLASS'] == "STAR"]
#    planck_galaxies = new_table[new_table['CLASS'] == "GALAXY"]
    
    return new_table

#%%

def cal_abs_mag_astro(table,f,f_plus):
    #ugriz  filters
    z_val = np.array(table["z"])
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

#%%
#
##SDSS data
#hdulist_SDSS = fits.open('/Users/daniel/Documents/Bsc_Project/\SDSS data dr14/DR14Q_v4_4.fits')
###hdu_Primary_SDSS = hdulist_SDSS[0]
#hdu_BinTable_SDSS = hdulist_SDSS[1]
###t_primary_SDSS  = Table(hdu_Primary_SDSS.data)
#t_data_SDSS = Table(hdu_BinTable_SDSS.data)
#colnames_SDSS = t_data_SDSS.colnames
#
#
##SUB - SDSS data
#hdulist_SDSS = fits.open('/Users/daniel/Documents/Bsc_Project/SDSS data dr14/subset_spec_2_5.fits')
###hdu_Primary_SDSS = hdulist_SDSS[0]
#hdu_BinTable_SDSS = hdulist_SDSS[1]
###t_primary_SDSS  = Table(hdu_Primary_SDSS.data)
#t_data_SDSS = Table(hdu_BinTable_SDSS.data)
#
#colnames_SDSS = t_data_SDSS.colnames
#
#reduced_table = t_data_SDSS['CLASS','Z','SPECTROFLUX','GroupID','GroupSize']
##create new table with the ojects within each compact source that has the largest flux
## create table of those with no tag
#
#def Select_by_flux(reduced_table,f):
#    mask_reduced_no_tag = reduced_table['GroupID'] < 1
#    reduced_no_tag = reduced_table[mask_reduced_no_tag]#CORRECT
#    #print(len(reduced_no_tag))
#    
#    mask_reduced_groups = reduced_table['GroupID'] > 0
#    reduced_groups = reduced_table[mask_reduced_groups]#CORRECT
#    #print(len(reduced_groups))
#    
#    
#
#    for i in np.arange(1,max(reduced_groups['GroupID'])+1):
#        group_mask  = reduced_groups['GroupID'] == i
#        red_group_i = reduced_groups[group_mask]
#        
#        flux_list_filter = []
#        for v in np.arange(red_group_i['GroupSize'][0]):
#            fluxes = np.asarray(red_group_i['SPECTROFLUX'][v])
#            flux_filter = fluxes[f]
#            flux_list_filter.append(flux_filter)
#        largest_flux = max(flux_list_filter)
#        index_flux = flux_list_filter.index(largest_flux)
#        largest_flux_object = red_group_i[index_flux]
#        reduced_no_tag.add_row(largest_flux_object)
#    
#    galactic_stars = reduced_no_tag[reduced_no_tag['CLASS'] == "STAR"]
#    planck_galaxies = reduced_no_tag[reduced_no_tag['CLASS'] == "GALAXY"]
#    
#    print(planck_galaxies)
#    return planck_galaxies,galactic_stars
#
##%%
#
#def cal_abs_mag(table,f):
#    #ugriz  filters
#    array_table_flux = np.array(table["SPECTROFLUX"])
#    
#    array_flux_f = np.abs(array_table_flux[:,f])
#    
#    #print(array_flux_f[array_flux_f < 0])#WEIRD SHIT WITH NEGATIVE FLUXES
#    
#    array_flux_f_plus = np.abs(array_table_flux[:,f-1])
#    
#    z_val = np.array(table["Z"])
#    
#    filter_m_f = 22.5 - 2.5*np.log10(np.abs(array_flux_f))
#    
#    filter_m_f_plus = 22.5 - 2.5*np.log10(np.abs(array_flux_f_plus))
#    
#    colour = filter_m_f_plus - filter_m_f
##    print(min(4.28e8*z_val))#NEGATIVE VALUES
#    Abs_mag = filter_m_f - 5*np.log10(4.28e8*z_val) + 5
#    
#    return Abs_mag, colour
#
##%%
###ugriz  filters
##    
##planck_galaxies_2,galactic_stars_2 = Select_by_flux(reduced_table, 2)
##a , b = cal_abs_mag(planck_galaxies_2,2)
##
##plt.figure(1)
##plt.title("Colour Magnitude plot")
##plt.xlabel("Absolute magnitude (g)")
##plt.ylabel("Colour (g-r)")
##plt.scatter(a,b,s = 0.45)
##
##planck_galaxies_1, galactic_stars_1 = Select_by_flux(reduced_table, 1)
##c , d = cal_abs_mag(planck_galaxies_1,1)
##plt.figure(2)
##plt.title("Colour Magnitude plot")
##plt.xlabel("Absolute magnitude (u)")
##plt.ylabel("Colour (u-g)")
##plt.scatter(c,d,s = 0.45)
##
##end = timeit.default_timer()
##print(end-start)
#
##plt.figure()
##plt.hist(b)