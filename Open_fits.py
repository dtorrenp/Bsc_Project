#%%
from astropy.io import fits
from astropy.table import Table
import timeit
import numpy as np
import matplotlib.pyplot as plt
#%%
"""Open fits file and create "Table" object with which actual data can be accessed"""
#start = timeit.default_timer()
##Planck data
#hdulist_PLANCK = fits.open('/Users/daniel/Documents/Bsc_Project/COM_PCCS-Catalogs_vPR2/COM_PCCS_857-excluded_R2.01.fits')
##hdu_Primary_PLANCK = hdulist_PLANCK[0]
#hdu_BinTable_PLANCK = hdulist_PLANCK[1]
##t_primary_PLANCK  = Table(hdu_Primary_PLANCK.data)
#t_data_PLANCK = Table(hdu_BinTable_PLANCK.data)
#colnames_PLANCK = t_data_PLANCK.colnames
#
##SDSS data
#hdulist_SDSS = fits.open('/Users/daniel/Documents/Bsc_Project/\SDSS data dr14/DR14Q_v4_4.fits')
###hdu_Primary_SDSS = hdulist_SDSS[0]
#hdu_BinTable_SDSS = hdulist_SDSS[1]
###t_primary_SDSS  = Table(hdu_Primary_SDSS.data)
#t_data_SDSS = Table(hdu_BinTable_SDSS.data)
#colnames_SDSS = t_data_SDSS.colnames
#

#SUB - SDSS data
hdulist_SDSS = fits.open('/Users/daniel/Documents/Bsc_Project/SDSS data dr14/subset_spec_2_5.fits')
##hdu_Primary_SDSS = hdulist_SDSS[0]
hdu_BinTable_SDSS = hdulist_SDSS[1]
##t_primary_SDSS  = Table(hdu_Primary_SDSS.data)
t_data_SDSS = Table(hdu_BinTable_SDSS.data)
colnames_SDSS = t_data_SDSS.colnames

#end = timeit.default_timer()
#print(end-start)



#%%
#def select_largest_mag(table,f):
#    
#    
#    
#    
#    array_table_flux = np.array(table["SPECTROFLUX"])
#    
#    
#    
#    return sub_table

def cal_abs_mag(table,f):
    #ugriz  filters
    array_table_flux = np.array(table["SPECTROFLUX"])
    
    array_flux_f = np.abs(array_table_flux[:,f])
#    print("array_flux_f")
#    print(array_flux_f)
    #print(array_flux_f[array_flux_f < 0])#WEIRD SHIT WITH NEGATIVE FLUXES
    
    array_flux_f_plus = np.abs(array_table_flux[:,f-1])
#    print("array_flux_f_plus")
#    print(array_flux_f_plus)
    
    z_val = np.array(table["Z"])
#    print("z_val")
#    print(z_val)
    
    filter_m_f = 22.5 - 2.5*np.log10(np.abs(array_flux_f))
    
    filter_m_f_plus = 22.5 - 2.5*np.log10(np.abs(array_flux_f_plus))
#    print("filter_m")
#    print(filter_m_f)
#    print("filter_m_plus")
#    print(filter_m_f_plus)
    
    colour = filter_m_f_plus - filter_m_f
#    print("colour")
#    print(colour)
    
    Abs_mag = filter_m_f - 5*np.log10(4.28e8*z_val) + 5
#    print("Abs_mag")
#    print(Abs_mag)
    
    return Abs_mag, colour
#%%
#ugriz  filters
a , b = cal_abs_mag(t_data_SDSS,2)
plt.title("Colour Magnitude plot")
plt.xlabel("Absolute magnitude(g)")
plt.ylabel("Colour (g-r)")
plt.scatter(a,b,s = 0.3)
#
#plt.figure()
#plt.hist(b)

#extended_sub_table = subset_table(t_data,"EXTENDED",0)#subset with only extened flag equal to 0
#one_extended_sub_table = subset_table(t_data,"EXTENDED",1)
#
#extended_external_sub_sub_table = subset_table(extended_sub_table,"EXT_VAL",3)#subset within reduced table with ext_val equal to 3
