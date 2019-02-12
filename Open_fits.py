#%%
from astropy.io import fits
from astropy.table import Table
import timeit
#%%
"""Open fits file and create "Table" object with which actual data can be accessed"""
start = timeit.default_timer()
#Planck data
#hdulist_PLANCK = fits.open('/Users/daniel/Documents/Bsc_Project/COM_PCCS-Catalogs_vPR2/COM_PCCS_857-excluded_R2.01.fits')
#hdu_Primary_PLANCK = hdulist_PLANCK[0]
#hdu_BinTable_PLANCK = hdulist_PLANCK[1]
#t_primary_PLANCK  = Table(hdu_Primary_PLANCK.data)
#t_data_PLANCK = Table(hdu_BinTable_PLANCK.data)
#colnames_PLANCK = t_data_PLANCK.colnames
print("hi")
#SDSS data
hdulist_SDSS = fits.open('/Users/daniel/Documents/Bsc_Project/\SDSS data dr14/specObj-dr14.fits')
print("hi")
hdu_Primary_SDSS = hdulist_SDSS[0]
print("hi")
hdu_BinTable_SDSS = hdulist_SDSS[1]
print("hi")
t_primary_SDSS  = Table(hdu_Primary_SDSS.data)
print("hi")
t_data_SDSS = Table(hdu_BinTable_SDSS.data)
print("hi")
colnames_SDSS = t_data_SDSS.colnames
print("hi")
end = timeit.default_timer()
print(end-start)
#%%
def subset_table(table,column,condition):
    mask = table[column] == condition
    print(len(mask))
    subset = table[mask]
    #print(subset)
    return subset
#%%
#extended_sub_table = subset_table(t_data,"EXTENDED",0)#subset with only extened flag equal to 0
#one_extended_sub_table = subset_table(t_data,"EXTENDED",1)
#
#extended_external_sub_sub_table = subset_table(extended_sub_table,"EXT_VAL",3)#subset within reduced table with ext_val equal to 3
