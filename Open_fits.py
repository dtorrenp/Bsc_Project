#%%
from astropy.io import fits
from astropy.table import Table
#%%
"""Open fits file and create "Table" object with which actual data can be accessed"""

hdulist = fits.open('/Users/daniel/Documents/Bsc_Project/COM_PCCS-Catalogs_vPR2/COM_PCCS_857-excluded_R2.01.fits')
hdu_Primary = hdulist[0]
hdu_BinTable = hdulist[1]

t_data = Table(hdu_BinTable.data)
colnames = t_data.colnames
#%%
def subset_table(table,column,condition):
    subset = table[table[column] == condition]
    print(subset)
    return subset
#%%
extended_sub_table = subset_table(t_data,"EXTENDED",0)
extended_external_sub_sub_table = subset_table(extended_sub_table,"EXT_VAL",3)

