#%%
from astropy.io import fits
from astropy.table import Table
#%%
"""Open fits file and create "Table" object with which actual data can be accessed"""

hdulist = fits.open('/Users/daniel/Documents/Bsc_Project/COM_PCCS-Catalogs_vPR2/COM_PCCS_857-excluded_R2.01.fits')
                 
print("FITS FILE INFO")
hdulist.info()
hdu_Primary = hdulist[0]
hdu_BinTable = hdulist[1]
#print("PRIMARY SHAPE")
#print(hdu_Primary.data)
#print("PRIMARY HEADER")
#print(hdu_Primary.header)
#print("BinTable SHAPE")
#print(hdu_BinTable.data.shape)
#print("BinTable HEADER")
#print(hdu_BinTable.header)

print("Table Data")
t_data = Table(hdu_BinTable.data)
print(t_data)
#%%

