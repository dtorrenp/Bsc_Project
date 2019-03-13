# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 11:41:36 2019

@author: daniel
"""
#%%
from astropy.io import fits
from astropy.table import Table
#%%
##Planck data
def open_planck_data():
    
    hdulist_PLANCK = fits.open('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/COM_PCCS_857_R2.01.fits')
    #hdulist_PLANCK = fits.open('/Users/daniel/Documents/Bsc_Project/COM_PCCS-Catalogs_vPR2/COM_PCCS_857_R2.01.fits')
    #hdulist_PLANCK = fits.open('C:/Users/DELL/Desktop/COM_PCCS_857_R2.01.fits')
    hdu_BinTable_PLANCK = hdulist_PLANCK[1]
    t_data_PLANCK = Table(hdu_BinTable_PLANCK.data)
    colnames_PLANCK = t_data_PLANCK.colnames
    
    return t_data_PLANCK,colnames_PLANCK 