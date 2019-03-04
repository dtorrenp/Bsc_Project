# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 11:41:37 2019

@author: daniel
"""
#%%

from astroquery.sdss import SDSS
from astropy.table import Table, unique
from astropy import coordinates as coords
from astropy.table import vstack
import astropy.units as u
import numpy as np

#%%

def loop_through_planck_petro(table,f,length,mag_limit,z_lim):
    
    new_table_galaxy = Table(names=(table.colnames))
    
    for i in np.arange(length):
        
        righta = table["RA"][i]
        declin = table["DEC"][i]
        
        if np.abs(table["GLAT"][i]) < 20:
            continue
        
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        
        xid = SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z','class'], photoobj_fields=['ra','dec','petroMag_u','petroMag_g','petroMag_r','petroMag_i'])
        
        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        
        mask_redshift = xid['Z'] < z_lim
        xid = xid[mask_redshift]
        
        if xid is None:
            continue
        
        xid.sort(f)
        
        for v in np.arange(len(xid)):
            
            max_mag_object = xid[v]
            max_mag_object = Table(max_mag_object)
            
            if max_mag_object[f] > mag_limit and v > 0:
                break
            
            if max_mag_object["class"] == "GALAXY":
                new_table_galaxy = vstack([new_table_galaxy, max_mag_object])
                    
    return [new_table_galaxy]

