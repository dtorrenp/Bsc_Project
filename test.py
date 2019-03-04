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

def loop_through_planck(table,f,length,mag_limit,z_lim):
    
    new_table_galaxy = Table(names=(table.colnames))
    new_table_galaxy.add_column(name = 'Aperture Type')
    
    new_table_star = Table(names=(table.colnames))
    new_table_other = Table(names=(table.colnames))
    
    sep_list = []
        
    mag_list = []
    mag_list_diff_one_two = []
    mag_list_diff_two_three = []
    
    for i in np.arange(length):
        
        righta = table["RA"][i]
        declin = table["DEC"][i]
        
        if np.abs(table["GLAT"][i]) < 20:
            continue
        
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        
        #objID not working correctly gives error "invalid load key, '\x00'." 
        #print("HI")
        #xid = SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z','class'], photoobj_fields=['ra','dec','u','g','r','i'])
        xid =  SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , fields=['Z','objID','class','ra','dec','u','g','r','i','deVlnL','explnL'])
        
        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        
        xid.sort(f)
        
        sub_length = len(xid)
        
        mag_list = np.concatenate((mag_list,xid[f]))
        
#        if len(xid) == 1:
#            mag_list_diff_one_two.append(xid[f][0])
        if len(xid) == 2:
            mag_list_diff_one_two.append(xid[f][1] - xid[f][0])
        elif len(xid) > 2:
            mag_list_diff_one_two.append(xid[f][1] - xid[f][0])
            mag_list_diff_two_three.append(xid[f][2] - xid[f][1])
            
        print("hi")
        mask_redshift = xid['Z'] < z_lim
        print(mask_redshift)
        xid = xid[mask_redshift]
        print(xid)
        if xid is None:
            continue
        
        for v in np.arange(sub_length):
            
            max_mag_object = xid[v]
            
            if max_mag_object[f] > mag_limit and v > 0:
                break
            
            if max_mag_object["class"] == "GALAXY":
                new_table_galaxy = vstack([new_table_galaxy, max_mag_object])
                righta_max = max_mag_object["ra"]
                declin_max = max_mag_object["dec"]
                sep = pos.separation(coords.SkyCoord( ra = righta_max*u.degree , dec = declin_max*u.degree))
                sep_list.append(sep.arcsecond)
                    
            elif max_mag_object["class"] == "STAR":
                new_table_star = vstack([new_table_star, max_mag_object])
            else:
                new_table_other = vstack([new_table_other, max_mag_object])
                
    return [new_table_galaxy,new_table_star,new_table_other,sep_list,mag_list,mag_list_diff_one_two,mag_list_diff_two_three]
