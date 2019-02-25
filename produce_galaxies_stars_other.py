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

def loop_through_planck(table,f):
    
    length = len(table)
    #length = 100
    
    new_table_galaxy = Table(names=(table.colnames))
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
        xid = SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z','class'], photoobj_fields=['ra','dec','u','g','r','i'])

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

        if sub_length == 1:
            number_objects_to_look_at = 1
        elif sub_length == 2:
            number_objects_to_look_at = 1
        else:
            number_objects_to_look_at = 1
        
        for v in np.arange(number_objects_to_look_at):
            max_mag_object = xid[v]
            
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