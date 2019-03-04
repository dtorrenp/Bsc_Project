# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:33:22 2019

@author: daniel
"""
#%%
from astroquery.ned import Ned
import astropy.units as u
from astropy import coordinates

import open_planck as op
#import produce_galaxies_stars_other as pg
#import abs_mag as ab
#import produce_plots as pp
#import top_ten_plate_pics as tt
import warnings
#import timeit

#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    righta = t_data_PLANCK["RA"][0]
    declin = t_data_PLANCK["DEC"][0]
            
    pos = coordinates.SkyCoord(ra=righta, dec=declin, unit=(u.deg, u.deg), frame='fk4')
    galatic_pos = pos.galactic
    result_table = Ned.query_region(galatic_pos, radius = 2.5*u.arcmin)
    #result_table = Ned.get_table("3C 273", table='positions')
    
    print(result_table)
