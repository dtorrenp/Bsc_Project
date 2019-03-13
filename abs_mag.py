# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 11:41:38 2019

@author: daniel
"""
#%%

import numpy as np
from astropy.table import Column
#%%
def cal_abs_mag_astro(table,f,f_plus):
    #ugriz  filters
    red_mask = table["Z"] > 0
    red_table = table[red_mask]
    z_val = np.array(red_table["Z"])
    
    colour = np.asarray(red_table[f_plus]) - np.asarray(red_table[f])
    Abs_mag = np.asarray(red_table[f]) - 5*np.log10(4.28e8*z_val) + 5 #NO PLUS FIVE
    
    colour_col = Column(data = colour,name='Colour ({} - {})'.format(f_plus,f))
    abs_col = Column(data = Abs_mag, name='Absolute Magnitude {}'.format(f))
    
    red_table.add_column(colour_col)
    red_table.add_column(abs_col)
    
    return red_table

def cal_abs_mag_astro_sdss(table,f,f_plus):
    #ugriz  filters
    red_mask = table["redshift"] > 0
    red_table = table[red_mask]
    z_val = np.array(red_table["redshift"])
    colour = np.asarray(red_table[f_plus]) - np.asarray(red_table[f])
    Abs_mag = np.asarray(red_table[f]) - 5*np.log10(4.28e8*z_val) + 5
    colour_col = Column(data = colour, name='Colour ({} - {})'.format(f_plus,f))
    abs_col = Column(data = Abs_mag, name='Absolute Magnitude {}'.format(f))
    red_table.add_column(colour_col)
    red_table.add_column(abs_col)#NO PLUS FIVE

    return red_table