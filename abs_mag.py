# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 11:41:38 2019

@author: daniel
"""
#%%

import numpy as np

#%%
def cal_abs_mag_astro(table,f,f_plus):
    #ugriz  filters
    red_mask = table["Z"] > 0
    red_table = table[red_mask]
    z_val = np.array(red_table["Z"])
    colour = np.asarray(red_table[f_plus]) - np.asarray(red_table[f])
    Abs_mag = np.asarray(red_table[f]) - 5*np.log10(4.28e8*z_val) + 5
    
    # ABSOLUTE MAGNITUDE NOT FEASBLE
    #wathc out for infinity magnitude!!!!!
    return Abs_mag, colour