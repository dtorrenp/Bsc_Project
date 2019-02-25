# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 15:43:39 2019

@author: daniel
"""
#%%
from astropy import coordinates as coords
from astroquery.sdss import SDSS
import astropy.units as u
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
#%%
plt.ioff()#TURNS OFF PICTURES PRODUCED IN CONSOLE, LOOK IN YOUR FOLDER
#%%
def top_ten_pics(galaxies,filter_one):
    galaxies.sort(filter_one)
    top_ten = galaxies[:9]
    top_ten_pos = coords.SkyCoord( ra = top_ten['ra']*u.degree ,  dec =top_ten['dec']*u.degree)
    data_list_img = []
    
    for i in np.arange(len(top_ten)):
        imgs = SDSS.get_images(top_ten_pos[i])
        
        for v in np.arange(len(imgs)):
            a = fits.HDUList(imgs[v])
            b = a[0].data
            
            data_list_img.append(b)
            
            plt.figure()
            plt.title("Plate image _(%s); number (%s), plate (%s)"%(filter_one,i,v))
            plt.imshow(data_list_img[i])
            plt.savefig("Plate_Pics/Plate image _(%s); number (%s), plate (%s)"%(filter_one,i,v))
    return data_list_img