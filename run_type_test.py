# -*- coding: utf-8 -*-
#%%
import open_planck as op
import abs_mag as ab
import warnings
import timeit
import numpy as np
from astroquery.sdss import SDSS
from astropy.table import Column
from astropy.table import Table, unique
from astropy import coordinates as coords
from astropy.table import vstack
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from astropy.io import fits

#%%
start = timeit.default_timer()

#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
#ugriz  filters
filter_one = "modelMag_g"
filter_two = "modelMag_r"
dot_size = 0.35
bins = 100
magnitude_limit = 20
z_lim = 0.6
number_of_compact_sources = len(t_data_PLANCK)

#%%
t_data_PLANCK,colnames_PLANCK = op.open_planck_data()

#%%
def loop_through_planck(table,f,length,mag_limit,z_lim):
    
    new_table_galaxy = Table()
    col = Column(name='Aperture Type', length = len(new_table_galaxy))
    new_table_galaxy.add_column(col)
    
    for i in np.arange(length):
        
        righta = table["RA"][i]
        declin = table["DEC"][i]
        
        if np.abs(table["GLAT"][i]) < 20:
            continue
        
        pos = coords.SkyCoord( ra = righta*u.degree , dec = declin*u.degree)
        galatic_pos = pos.galactic
        
        xid = SDSS.query_region(galatic_pos, spectro = True , radius = 2.5*u.arcmin , specobj_fields=['Z','class', 'BestObjID'], photoobj_fields=['ra','dec','modelMag_u','modelMag_g','modelMag_r','modelMag_i','lnLExp_u','lnLDeV_u'])
        
        if xid is None:
            continue
        
        xid = unique(xid,keys=['ra', 'dec'])
        
        xid.sort(f)
            
        mask_redshift_gal = (xid['Z'] < z_lim) & (xid["class"] == "GALAXY")
        xid = xid[mask_redshift_gal]
        
        if len(xid) == 0:
            continue
        
        xid.sort(f)
        
        data_mask = xid[f] > mag_limit
        data = xid[data_mask]
        
        if len(data) == 0:
            data = Table(xid[0])
#        else:
#            new_table_galaxy = vstack([new_table_galaxy, data])
        
        deV_data_mask = data['lnLDeV_u'] > data['lnLExp_u']
        exp_data_mask = np.invert(deV_data_mask)
        
        deV_data = data[deV_data_mask]
        exp_data = data[exp_data_mask]
        
        if len(deV_data) != 0:
            deV_data['Aperture Type'] = 0
            new_table_galaxy = vstack([new_table_galaxy, deV_data])
            
        if len(exp_data) != 0:
            exp_data['Aperture Type'] = 1
            new_table_galaxy = vstack([new_table_galaxy, exp_data])
        print(len(new_table_galaxy))
    return [new_table_galaxy]

#%%
    

#hdulist_PLANCK = fits.open('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/COM_PCCS_857_R2.01.fits')
#hdulist_SDSS = fits.open('/Users/daniel/Documents/Bsc_Project/random_galaxies_model_mag.fits')
hdulist_SDSS = fits.open('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/local code/fifty_thousand_redshift_lim.fits')
#hdulist_PLANCK = fits.open('C:/Users/DELL/Desktop/COM_PCCS_857_R2.01.fits')
hdu_BinTable_SDSS = hdulist_SDSS[1]
t_data_SDSS = Table(hdu_BinTable_SDSS.data)
print("SDSS OPEN")

#%%
#ugriz  filters
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    
    #SDSS_DATA = loop_through_planck(t_data_PLANCK,filter_one,number_of_compact_sources,magnitude_limit,z_lim)
    #length = len(SDSS_DATA[0])
    #SDSS_DATA[0].write("Data/sdss_galaxies_data_type_test_{}___lengthGalaxies_{}".format(filter_one,length), format = 'ascii')
    
    SDSS_DATA = [Table.read('//icnas2.cc.ic.ac.uk/dt1716/Year 3/Bsc Project/local code/Data/sdss_galaxies_data_type_test_modelMag_g___lengthGalaxies_2839', format = 'ascii')]
    length = len(SDSS_DATA[0])
    print("PLANCK CROSS MATCH")
    
#%%
    
galaxies = SDSS_DATA[0]
exp_mask = galaxies["Aperture Type"] == 1
de_v_mask = np.invert(exp_mask)

exp_galaxies = galaxies[exp_mask]
de_v_galaxies = galaxies[de_v_mask]

exp_table = ab.cal_abs_mag_astro(exp_galaxies,filter_one,filter_two)
de_v_table = ab.cal_abs_mag_astro(de_v_galaxies,filter_one,filter_two)
#data_table = ab.cal_abs_mag_astro(galaxies,filter_one,filter_two)

SDSS_sample =  ab.cal_abs_mag_astro_sdss(t_data_SDSS ,filter_one,filter_two)

print("ABSOLUTE AND COLOUR CALC")

#%%
#CONVERT TO DATAFRAME
df_exp = pd.DataFrame(columns=['colour', 'abs mag'])
df_exp['abs mag'] = exp_table['Absolute Magnitude {}'.format(filter_one)]
df_exp['colour'] = exp_table['Colour ({} - {})'.format(filter_two,filter_one)]

df_de_v = pd.DataFrame(columns=['colour', 'abs mag'])
df_de_v['abs mag'] = de_v_table['Absolute Magnitude {}'.format(filter_one)]
df_de_v['colour'] = de_v_table['Colour ({} - {})'.format(filter_two,filter_one)]

df_sdss = pd.DataFrame(columns=['colour', 'abs mag'])
df_sdss['abs mag'] = SDSS_sample['Absolute Magnitude {}'.format(filter_one)]
df_sdss['colour'] = SDSS_sample['Colour ({} - {})'.format(filter_two,filter_one)]

print("DATA FRAME COMPLETE")

#%%

def produce_colour_magnitude_ap_type(exp_table,de_v_table,filter_one,filter_two,dot_size):
    plt.figure()
    plt.title("Colour Magnitude plot aperture type")
    plt.xlabel("Absolute magnitude ({})".format(filter_one))
    plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
    plt.scatter(df_exp['abs mag'],df_exp['colour'],s = dot_size, label = "Exp")
    plt.scatter(df_de_v['abs mag'],df_de_v['colour'],s = dot_size, label = "De V")
    plt.legend()
    plt.savefig("Plots/Colour_Magnitude_plot_aperture_type_{}_{}___length_exp_{}_length_de_V_{}.png".format(filter_two,filter_one,len(exp_table),len(de_v_table)))
    
    plt.figure()
    plt.title("sdss sample colour contours Histogram")
    plt.xlabel("Colour ({}-{})".format(filter_one,filter_two))
    plt.ylabel("Count")
    plt.hist(df_sdss['colour'], bins = bins,alpha = 0.5)
    plt.savefig("Plots/Colour_Magnitude_sdss_sample-Histogram_test_contours_{}_{}___length_{}.png".format(filter_one,filter_two,length))
    
    plt.figure()
    plt.title("Colour Magnitude plot aperture type")
    plt.xlabel("Absolute magnitude ({})".format(filter_one))
    plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
    ax_one = sns.kdeplot(df_sdss['abs mag'], df_sdss['colour'], cmap='cubehelix', shade=False, shade_lowest= False, cbar=True, n_levels = 8)
    ax_one.scatter(df_exp['abs mag'],df_exp['colour'],s = dot_size, label = "Exp")
    ax_one.scatter(df_de_v['abs mag'],df_de_v['colour'],s = dot_size, label = "De V")
    plt.legend()
    plt.savefig("Plots/Colour_Magnitude_plot_aperture_SDSS_sample_type_{}_{}___length_exp_{}_length_de_V_{}.png".format(filter_two,filter_one,len(exp_table),len(de_v_table)))
    
    plt.figure()
    plt.title("Type galaxy contours")
    plt.xlabel("Absolute magnitude ({})".format(filter_one))
    plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
    ax_one = sns.kdeplot(df_exp['abs mag'], df_exp['colour'], cmap='Reds', shade=False, shade_lowest= False, cbar=True, n_levels = 15 )
    ax_two = sns.kdeplot(df_de_v['abs mag'], df_de_v['colour'], cmap='Greens', shade=False, shade_lowest= False, cbar=True, n_levels = 15 )
    ax_one.scatter(df_exp['abs mag'],df_exp['colour'], s = dot_size, label = "Exp")
    ax_two.scatter(df_de_v['abs mag'],df_de_v['colour'], s = dot_size, label = "De V")
    plt.legend()
    plt.savefig("Plots/Colour_Magnitude_type_test_contours_{}_{}___length_{}.png".format(filter_one,filter_two,length))
    
    plt.figure()
    plt.title("Type galaxy contours")
    plt.xlabel("Absolute magnitude ({})".format(filter_one))
    plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
    ax_one = sns.kdeplot(df_exp['abs mag'], df_exp['colour'], cmap='Reds', shade=False, shade_lowest= False, cbar=True, n_levels = 15 )
    ax_one.scatter(df_exp['abs mag'],df_exp['colour'], s = dot_size, label = "Exp")
    plt.legend()
    plt.savefig("Plots/Colour_Magnitude_type__exp_test_contours_{}_{}___length_{}.png".format(filter_one,filter_two,length))
    
    plt.figure()
    plt.title("Type galaxy contours")
    plt.xlabel("Absolute magnitude ({})".format(filter_one))
    plt.ylabel("Colour ({}-{})".format(filter_one,filter_two))
    ax_two = sns.kdeplot(df_de_v['abs mag'], df_de_v['colour'], cmap='Greens', shade=False, shade_lowest= False, cbar=True, n_levels = 15 )
    ax_two.scatter(df_de_v['abs mag'],df_de_v['colour'], s = dot_size, label = "De V")
    plt.legend()
    plt.savefig("Plots/Colour_Magnitude_type_De_V_test_contours_{}_{}___length_{}.png".format(filter_one,filter_two,length))
    
    plt.figure()
    plt.title("Type galaxy contours Histogram")
    plt.xlabel("Colour ({}-{})".format(filter_one,filter_two))
    plt.ylabel("Count")
    plt.hist(df_exp['colour'], bins = bins, alpha = 0.5, label = "Exp")
    plt.hist(df_de_v['colour'], bins = bins,alpha = 0.5, label = "De V")
    plt.legend()
    plt.savefig("Plots/Colour_Magnitude_type-Histogram_test_contours_{}_{}___length_{}.png".format(filter_one,filter_two,length))
    
    plt.figure()
    plt.title("Type galaxy contours Histogram")
    plt.xlabel("Colour ({}-{})".format(filter_one,filter_two))
    plt.ylabel("Count")
    plt.hist(df_exp['colour'], bins = bins, alpha = 0.5, label = "Exp")
    plt.hist(df_de_v['colour'], bins = bins,alpha = 0.5, label = "De V")
    plt.hist(df_sdss['colour'], bins = bins,alpha = 0.5, label = "SDSS sample")
    plt.legend()
    plt.savefig("Plots/Colour_Magnitude_type_sdss_toooo____Histogram_test_contours_{}_{}___length_{}.png".format(filter_one,filter_two,length))
    
#%%
produce_colour_magnitude_ap_type(exp_table,de_v_table,filter_one,filter_two,dot_size)

#%%
end = timeit.default_timer()
print("time taken : {}s".format(round(end-start,2)))