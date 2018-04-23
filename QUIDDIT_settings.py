# -*- coding: utf-8 -*-
"""
Created on Mon Sep 04 16:10:39 2017
This file contains a number of settings used across QUIDDIT
@author: ls13943
"""
import os
import numpy as np
import matplotlib as mpl

home = os.getcwd()
version='2.0'

# path to type IIa spectum (CSV)
#IIa_path = 'C:\FTIR/typeIIa.csv'
IIa_path = os.getcwd() + '/typeIIa.csv'


#path to file with standard spectra of N components (CSV)
std_path = os.getcwd() + '/CAXBD.csv'
std = np.loadtxt(std_path, delimiter = ',')     # read CAXBD spectra

#standard first guess for platelet fit (p_x0, p_I, p_HWHM_l, p_HWHM_r, p_sigma, 
#H1405_x0, H1405_I, H1405_HWHM_l, H1405_HWHM_r, H1405_sigma, 
#B_x0, B_I, B_HWHM_l, B_HWHM_r, B_sigma,
#const)
pp_res_prev = (1370, 0, 5, 5, 1, 
               1405, 0, 5, 5, 1, 
               1332, 0, 5, 5, 0, 
               1)

N_comp = np.array((0, 1, 0, 1, 1, 1))

#pp_res_prev = (1365, 1.2, 3, 3, 1, 
#               1405, 0, 5, 5, 1, 
#               1332, 0, 5, 5, 0, 
#              1)

#*this set of parameters is only used when the first attempt of fitting fails.


# settings for plotting

colors=['blue','green','red']
levels=[0,1]

cmap, norm = mpl.colors.from_levels_and_colors([1,63,76,79], ['blue','green','red'])

f=16
l=2
m=3