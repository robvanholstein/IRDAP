#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This file contains utilities to get the star flux density in Jy.

Author: Julien Milli
'''     

from astropy.io import ascii,fits
import numpy as np
import os #,sys
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy import constants as const
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
path_static_calib_dir = os.path.join(os.path.dirname(__file__), 'static_calibs')

###############################################################################
# 
###############################################################################


def sphere_transmission(BB_filter='B_H', NDset=0.):
    """
    
    Input:
        - BB_filter: name of the broad band filter (among 'B_H',
                    'B_Hcmp2', 'B_J', 'B_Ks', 'B_ND-H', 'B_Y'). By default, assumes 'B_H'
        - NDset: ND filter (float) among 0., 1., 2. or 3.5. By default 0.
    """
    # BB filter
    data_bb = ascii.read(os.path.join(path_static_calib_dir,'SPHERE_IRDIS_'+BB_filter+'.txt'))
    w_bb = data_bb['col1']
    t_bb = data_bb['col2']
    # ND CPI
    data_nd = ascii.read(os.path.join(path_static_calib_dir,'SPHERE_CPI_ND.txt'))
    w_nd = data_nd['col1']
    if float(NDset) == 0.:
        t_nd = data_nd['col2']
    elif float(NDset) == 1.:
        t_nd = data_nd['col3']
    elif float(NDset) == 2.:
        t_nd = data_nd['col4']
    elif float(NDset) == 3.5:
        t_nd = data_nd['col5']
    else:
        print('I did not understand your choice of ND filter: {0:3.2f}'.format(NDset))
        return
    # interpolation
    lambdainterp  = np.arange(900,2401,1)

    interp_function_bb = interp1d(w_bb,t_bb)
    t_bb_interp = interp_function_bb(lambdainterp)
    t_bb_interp[t_bb_interp<0.]=0.

    interp_function_nd = interp1d(w_nd,t_nd)
    t_nd_interp = interp_function_nd(lambdainterp)
    t_nd_interp[t_nd_interp<0.]=0.

    t_final = np.sum(t_bb_interp)
    t_final_nd = np.sum(t_bb_interp *t_nd_interp)
    
    return np.asarray(t_final_nd / t_final)
