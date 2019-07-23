'''
This file contains all functions used by IRDAP.

IRDAP is a Python package to accurately reduce SPHERE-IRDIS polarimetric data.
Copyright (C) 2019 R.G. van Holstein

Full documentation: https://irdap.readthedocs.io
Feedback, questions, comments: vanholstein@strw.leidenuniv.nl

When publishing data reduced with IRDAP, please cite van Holstein et al. 
(2019): <ADS link>. 
For data in pupil-tracking mode please additionally cite van Holstein et al. 
(2017): https://ui.adsabs.harvard.edu/abs/2017SPIE10400E..15V.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''   

###############################################################################
# Import packages
###############################################################################

import os
import sys
import glob
import time
import datetime
import warnings
import shutil
import configparser
import textwrap
import urllib
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from ast import literal_eval
from scipy.ndimage.interpolation import rotate
from scipy import interpolate
from scipy import optimize
from scipy import ndimage
from scipy.interpolate import interp1d
from scipy.stats import trim_mean
from scipy.stats import sigmaclip
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
from skimage.transform import rotate as rotateskimage
from skimage.feature import register_translation
from .version import __version__

###############################################################################
# read_config_file
###############################################################################

def read_config_file(path_config_file):
    '''
    Read the configuration file with input parameters
    
    Input:
        path_config_file: string specifying path of configuration file
            
    Output:
        all input parameters from the configuration file
    
    File written by Rob van Holstein
    Function status: verified       
    '''
    
    def config_true_false(x):
        if x in ['True', 'False']:
            return literal_eval(x)
        else:
            return x

    def config_list_tuple(x):
        if '(' in x or '[' in x:
            return literal_eval(x)
        else:
            return x

    # Create a configparser object
    config = configparser.ConfigParser()
    
    # Read the configuration file
    config_read = config.read(path_config_file)
    
    # Raise error if configuration file does not exist
    if len(config_read) == 0:
        raise IOError('\n\nThere is no valid configuration file ' + path_config_file + '.')
          
    # Get parameters from [Basic pre-processing options] section
    sigma_filtering         = config_true_false(config.get('Basic pre-processing options', 'sigma_filtering'))    
    object_collapse_ndit    = config_true_false(config.get('Basic pre-processing options', 'object_collapse_ndit'))
    object_centering_method = config.get('Basic pre-processing options', 'object_centering_method')
    skip_preprocessing      = config_true_false(config.get('Basic pre-processing options', 'skip_preprocessing'))
    frames_to_remove        = literal_eval(config.get('Basic pre-processing options', 'frames_to_remove'))
      
    # Get parameters from [Basic post-processing options] section
    annulus_star                    = config_list_tuple(config.get('Basic post-processing options', 'annulus_star'))
    annulus_background              = config_list_tuple(config.get('Basic post-processing options', 'annulus_background'))
    normalized_polarization_images  = config_true_false(config.get('Basic post-processing options', 'normalized_polarization_images'))

    # Get parameters from [Advanced pre-processing options] section
    center_subtract_object    = config_true_false(config.get('Advanced pre-processing options', 'center_subtract_object'))
    center_param_centering    = literal_eval(config.get('Advanced pre-processing options', 'center_param_centering'))
    object_center_coordinates = config_list_tuple(config.get('Advanced pre-processing options', 'object_center_coordinates'))
    object_param_centering    = literal_eval(config.get('Advanced pre-processing options', 'object_param_centering'))
    flux_centering_method     = config.get('Advanced pre-processing options', 'flux_centering_method')
    flux_center_coordinates   = literal_eval(config.get('Advanced pre-processing options', 'flux_center_coordinates'))
    flux_param_centering      = literal_eval(config.get('Advanced pre-processing options', 'flux_param_centering'))
    flux_annulus_background   = config_list_tuple(config.get('Advanced pre-processing options', 'flux_annulus_background'))
    flux_annulus_star         = config_list_tuple(config.get('Advanced pre-processing options', 'flux_annulus_star'))
    
    # Get parameters from [Advanced post-processing options] section
    double_difference_type          = config.get('Advanced post-processing options', 'double_difference_type')
    single_posang_north_up          = config_true_false(config.get('Advanced post-processing options', 'single_posang_north_up'))
    convert_in_contrast_per_arcsec2 = config_true_false(config.get('Advanced post-processing options', 'convert_in_contrast_per_arcsec2'))

    return sigma_filtering, \
           object_collapse_ndit, \
           object_centering_method, \
           skip_preprocessing, \
           frames_to_remove, \
           annulus_star, \
           annulus_background, \
           normalized_polarization_images, \
           center_subtract_object, \
           center_param_centering, \
           object_center_coordinates, \
           object_param_centering, \
           flux_centering_method, \
           flux_center_coordinates, \
           flux_param_centering, \
           flux_annulus_background, \
           flux_annulus_star, \
           double_difference_type, \
           single_posang_north_up, \
           convert_in_contrast_per_arcsec2

###############################################################################
# wrapstr
###############################################################################

def wrapstr(string):
    '''
    Wrap a string to a maximum of 80 characters
    
    Input:
        string: string to be wrapped
        
    Return:
        string: wrapped string
    
    File written by Rob van Holstein
    Function status: verified  
    '''
    
    return textwrap.fill(string, width=80, 
                         replace_whitespace=False,
                         drop_whitespace=True,
                         break_long_words=False)
    
###############################################################################
# print_wrap
###############################################################################

def print_wrap(string):
    '''
    Print a string that is wrapped to a maximum of 80 characters
    
    Input:
        string: string to print
    
    File written by Rob van Holstein
    Function status: verified  
    '''
    
    print(wrapstr(string))

###############################################################################
# input_wrap
###############################################################################

def input_wrap(string):
    '''
    Redefine input() function so that it prints the string wrapped to a maximum 
    of 80 characters
    
    Input:
        string: string to use for input() function
    
    Output:
        result from input() function
        
    File written by Rob van Holstein
    Function status: verified  
    '''
    
    return input(textwrap.fill(string, width=80, 
                               replace_whitespace=False,
                               drop_whitespace=False,
                               break_long_words=False))

###############################################################################
# printandlog
###############################################################################

def printandlog(single_object, wrap=True):
    '''
    Print a single object (string) on screen and save it to a log file. The log
    file is located at path_log_file, which is a global variable to the 
    funtion.
    
    Input:
        single_object: single object (string) to be printed and logged
        wrap: if True, wrap print statements to a maximum of 80 character. If 
            False, do not wrap print statements (default = True).
           
    File written by Rob van Holstein
    Function status: verified
    '''
            
    if not os.path.exists(path_log_file):
        # Create log file
        open(path_log_file, 'w+')
    
    # Wrap string to not exceed 79 characters
    if type(single_object) == str and wrap == True:
        single_object = wrapstr(single_object)
    
    # Print object in log file and on screen
    print(single_object, file=open(path_log_file, 'a'))
    print(single_object) 

###############################################################################
# create_overview_headers
###############################################################################

def create_overview_headers(path_raw_dir, path_overview, log=True):
    '''
    Create an overview of relevant FITS-headers and write it to a text-file
    
    Input:
        path_raw_dir: string specifying path of raw directory
        path_overview: string specifying path of header overview to be created
        log: if True print and log statement that overview has been created, if
            False do not print anything
            
    File written by Rob van Holstein
    Function status: verified
    '''

    # Extract paths to FITS-files in raw directory
    path_raw_files = glob.glob(os.path.join(path_raw_dir,'*.fits'))
    
    # Extract headers
    header = [pyfits.getheader(x) for x in path_raw_files]
     
    # Sort raw files and headers based on observation date in headers
    date_obs = [x['DATE-OBS'] for x in header]
    sort_index = list(np.argsort(date_obs))
    path_raw_files = [path_raw_files[i] for i in sort_index]
    header = [header[i] for i in sort_index]
    
    # Define headers to be included in overview
    header_names = ['ESO OBS TARG NAME', 
                    'ESO DPR CATG', 
                    'ESO DPR TYPE', 
                    'ESO OCS DPI H2RT STOKES', 
                    'EXPTIME', 
                    'ESO DET NDIT', 
                    'ESO INS4 FILT2 NAME', 
                    'ESO INS1 FILT ID', 
                    'ESO INS COMB ICOR', 
                    'ESO INS4 DROT2 MODE', 
                    'ESO INS4 DROT3 MODE', 
                    'ESO TEL PARANG START', 
                    'ESO TEL ALT', 
                    'ESO INS4 DROT2 BEGIN', 
                    'ESO INS4 DROT2 POSANG', 
                    'ESO INS4 DROT3 BEGIN', 
                    'ESO INS4 DROT3 POSANG', 
                    'ESO INS4 DROT3 GAMMA', 
                    'ESO TEL AMBI FWHM START', 
                    'ESO TEL AMBI TAU0', 
                    'ESO INS1 OPTI2 NAME', 
                    'ESO INS4 OPTI8 NAME', 
                    'DATE-OBS']

    # Define the minimum separation between the column strings in the overview
    min_separation = 4
    
    # Create empty array to store overview in    
    m = len(header)
    n = len(header_names)
    print_array = np.empty((m+1, n+2), dtype=object)
    
    # Define header names and put in first line of overview
    header_names_print = [x.replace('ESO ', '').replace(' ', '.') for x in header_names]
    print_array[0, :] = ['FILE', 'NAME'] + header_names_print
    
    # Include file number and file name in first two columns of overview
    print_array[1:, 0] = [str(x) for x in range(1, m + 1)]
    print_array[1:, 1] = [os.path.basename(x) for x in path_raw_files]
    
    # Iterate over headers and header names to fill overview
    for i, header_sel in enumerate(header):   
        for j, header_name_sel in enumerate(header_names):
            if header_name_sel in header_sel:
                # If header exists, include it in overview
                print_array[i+1, j+2] = header_sel[header_name_sel]
            else:
                # If header does note exist, keep element of overview empty
                print_array[i+1, j+2] = ''
            
            if type(print_array[i+1, j+2]) is float:
                # Include string with 4 decimal places
                print_array[i+1, j+2] = '%.4f' % print_array[i+1, j+2]
            elif type(print_array[i+1, j+2]) is int:
                # Include string as an integer
                print_array[i+1, j+2] = '%s' % print_array[i+1, j+2]
    
    # Vectorize the function len()
    lenvec = np.vectorize(len)
    
    # Determine the length of each element of the overview and the maximum length in each column
    length_print = lenvec(print_array)
    max_length_print = np.max(length_print, axis=0)
    
    # Define a function to make a string of empty space to have each element at the right separation
    def make_separator(x): return x*' '
    
    # Vectorize the function and add empty spaces to each element in the overview
    make_separator_vec = np.vectorize(make_separator)
    separator_array = make_separator_vec(max_length_print - length_print + min_separation)
    print_array += separator_array
    
    # Save the overview to a text file
    np.savetxt(path_overview, print_array, fmt = '%s', newline= '\r\n')

    # Save the overview to a csv file
#    df_array = pd.DataFrame(print_array[1:,1:],columns=print_array[0,1:],index=print_array[1:, 0])
    df_array = pd.DataFrame(print_array[1:,1:],columns=print_array[0,1:],index=print_array[1:, 0])
    df_array.index.name = print_array[0, 0]
    df_array.to_csv(path_overview.replace('.txt','.csv'))

    if log:
        printandlog('\nWrote files ' + path_overview + ' and ' + \
                    path_overview.replace('.txt','.csv') + ' showing overviews of the relevant header keywords for each file in the raw directory.')

###############################################################################
# check_own_programs
###############################################################################

def check_own_programs(header):
    '''
    Check if data is from one of my own programs to print message
    
    Input:
        header: header or list of FITS-headers of raw files
    
    File written by Rob van Holstein
    Function status: verified
    '''
    
    own_programs = ['098.C-0790(A)',
                    '0101.C-0502(A)', 
                    '0101.C-0502(B)',
                    '0102.C-0871(A)',
                    '0102.C-0916(A)', 
                    '0102.C-0916(B)', 
                    '0102.C-0916(C)', 
                    '2102.C-5016(B)']
    
    program_id = [x['ESO OBS PROG ID'] for x in header]
    
    if any([x in own_programs for x in program_id]):
        printandlog('\nTerminating reduction.')
        for i in range(20):
            printandlog('\n')
        printandlog('\n###############################################################################')
        printandlog('Please note')                                
        printandlog('###############################################################################')
        printandlog('\nDear colleague,')      
        printandlog('\nThe data you are trying to reduce comes from one of my observing programs. ' + 
                    'In case you are interested in this data, please contact me at \'vanholstein@strw.leidenuniv.nl\'. ' +
                    'I am always more than happy to collaborate on projects!')
        printandlog('\nPlease note that the calibration of the IRDIS ' +
                    'polarimetric mode (which is the core of IRDAP) and the development of IRDAP itself have taken 4 years ' +
                    '(started early 2015). I would therefore very much appreciate if you could take this into ' +
                    'consideration. Thank you for your understanding.')
        printandlog('\nKind regards,')
        printandlog('\nRob van Holstein')
        printandlog('\n###############################################################################')
        for i in range(20):
            printandlog('\n')           
                    
        sys.exit()
        
###############################################################################
# check_sort_data_create_directories
###############################################################################

def check_sort_data_create_directories(frames_to_remove=[], 
                                       combination_method_polarization='least squares', 
                                       object_centering_method='automatic', 
                                       save_preprocessed_data=True, 
                                       show_images_center_coordinates=True):
    '''
    Check the FITS-headers of the data in the raw directory, remove files and 
    frames as specified by the user, sort the data and create directories to 
    write processed data to.
    
    Input:
        frames_to_remove: list of integers and length-2-tuples of integers 
            indicating which files and frames to remove (0-based). A complete
            file can be removed by specifying its integer index, while a frame
            of specific file can be removed by specifying a tuple
            (file_index, frame_index) (default = []).
        combination_method_polarization: method to be used to produce the 
            incident Q- and U-images, 'least squares', 'trimmed mean' or 'median'. 
            In this function, if combination_method_polarization is forced 
            to 'least squares' in case there is an unequal number of Q- and 
            U-measurements (default = 'least squares').
        object_centering_method: method to center the OBJECT-frames. In this 
            function, if object_centering_method is 'automatic', it is set to
            'center frames' if there are CENTER-files, and is set to '
            gaussian' if there are no CENTER-files. If object_centering_method
            is 'center frames' or 'manual', use fixed coordinates as provided by 
            center_coordinates. If 'gaussian', fit a 2D Gaussian to each frame. 
            If 'cross-correlation', fit a 2D Gaussian to the first frame and then
            use cross-correlation to align (register) the other frames onto the 
            centered first  frame. For 'gaussian' and 'cross-correlation' 
            center_coordinates is used as initial guess of the center 
            coordinates and the determined center coordinates are plotted for
            each image (default = 'automatic').
        save_preprocessed_data: If True, save preprocessed cubes of single-sum 
            and single-difference images in the 'preprocessed' folder so that 
            the preprocessing can be skipped when re-running the pipeline
            (default = True).
        show_images_center_coordinates: If True, plot the sub-images showing the 
            center coordinates for each frame. The plots allow for checking
            whether the centering is correct and to scan the data for frames 
            with bad quality (default = True). 
    
    Note that combination_method_polarization, object_centering_method, 
    save_preprocessed_data and show_images_center_coordinates are input to this 
    function as they are required for the sorting of the data or preparing the
    pre-processing, not to actually perform centering for example.
    
    Note that path_raw_dir, path_flat_dir, path_bpm_dir, path_sky_dir, 
    path_center_dir, path_flux_dir, path_sky_flux_dir, path_preprocessed_dir, 
    path_reduced_dir, path_reduced_star_pol_subtr_dir and path_overview are 
    global variables to the function.
    
    Output:
        path_dark_files: list of paths to raw DARK(,BACKGROUND)-files
        path_flat_files: list of paths to raw FLAT-files
        path_object_files: list of paths to raw OBJECT-files
        path_sky_files: list of paths to raw SKY-files for OBJECT
        path_center_files: list of paths to raw CENTER-files
        path_object_center_files: list of paths to raw OBJECT-files to be 
            subtracted from raw CENTER-files
        path_flux_files: list of paths to raw FLUX-files
        path_sky_flux_files: list of paths to raw SKY-files for FLUX
        indices_to_remove_object: list of 1-D arrays with indices of frames to 
            be removed for each OBJECT-file. If no frames are to be removed the 
            array is empty. 
        indices_to_remove_sky: list of 1-D arrays with indices of frames to 
            be removed for each SKY-file. If no frames are to be removed the 
            array is empty. 
        indices_to_remove_center: list of 1-D arrays with indices of frames to 
            be removed for each CENTER-file. If no frames are to be removed the 
            array is empty. 
        indices_to_remove_object_center: list of 1-D arrays with indices of frames 
            to be removed for each OBJECT-file to be subtracted from a CENTER-file. 
            If no frames are to be removed the array is empty. 
        indices_to_remove_flux: list of 1-D arrays with indices of frames to 
            be removed for each FLUX-file. If no frames are to be removed the 
            array is empty. 
        indices_to_remove_sky_flux: list of 1-D arrays with indices of frames to 
            be removed for each SKY-file to be subtracted from the FLUX-file(s). 
            If no frames are to be removed the array is empty.        
        file_index_object: list of file indices of OBJECT-files (0-based)
        file_index_flux: list of file indices of FLUX-files (0-based)
        object_centering_method: method to center the OBJECT-frames, 'center frames', 
            'gaussian', cross-correlation' or 'manual'.
        combination_method_polarization: method to be used to produce the 
            incident Q- and U-images, 'least squares', 'trimmed mean' or 'median'

    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: verified
    '''

    ###############################################################################
    # Obtain raw files and headers
    ###############################################################################
   
    # Extract paths to FITS-files in raw directory
    printandlog('\nReading raw directory ' + path_raw_dir + '.')
    path_raw_files = glob.glob(os.path.join(path_raw_dir,'*.fits'))

    # Create overview of relevant FITS-headers in a text file
    create_overview_headers(path_raw_dir, path_overview, log=True)
    
    # Extract headers
    header = [pyfits.getheader(x) for x in path_raw_files]
    check_own_programs(header)
    
    # Sort raw files and headers based on observation date in headers
    date_obs = [x['DATE-OBS'] for x in header]
    sort_index = list(np.argsort(date_obs))
    path_raw_files = [path_raw_files[i] for i in sort_index]
    header = [header[i] for i in sort_index]

    ###############################################################################
    # Remove files as specified by the user
    ###############################################################################

    # Define file indices
    file_index = [x for x in range(0, len(path_raw_files))]

    # Extract list of complete files to remove            
    files_to_remove  = [x for x in frames_to_remove if type(x) is int]

    # Raise error if files to remove are listed more than once
    if len(files_to_remove) != len(set(files_to_remove)):
        raise ValueError('\n\nOne or more files to be removed are listed more than once.')
        
    # Raise error if files to remove do not exist
    if not all([x in file_index for x in files_to_remove]):
        raise ValueError('\n\nOne or more files to be removed do not exist.')
    
    # Print which files will be removed as specified by the user; do not change len() != 0 to any() because any([0]) = False and then it doesn't work
    if len(files_to_remove) != 0:
        printandlog('\nAs requested by the user, the following FITS-file(s) will not be used:')
        printandlog('File     Name')
        for file_to_remove in files_to_remove:
            separation = 9 - len(str(file_index[file_to_remove]))
            printandlog(str(file_index[file_to_remove] + 1) + separation*' ' + os.path.basename(path_raw_files[file_to_remove]), wrap=False)

    # Remove files as specified by the user
    path_raw_files = [x for i,x in enumerate(path_raw_files) if i not in files_to_remove]  
    header = [x for i,x in enumerate(header) if i not in files_to_remove]  
    file_index = [x for i,x in enumerate(file_index) if i not in files_to_remove]  
    
    ###############################################################################
    # Check the headers of all files
    ###############################################################################
    
    # Perform checks on header values that apply to all data
    printandlog('\nChecking the headers of all files.')

    if not all([x['INSTRUME'] == 'SPHERE' for x in header]):
        raise IOError('\n\nOne or more files are not taken with SPHERE.')

    if not all([x['ESO DET NAME'] == 'IRDIS' for x in header]):
        raise IOError('\n\nOne or more files are not taken with IRDIS.')

    if not all([x['ESO DPR CATG'] in ['SCIENCE', 'CALIB'] for x in header]):
        raise IOError('\n\nOne or more files are not of category SCIENCE or CALIB.')

    if not all([x['ESO DPR TYPE'] in ['OBJECT', 'SKY', 'OBJECT,CENTER', 'OBJECT,FLUX', \
                                      'DARK', 'DARK,BACKGROUND', 'FLAT,LAMP'] for x in header]):
        raise IOError('\n\nOne or more files are not of type OBJECT, SKY, OBJECT,CENTER, OBJECT,FLUX, DARK(,BACKGROUND) or FLAT,LAMP.')

    # Perform checks on header values that apply to all data except for DARK(,BACKGROUND)- and FLAT-files 
    header_on_sky = [x for x in header if x['ESO DPR TYPE'] in ['OBJECT', 'SKY', 'OBJECT,CENTER', 'OBJECT,FLUX']]
    
    if not any(header_on_sky):
        raise IOError('\n\nThere is no on-sky data.')
    
    if len(set([x['ESO OBS TARG NAME'] for x in header_on_sky])) != 1:
        different_targets = input_wrap('\nThe on-sky data provided have different targets. Continue anyway? (y/n) ')
        if different_targets == 'y':
            printandlog('\nWARNING, continuing reduction although the on-sky data provided have different targets.')
        elif different_targets == 'n':
            raise IOError('\n\nThe on-sky data provided have different targets.')
        else:
            raise IOError('\n\nThe provided input \'' + str(different_targets) + '\' is not valid.')
        
    if len(set([x['ESO INS4 COMB ROT'] for x in header_on_sky])) != 1:
        raise IOError('\n\nThe on-sky data provided use different tracking modes.')

    if not header_on_sky[0]['ESO INS4 COMB ROT'] in ['FIELD', 'PUPIL']:
        raise IOError('\n\nThe tracking mode of the on-sky data is not field-tracking or pupil-tracking.')
        
    if not all([x['ESO INS4 OPTI8 NAME'] == 'H_NIR' for x in header_on_sky]):
        raise IOError('\n\nOne or more files of the on-sky data do not have the NIR half-wave plate inserted.')

    if not all([x['ESO INS1 OPTI2 NAME'] == 'P0-90' for x in header_on_sky]):
        raise IOError('\n\nOne or more files of the on-sky data do not have the P0-90 polarizer set inserted.')
        
    if len(set([x['ESO INS1 FILT ID'] for x in header_on_sky])) != 1:
        raise IOError('\n\nThe on-sky data provided have different filters.')

    if header_on_sky[0]['ESO INS1 FILT ID'] in ['FILT_NBF_HeI', 'FILT_NBF_ContJ', 'FILT_NBF_PaB', 'FILT_NBF_ContH', 'FILT_NBF_FeII', 'FILT_NBF_ContK1', 'FILT_NBF_H2', 'FILT_NBF_BrG', 'FILT_NBF_CntK2', 'FILT_NBF_CO']:
        narrowband_continue = input_wrap('\nThe on-sky data uses a narrowband filter. The model-based correction method for the narrowband filters will be added in the near future. For now, the narrowband data can be reduced using the model of the broadband filter that encompasses the wavelength range of the narrowband filter. However, the correction of the instrumental polarization effects will not be as accurate as for the broadband filters. If continuing, note that for accurate flat fielding it is highly recommended to construct a master flat for the narrowband filter by including a set of FLAT- and DARK(,BACKGROUND)-files in the raw subdirectory. If not, a static flat of the corresponding broadband filter will be used. Finally, you may need to adapt the input parameters \'center_param_centering\', \'object_center_coordinates\' and \'object_param_centering\' to make IRDAP correctly center the data, and change the input parameter \'annulus_star\' to better position the annulus to compute the stellar polarization from. Reduce the narrowband data using the model of the corresponding broadband filter? (y/n) ')
        if narrowband_continue == 'y':
            printandlog('\nWARNING, continuing reduction of narrowband data although the results will not be as accurate as for the broadband filters.')
        elif narrowband_continue == 'n':
            raise IOError('\n\nAborting data reduction.')
        else:
            raise IOError('\n\nThe provided input \'' + str(narrowband_continue) + '\' is not valid.')
    elif not header_on_sky[0]['ESO INS1 FILT ID'] in ['FILT_BBF_Y', 'FILT_BBF_J', 'FILT_BBF_H', 'FILT_BBF_Ks']:
        raise IOError('\n\nThe filter used for the on-sky data is not broadband Y, J, H or Ks or any of the narrowband filters.')
     
    # Perform checks on header values that apply only to OBJECT files
    header_object = [x for x in header if x['ESO DPR TYPE'] == 'OBJECT']
        
    if not any(header_object):
        raise IOError('\n\nThere are no OBJECT-files.')
        
    if len(set([x['EXPTIME'] for x in header_object])) != 1:
        raise IOError('\n\nThe OBJECT-files have different exposure times.')
    
    if len(set([x['ESO INS4 FILT2 NAME'] for x in header_object])) != 1:
        raise IOError('\n\nThe OBJECT-files use different NIR neutral density filters.')
    
    if header_object[0]['ESO INS4 FILT2 NAME'] != 'OPEN':
        printandlog('\nWARNING, the OBJECT-files use a NIR neutral density filter. The final data product will be less accurate because the neutral density filters are known to have a depolarizing effect that is not calibrated.')
    
    if len(set([x['ESO INS COMB ICOR'] for x in header_object])) != 1:
        raise IOError('\n\nThe OBJECT-files use different coronagraph settings.')
    
    # Determine exposure time and NIR neutral density filter for OBJECT files  
    object_filter = header_object[0]['ESO INS1 FILT ID']
    object_exposure_time = header_object[0]['EXPTIME']
    object_nd_filter = header_object[0]['ESO INS4 FILT2 NAME']
    
    # Perform checks on header values that apply only to FLUX files
    header_flux = [x for x in header if x['ESO DPR TYPE'] == 'OBJECT,FLUX']
    
    if any(header_flux):
        if len(set([x['EXPTIME'] for x in header_flux])) != 1:
            raise IOError('\n\nThe FLUX-files have different exposure times.')
        
        if len(set([x['ESO INS4 FILT2 NAME'] for x in header_flux])) != 1:
            raise IOError('\n\nThe FLUX-files use different NIR neutral density filters.')
        
        # Determine exposure time and NIR neutral density filter for FLUX files     
        flux_exposure_time = header_flux[0]['EXPTIME']
        flux_nd_filter = header_flux[0]['ESO INS4 FILT2 NAME']
    
    else:
        # Set exposure time and NIR neutral density filter for FLUX files to None
        printandlog('\nWARNING, there are no FLUX-files.')
        flux_exposure_time = None
        flux_nd_filter = None
        
    # Perform checks on header values that apply only to SKY files
    header_sky = [x for x in header if x['ESO DPR TYPE'] == 'SKY']
    
    if any(header_sky):
        if not all([x['EXPTIME'] in [object_exposure_time, flux_exposure_time] for x in header_sky]):
            if any(header_flux):
                raise IOError('\n\nOne or more SKY-files have an exposure time different from that of the OBJECT- or FLUX-files.')
            else:
                raise IOError('\n\nOne or more SKY-files have an exposure time different from that of the OBJECT-files.')
    
        if not all([x['ESO INS4 FILT2 NAME'] in [object_nd_filter, flux_nd_filter] for x in header_sky]):
            if any(header_flux):
                raise IOError('\n\nOne or more SKY-files use a NIR neutral density filter different from that of the OBJECT- or FLUX-files.')
            else:
                raise IOError('\n\nOne or more SKY-files use a NIR neutral density filter different from that of the OBJECT-files.')
    
    if not any([x['ESO DPR TYPE'] == 'SKY' and x['EXPTIME'] == object_exposure_time and x['ESO INS4 FILT2 NAME'] == object_nd_filter for x in header]):
        printandlog('\nWARNING, there are no SKY-files to subtract from the OBJECT-files. While this is generally no problem for the final polarization images, the final total intensity images can show some artifacts.')    
    
    if any(header_flux):    
        if not any([x['ESO DPR TYPE'] == 'SKY' and x['EXPTIME'] == flux_exposure_time and x['ESO INS4 FILT2 NAME'] == flux_nd_filter for x in header]):
            printandlog('\nWARNING, there are no SKY-files to subtract from the FLUX-file(s). Although the background will be subtracted after determining it using the annulus as defined by the input variable \'flux_annulus_background\', the result will be less accurate than when subtracting a SKY-image.')
    
    # Perform checks on header values that apply only to CENTER files
    header_center = [x for x in header if x['ESO DPR TYPE'] == 'OBJECT,CENTER']
    
    if any(header_center):
        if not all([x['EXPTIME'] == object_exposure_time for x in header_center]):
            raise IOError('\n\nOne or more CENTER-files have an exposure time different from that of the OBJECT-files.')
        
        if not all([x['ESO INS4 FILT2 NAME'] == object_nd_filter for x in header_center]):
            raise IOError('\n\nOne or more CENTER-files use a NIR neutral density filter different from that of the OBJECT-files.')
         
    # Perform checks on header values that apply only to DARK(,BACKGROUND)- and FLAT-files
    header_dark_all = [x for x in header if x['ESO DPR TYPE'] in ['DARK', 'DARK,BACKGROUND']]
    header_dark_background = [x for x in header_dark_all if x['ESO DPR TYPE'] == 'DARK,BACKGROUND']
    header_flat = [x for x in header if x['ESO DPR TYPE'] == 'FLAT,LAMP']
               
    if any(header_dark_all) and any(header_flat):    
        
        if len(header_dark_all) != len(header_flat):
            raise IOError('\n\nThe number of DARK(,BACKGROUND)- and FLAT-files is unequal.')

        if len(header_flat) == 1:
            raise IOError('\n\nThere is only one pair of DARK(,BACKGROUND)- and FLAT-files. To make the bad pixel map we need at least two pairs with each a different exposure time.')

        if any([x['ESO INS1 FILT ID'] != object_filter for x in header_dark_background]):
            raise IOError('\n\nOne or more DARK,BACKGROUND-files have a different filter than the OBJECT-files.')
        
        if any([x['ESO INS1 FILT ID'] != object_filter for x in header_flat]):
            raise IOError('\n\nOne or more FLAT-files have a different filter than the OBJECT-files.')

        if len(set([x['ESO INS1 OPTI2 NAME'] for x in header_flat])) != 1:
            raise IOError('\n\nThe FLAT-files use different setups for the filter wheel containing the polarizer set.')
            
        if not header_flat[0]['ESO INS1 OPTI2 NAME'] in ['P0-90', 'CLEAR']:
            raise IOError('\n\nThe FLAT-files do not use the P0-90 polarizer set or no filter.')

        if header_flat[0]['ESO INS1 OPTI2 NAME'] != 'P0-90':
            printandlog('\nWARNING, the FLAT-files do not have the P0-90 polarizer set inserted. The polarizer set causes strong vignetting at the edges of the field of view, which will not be corrected with these flats.')

        if set([x['EXPTIME'] for x in header_dark_all]) != set([x['EXPTIME'] for x in header_flat]):
            raise IOError('\n\nThe exposure time of one or more FLAT-files cannot be matched to the exposure time of a DARK(,BACKGROUND)-file or vice versa.')

        if any(np.sort([x['EXPTIME'] for x in header_dark_all]) != np.sort([x['EXPTIME'] for x in header_flat])):
            raise IOError('\n\nFor one or more exposure times, the number of DARK(,BACKGROUND)- and FLAT-files are different.')
            
        mjd_dark_background = [x['MJD-OBS'] for x in header_dark_background]
        mjd_flat = [x['MJD-OBS'] for x in header_flat] 
        
        if header_flat[0]['ESO INS1 FILT ID'] in ['FILT_BBF_Ks', 'FILT_NBF_ContK1', 'FILT_NBF_H2', 'FILT_NBF_BrG', 'FILT_NBF_CntK2', 'FILT_NBF_CO']:
            if max(mjd_flat) - min(mjd_dark_background) > 1/24 or max(mjd_dark_background) - min(mjd_flat) > 1/24:
                printandlog('\nWARNING, (some of) the DARK,BACKGROUND- and FLAT-files are taken more than 1 hour apart and use the Ks-band filter or one of the narrowband filters that are close in wavelength to Ks. The background in these files in Ks-band is known to generally be strongly time-varying.')

    if any(header_dark_all) and not any(header_flat):
        raise IOError('\n\nThere are DARK(,BACKGROUND)-files but no FLAT-files.')
        
    if any(header_flat) and not any(header_dark_all):
        raise IOError('\n\nThere are FLAT-files but no DARK(,BACKGROUND)-files.')
        
    # Print that headers have been checked and passed all checks
    printandlog('\nThe FITS-headers of the raw data have passed all checks.')
        
    ###############################################################################
    # Define list of indices to remove frames
    ###############################################################################

    # Create list of files to remove frames from and frames to remove for each of these files    
    files_to_remove_frames_from = np.array([x[0] for x in frames_to_remove if type(x) is tuple])
    frames_to_remove_per_file = np.array([x[1] for x in frames_to_remove if type(x) is tuple])

    # Raise error if frames to be removed are listed more than once
    frames_to_remove_tuples = [x for x in frames_to_remove if type(x) is tuple]
    if len(frames_to_remove_tuples) != len(set(frames_to_remove_tuples)):
        raise ValueError('\n\nOne or more frames to be removed are listed more than once.')

    # Raise error if frames need to be removed from a file that is already completely removed
    if any([x in files_to_remove for x in files_to_remove_frames_from]):
        raise ValueError('\n\nOne or more files to remove frames from are already completely removed.')  

    # Raise error if files to remove do not exist
    if not all([x in file_index for x in files_to_remove_frames_from]):
        raise ValueError('\n\nOne or more files to remove frames from do not exist.')

    # Raise error if for one or more files the indices of frames to remove are negative or zero
    if any([x < 0 for x in frames_to_remove_per_file]):
        raise ValueError('\n\nOne or more of the frame numbers to be removed are negative or zero.')

    # Create list of arrays with indices of frames to remove for each file
    indices_to_remove = [np.sort(frames_to_remove_per_file[files_to_remove_frames_from == x]).astype(np.int) for x in file_index]

    # Extract NDIT of each file
    NDIT = [x['ESO DET NDIT'] for x in header]
    
    # Raise error if one or more files have all the frames removed
    if any([np.array_equal(x, np.arange(0, y)) for x,y in zip(indices_to_remove, NDIT)]):
        raise ValueError('\n\nOne or more files have all frames removed. If you want to remove the complete file, please specify just the file number in frames_to_remove and not a tuple of frames.')
    
    # Raise error if for one or more files the indices of frames to remove are outside of the NDIT
    if any([any(x >= y) for x,y in zip(indices_to_remove, NDIT)]):
        raise ValueError('\n\nOne or more of the frame numbers to be removed are higher than the NDIT of the corresponding file.')

    # Print which frames will be removed as specified by the user; do not change len() != 0 to any() because any([0]) = False and then it doesn't work
    if len(files_to_remove_frames_from) != 0:
        printandlog('\nAs requested by the user, the following frame(s) will not be used:')
        separation_title = max(list(map(len, [os.path.basename(x) for x,y in zip(path_raw_files, indices_to_remove) if y.size > 0]))) + 1
        printandlog('File     Name' + separation_title*' ' + 'Frame')
        for file_index_sel, path_sel, indices_sel in zip(file_index, path_raw_files, indices_to_remove):
            if indices_sel.size > 0:
                separation_1 = 9 - len(str(file_index_sel + 1))
                separation_2 = separation_title + 4 - len(os.path.basename(path_sel))
                printandlog(str(file_index_sel + 1) + separation_1*' ' + os.path.basename(path_sel) + separation_2*' ' + ', '.join(map(str, indices_sel + 1)), wrap=False)

    ###############################################################################
    # Sort file paths and indices of frames to be removed according to file type
    ###############################################################################

    # Create empty lists
    path_dark_files = []
    path_flat_files = []
    path_object_files = []
    path_sky_files = []
    path_center_files = []
    path_flux_files = []
    path_sky_flux_files = [] 
    path_imcompatible_files = []
    indices_to_remove_dark = []
    indices_to_remove_flat = []    
    indices_to_remove_object = []
    indices_to_remove_sky = []    
    indices_to_remove_center = []
    indices_to_remove_flux = []
    indices_to_remove_sky_flux = []    
    file_index_object = []
    file_index_flux = []    
    stokes_parameter = []
    mjd_half_object = []    
    mjd_half_center = []    

    # Sort file paths and indices of frames to be removed according to file type 
    for file_sel, header_sel, file_index_sel, NDIT_sel, indices_sel in zip(path_raw_files, header, file_index, NDIT, indices_to_remove):
    
        if header_sel['ESO DPR TYPE'] in ['DARK', 'DARK,BACKGROUND']:
            path_dark_files.append(file_sel)
            indices_to_remove_dark.append(indices_sel)

        elif header_sel['ESO DPR TYPE'] == 'FLAT,LAMP':
            path_flat_files.append(file_sel)
            indices_to_remove_flat.append(indices_sel)
            
        elif header_sel['ESO DPR TYPE'] == 'OBJECT':
            path_object_files.append(file_sel)
            indices_to_remove_object.append(indices_sel)
            file_index_object.append(file_index_sel)           

            # Append Stokes parameter to list
            stokes_parameter.append(header_sel['ESO OCS DPI H2RT STOKES'])   

            # Calculate mean Julian date halfway the exposure 
            mjd = header_sel['MJD-OBS']
            file_execution_time = NDIT_sel * (0.938 + object_exposure_time) + 2.4
            mjd_half_object.append(mjd + 0.5 * file_execution_time / msd)    
            
        elif header_sel['ESO DPR TYPE'] == 'SKY' and \
           header_sel['EXPTIME'] == object_exposure_time and \
           header_sel['ESO INS4 FILT2 NAME'] == object_nd_filter:
            path_sky_files.append(file_sel)          
            indices_to_remove_sky.append(indices_sel)
            
        elif header_sel['ESO DPR TYPE'] == 'OBJECT,CENTER' and \
           header_sel['EXPTIME'] == object_exposure_time and \
           header_sel['ESO INS4 FILT2 NAME'] == object_nd_filter:
            path_center_files.append(file_sel)          
            indices_to_remove_center.append(indices_sel)
            
            # Calculate mean Julian date halfway the exposure 
            mjd = header_sel['MJD-OBS']
            file_execution_time = NDIT_sel * (0.938 + object_exposure_time) + 1.4
            mjd_half_center.append(mjd + 0.5 * file_execution_time / msd)

        elif header_sel['ESO DPR TYPE'] == 'OBJECT,FLUX':
            path_flux_files.append(file_sel)
            indices_to_remove_flux.append(indices_sel)
            file_index_flux.append(file_index_sel)
            
        elif header_sel['ESO DPR TYPE'] == 'SKY' and \
           header_sel['EXPTIME'] == flux_exposure_time and \
           header_sel['ESO INS4 FILT2 NAME'] == flux_nd_filter:
            path_sky_flux_files.append(file_sel)          
            indices_to_remove_sky_flux.append(indices_sel)
            
        else:
            path_imcompatible_files.append(file_sel)
   
    # Sort DARK(,BACKGROUND)- and FLAT-files based on exposure time in headers
    exposure_time_dark = [pyfits.getheader(x)['EXPTIME'] for x in path_dark_files]
    sort_index = list(np.argsort(exposure_time_dark))
    path_dark_files = [path_dark_files[i] for i in sort_index]
    
    exposure_time_flat = [pyfits.getheader(x)['EXPTIME'] for x in path_flat_files]
    sort_index = list(np.argsort(exposure_time_flat))    
    path_flat_files = [path_flat_files[i] for i in sort_index]
    
    # Find object-files that are closest in time to center-files
    path_object_center_files = []
    indices_to_remove_object_center = []
    
    for mjd_half_center_sel in mjd_half_center:
        index = np.argmin(np.abs(np.array(mjd_half_object) - np.array(mjd_half_center_sel)))
        path_object_center_files.append(path_object_files[index])
        indices_to_remove_object_center.append(indices_to_remove_object[index])

    ###############################################################################
    # Remove files that lack a Q/U^+/- counterpart
    ###############################################################################

    # Find object-files that do not form a pair with their Q/U^+/- counterpart
    files_to_remove_stokes = []
    n = len(stokes_parameter)
    
    for i in range(n):
        if stokes_parameter[i] == 'Qplus':
            if i + 1 == n:
                files_to_remove_stokes.append(i) 
            elif stokes_parameter[i + 1] != 'Qminus':        
                files_to_remove_stokes.append(i)
        if stokes_parameter[i] == 'Qminus':
            if stokes_parameter[i - 1] != 'Qplus':        
                files_to_remove_stokes.append(i)
        if stokes_parameter[i] == 'Uplus':
            if i + 1 == n:
                files_to_remove_stokes.append(i)
            elif stokes_parameter[i + 1] != 'Uminus':        
                files_to_remove_stokes.append(i)
        if stokes_parameter[i] == 'Uminus':
            if stokes_parameter[i - 1] != 'Uplus':        
                files_to_remove_stokes.append(i)

    # Print which object-files will be removed because they lack a Q/U^+/- counterpart; do not change len() != 0 to any() because any([0]) = False and then it doesn't work
    if len(files_to_remove_stokes) != 0:
        printandlog('\nWARNING, the following FITS-file(s) will not be used because of a missing Q^+/- or U^+/-counterpart:')
        printandlog('File     Name')
        for file_to_remove in files_to_remove_stokes:
            separation = 9 - len(str(file_index_object[file_to_remove] + 1))
            printandlog(str(file_index_object[file_to_remove] + 1) + separation*' ' + os.path.basename(path_object_files[file_to_remove]), wrap=False)
 
    # Print warning if one or more files are removed that would have frames removed
    if any([len(x) != 0 for i,x in enumerate(indices_to_remove_object) if i in files_to_remove_stokes]):
        printandlog('\nWARNING, the user has specified frames to be removed for one or more FITS-files that are removed because of a missing Q^+/- or U^+/-counterpart.')
                    
    # Remove files that lack a Q/U^+/- counterpart
    path_object_files = [x for i,x in enumerate(path_object_files) if i not in files_to_remove_stokes]  
    stokes_parameter = [x for i,x in enumerate(stokes_parameter) if i not in files_to_remove_stokes]  
    indices_to_remove_object = [x for i,x in enumerate(indices_to_remove_object) if i not in files_to_remove_stokes] 
    file_index_object = [x for i,x in enumerate(file_index_object) if i not in files_to_remove_stokes]    

    # Check if there are both Q- and U-measurements
    if not all([x in stokes_parameter for x in ['Qplus', 'Qminus', 'Uplus', 'Uminus']]):
        if any([x in stokes_parameter for x in ['Qplus', 'Qminus']]):
            raise IOError('\n\nThe data has no U-measurements and therefore cannot be reduced.')
        if any([x in stokes_parameter for x in ['Uplus', 'Uminus']]):
            raise IOError('\n\nThe data has no Q-measurements and therefore cannot be reduced.')
        
    # Force 'least squares' for combining the polarization images if the number of Q- and U-measurements is unequal
    if combination_method_polarization != 'least squares':
        if stokes_parameter.count('Qplus') != stokes_parameter.count('Uplus'):
            combination_method_polarization = 'least squares'
            printandlog('\nWARNING, combination_method_polarization is forced to \'least squares\' because the number of Q- and U-measurements are unequal. The plots showing the measured star polarization as a function of HWP cycle number will only show the complete HWP cycles.')
 
    ###############################################################################
    # Print number of files for each file type
    ###############################################################################

    # Print number of files for each file type
    printandlog('\nNumber of files for each file type:')
    printandlog('DARK(,BACKGROUND):    ' + str(len(path_dark_files)))    
    printandlog('FLAT,LAMP:            ' + str(len(path_flat_files)))
    printandlog('OBJECT (O):           ' + str(len(path_object_files)))
    printandlog('SKY (S) for OBJECT:   ' + str(len(path_sky_files)))
    printandlog('CENTER (C):           ' + str(len(path_center_files)))
    printandlog('FLUX (F):             ' + str(len(path_flux_files)))
    printandlog('SKY (S) for FLUX:     ' + str(len(path_sky_flux_files)))
    
    # Print warning when one or more raw files do not fall under any of the categories above
    if len(path_imcompatible_files) != 0:
        printandlog('\nWARNING, one or more files do not fall under any of the file type categories listed above and will be ignored:')
        for file_sel in path_imcompatible_files:
            printandlog('{0:s}'.format(file_sel))
    
    # Raise error if there are no OBJECT-files to reduce
    if not any(path_object_files):
        raise IOError('\n\nThere are no OBJECT-files to reduce.')
    
    # If object_centering_method is 'automatic', set to 'center frames' if there are CENTER-files and otherwise to 'gaussian'
    if object_centering_method == 'automatic':
        if any(path_center_files):
            printandlog('\nobject_centering_method is \'automatic\': changing it to \'center frames\', because there are CENTER-files.')
            object_centering_method = 'center frames'
        else:
            printandlog('\nobject_centering_method is \'automatic\': changing it to \'gaussian\', because there are no CENTER-files.')
            object_centering_method = 'gaussian'
              
    # Raise error when there are no center files, but they are required by the selected centering method
    if object_centering_method == 'center frames' and not any(path_center_files):
        raise IOError('\n\nobject_centering_method = \'{0:s}\' (or \'automatic\'), but there are no CENTER-files provided.'.format(object_centering_method))

    # Raise warning when there are no center files, but the data is coronagraphic
    coronagraph_used = header_object[0]['ESO INS COMB ICOR']
    
    if not any(path_center_files) and coronagraph_used != 'N_NS_CLEAR':
        printandlog('\nWARNING, the data uses a coronagraph, but there are no CENTER-files. Although challenging, in most cases the data can still be centered with some accuracy by manually defining the centers of the left and right images using the input parameter \'object_center_coordinates\', and by setting a sub-image crop radius of around 20 pixels with the input parameter \'object_param_centering\'.')
        
    ###############################################################################
    # Create directories to write processed data to
    ###############################################################################

    # Create directories to write processed data to
    directories_created = []
    directories_already_existing = []

    if any(path_flat_files):
        if not os.path.exists(path_flat_dir):
            os.makedirs(path_flat_dir)
            directories_created.append(path_flat_dir)
        else:
            directories_already_existing.append(path_flat_dir)
        if not os.path.exists(path_bpm_dir):
            os.makedirs(path_bpm_dir)
            directories_created.append(path_bpm_dir)
        else:
            directories_already_existing.append(path_bpm_dir)
       
    if any(path_sky_files):
        if not os.path.exists(path_sky_dir):
            os.makedirs(path_sky_dir)
            directories_created.append(path_sky_dir)
        else:
            directories_already_existing.append(path_sky_dir)
    
    if any(path_center_files) and object_centering_method == 'center frames':
        if not os.path.exists(path_center_dir):
            os.makedirs(path_center_dir)
            directories_created.append(path_center_dir) 
        else:
            directories_already_existing.append(path_center_dir)
        
    if any(path_flux_files):
        if not os.path.exists(path_flux_dir):
            os.makedirs(path_flux_dir)
            directories_created.append(path_flux_dir)
        else:
            directories_already_existing.append(path_flux_dir)
        
    if any(path_sky_flux_files):
        if not os.path.exists(path_sky_flux_dir):
            os.makedirs(path_sky_flux_dir)
            directories_created.append(path_sky_flux_dir)
        else:
            directories_already_existing.append(path_sky_flux_dir)
       
    if save_preprocessed_data == True or show_images_center_coordinates == True or \
       object_centering_method in ['gaussian', 'cross-correlation']:
        if not os.path.exists(path_preprocessed_dir):
            os.makedirs(path_preprocessed_dir)
            directories_created.append(path_preprocessed_dir)
        else:
            directories_already_existing.append(path_preprocessed_dir)
    
    if not os.path.exists(path_reduced_dir):  
        os.makedirs(path_reduced_dir)
        directories_created.append(path_reduced_dir)
    else:
        directories_already_existing.append(path_reduced_dir)   
        
    if not os.path.exists(path_reduced_star_pol_subtr_dir):  
        os.makedirs(path_reduced_star_pol_subtr_dir)
        directories_created.append(path_reduced_star_pol_subtr_dir)
    else:
        directories_already_existing.append(path_reduced_star_pol_subtr_dir)
    
    # Print which directories have been created and which already existed
    if any(directories_created):
        printandlog('\nThe following directories have been created:')
        for directory_sel in directories_created:
            printandlog('{0:s}'.format(directory_sel))
    
    if any(directories_already_existing):
        printandlog('\nThe following directories already exist. Data in these directories will be overwritten:')
        for directory_sel in directories_already_existing:
            printandlog('{0:s}'.format(directory_sel), wrap=False)
    
    return path_dark_files, path_flat_files, path_object_files, path_sky_files, path_center_files, \
           path_object_center_files, path_flux_files, path_sky_flux_files, \
           indices_to_remove_dark, indices_to_remove_flat, indices_to_remove_object, \
           indices_to_remove_sky, indices_to_remove_center, indices_to_remove_object_center, \
           indices_to_remove_flux, indices_to_remove_sky_flux, file_index_object, \
           file_index_flux, object_centering_method, combination_method_polarization   
    
###############################################################################
# read_fits_files
###############################################################################

def read_fits_files(path, silent=False):
    ''' 
    Read FITS-files
    
    Input:
        path: string or list of strings specifying paths to read FITS-files 
              and its headers from 
        silent: if True do not output message stating which file is read 
                (default = False)
    
    Output:
        data: image data cube (when 1 path specified) or list of image data 
              cubes (when more than 1 path specified) read from FITS-files
        header: headers (when 1 path specified) or list of headers 
                (when more than 1 path specified) read from FITS-files
    
    Note:
        If the data read is a single frame, an additional dimension is added to
        make it a 3D image data cube (third dimension of length 1).
                
    File written by Rob van Holstein
    Function status: verified
    '''
    
    if type(path) == str:
        # Turn string into list with single string
        path = [path]   

    # Read data and header        
    data = []
    header = []
    for path_sel in path:
        data_read, header_read = pyfits.getdata(path_sel, header=True)
        if data_read.ndim == 2:
            data_read = np.expand_dims(data_read, axis=0)
        elif data_read.ndim != 3:
            raise IOError('\n\n{0:s} is neither a cube nor a frame'.format(path_sel))
        data.append(data_read)
        header.append(header_read) 
        if silent is not True:
            printandlog('Read file ' + path_sel + '.', wrap=False)

    if len(data) == 1:
        # Remove data and header from list
        data = data[0]
        header = header[0]
        
    return data, header

###############################################################################
# write_fits_files
###############################################################################

def write_fits_files(data, path, header=False, silent=False):
    ''' 
    Write FITS-files
     
    Input:
        data: image data cube or list of image data cubes to write to FITS-files
        path: string or list of strings specifying paths to write FITS-files to
        header: headers or list of headers to write to FITS-files; if set to False
              no header will be written (default=False)
        silent: if True do not output message stating which file is written 
                (default = False)

    Note:
        If data, header and path are lists of unequal length, the redundant 
        list elements are ignored, meaning that files may not be written or
        that headers may not be used.
               
    File written by Rob van Holstein
    Function status: verified             
    '''

    # Prepare input to be processed
    if type(path) == str:
        path = [path]
    
    if type(data) == np.ndarray:
        data = [data]

    if type(header) == pyfits.header.Header:
        header = [header]

    if header == False:
        # Write FITS-files without headers
        for data_sel, path_sel in zip(data, path):
            hdu = pyfits.PrimaryHDU()
            hdu.data = data_sel.astype(np.float32)
            hdu.writeto(path_sel, overwrite=True, output_verify='silentfix+ignore')
            if silent is not True:
                printandlog('Wrote file ' + path_sel + '.', wrap=False)
    else:
        # Write FITS-files with headers
        for data_sel, header_sel, path_sel in zip(data, header, path):
            hdu = pyfits.PrimaryHDU()
            hdu.data = data_sel.astype(np.float32)
            hdu.header = header_sel
            hdu.writeto(path_sel, overwrite=True, output_verify='silentfix+ignore')
            if silent is not True:
                printandlog('Wrote file ' + path_sel + '.', wrap=False)

###############################################################################
# remove_bad_pixels
###############################################################################

def remove_bad_pixels(cube, frame_master_bpm, sigma_filtering=True):
    ''' 
    Remove bad pixels from an image cube or frame using the bad pixel map
    followed by optional repeated sigma-filtering
    
    Input:
        cube: image data cube or frame to filtered for bad pixels
        frame_master_bpm: frame indicating location of bad pixels with 0's and good
            pixels with 1's
        sigma_filtering: if True remove bad pixels remaining after applying
            master bad pixel map using repeated sigma-filtering (default = True)
    
    Output:
        cube_filtered: image data cube with bad pixels removed
    
    File written by Rob van Holstein
    Function status: verified
    '''

    # If the input is a frame turn it into a cube
    cube_ndim = cube.ndim
    if cube_ndim == 2:
        cube = np.expand_dims(cube, axis=0)
        
    # Define size of side of kernel for median filter
    filter_size_median = 5

    # Round filter size up to nearest odd number for a symmetric filter kernel
    filter_size_median = 2*(filter_size_median // 2) + 1
    
    # Remove bad pixels using the bad pixel map
    cube_median = ndimage.filters.median_filter(cube, size=(1, filter_size_median, \
                                                            filter_size_median))
    cube_filtered = cube_median + frame_master_bpm * (cube - cube_median)
        
    if sigma_filtering == True:
        # Define threshold factor for sigma filtering
        factor_threshold = 5
        
        # Define maximum number of iterations and counters for while-loop
        maximum_iterations = 10
        number_pixels_replaced = 1
        iteration_counter = 0
        
        # Prepare weights to compute mean without central pixel using convolution       
        filter_size = 7
        kernel = np.ones((1, filter_size, filter_size)) / (filter_size**2 - 9)
        for i in range(-1, 2):
            for j in range(-1, 2):
                kernel[0, filter_size//2 + i, filter_size//2 + j] = 0

        while number_pixels_replaced > 0 and iteration_counter < maximum_iterations:                       
            # Calculate local standard deviation using convolution
            cube_mean = ndimage.filters.convolve(cube_filtered, kernel)
            cube_EX2 = ndimage.filters.convolve(cube_filtered**2, kernel)
            cube_std = np.sqrt(cube_EX2 - cube_mean**2)
            
            # Compute threshold map for removal of pixels
            cube_threshold = factor_threshold*cube_std                             
            
            # Determine difference of image data with the local median
            cube_difference = np.abs(cube_filtered - cube_median)
            
            # Replace bad pixels by values in median filtered images
            pixels_to_replace = cube_difference > cube_threshold
            cube_filtered[pixels_to_replace] = cube_median[pixels_to_replace]
            
            # Compute number of pixels replaced and update number of iterations
            number_pixels_replaced = np.count_nonzero(pixels_to_replace)  
            iteration_counter += 1
    
    # If the input is a frame turn the resulting cube back into a frame
    if cube_ndim == 2:
        cube_filtered = np.squeeze(cube_filtered)
       
    return cube_filtered

###############################################################################
# process_sky_frames
###############################################################################

def process_sky_frames(path_sky_files, indices_to_remove_sky):
    '''
    Create a master sky frame from the SKY-files
    
    Input:
        path_sky_files: string or list of strings specifying paths of SKY-files
        indices_to_remove_sky: list of arrays with indices of frames to remove for each SKY-file

    Output:
        frame_master_sky: master sky frame
    
    File written by Rob van Holstein; based on a function by Christian Ginski
    Function status: verified       
    '''
      
    # Read sky frames
    list_cube_sky_raw = []
        
    for path_sel, indices_sel in zip(path_sky_files, indices_to_remove_sky):
        # Read data from file
        cube_sel = read_fits_files(path=path_sel, silent=True)[0]

        # Remove frames based on list of indices provided by the user
        cube_sel = np.delete(cube_sel, indices_sel, axis=0)

        # Append cube to list
        list_cube_sky_raw.append(cube_sel)
    
    # Make a single image cube from list of image cubes or frames
    cube_sky_raw = np.vstack(list_cube_sky_raw)

    # Compute median of sky frames
    frame_master_sky = np.median(cube_sky_raw, axis=0)
            
    printandlog('\nCreated master sky frame out of ' + str(len(path_sky_files)) + 
                ' raw SKY-file(s) comprising a total of ' + str(cube_sky_raw.shape[0]) + ' frame(s).')
    
    return frame_master_sky

###############################################################################
# create_bpm_darks
###############################################################################

def create_bpm_darks(list_frame_dark):
    '''
    Create a bad pixel map from DARK(,BACKGROUND)-files based on bias offsets and
    sigma filtering.
    
    Input:
        list_frame_dark: list of (mean-combined) DARK(,BACKGROUND)-frames
    
    Output:
        frame_bpm_dark: bad pixel map created from darks
        
    Function written by Rob van Holstein; constructed from functions by Christian Ginski
    Function status: verified
    '''

    # Create initial bad pixel map with only 1's
    frame_bpm_dark = np.ones(list_frame_dark[0].shape)

    for frame_dark in list_frame_dark:  
        # Remove outliers from dark frame and compute median and standard deviation    
        frame_dark_cleaned = sigmaclip(frame_dark, 5, 5)[0]
        stddev = np.nanstd(frame_dark_cleaned)
        median = np.nanmedian(frame_dark_cleaned)
    
        # Subtract median from dark frame and take absolute value
        frame_dark = np.abs(frame_dark - median)    
        
        # Initialize a bad pixel array with 1 as default pixel value
        frame_bpm = np.ones(frame_dark.shape)
        
        # Set pixels that deviate by more than 3.5 sigma from the frame median value 
        # to 0 to flag them as bad pixels
        frame_bpm[frame_dark > 3.5*stddev] = 0    
        
        # Add bad pixels found to master bad pixel map 
        frame_bpm_dark *= frame_bpm
        
    return frame_bpm_dark
        
###############################################################################
# create_bpm_flats
###############################################################################

def create_bpm_flats(list_frame_flat, list_exptime_flat):
    '''
    Create a bad pixel map from FLAT-files by finding all non-linear responding 
    pixels
    
    Input:
        list_frame_flat: list of (mean-combined) FLAT-frames
        list_exptime_flat: list of corresponding exposure times of flat
    
    Output:
        frame_bpm_flat: bad pixel map created from flats
        
    Function written by Christian Ginski; adapted by Rob van Holstein
    Function status: verified
    '''
    
    # Divide the flats with the longest and shortest exposure time
    comparison = list_frame_flat[-1] / list_frame_flat[0]

    # Compute the factor by which exposure time and thus pixel counts should have increased
    factor = list_exptime_flat[-1] / list_exptime_flat[0]
    
    # Set all pixels that were undefined to 1
    comparison[np.isinf(comparison)] = 1
    comparison[np.isnan(comparison)] = 1      
   
    # Create an array for the left and right detector side that does not contain strong outliers
    left = sigmaclip(comparison[15:1024, 36:933], 5, 5)[0]
    right = sigmaclip(comparison[5:1018, 1062:1958], 5, 5)[0]
    
    # Determine the standard deviation in the sigma filtered array
    stddev_left = np.nanstd(left)
    stddev_right = np.nanstd(right) 
    
    # Subtract the factor from the comparison array
    # The resulting array should contain 0 for all well responding pixels
    index_array_left = np.abs(comparison[15:1024, 36:933] - factor)
    index_array_right = np.abs(comparison[5:1018, 1062:1958] - factor)    
    
    # Identify pixels that deviate by more than 3.5 sigma from linear response
    mask_left = index_array_left > 3.5*stddev_left        
    mask_right = index_array_right > 3.5*stddev_right    

    # All good pixels are set to 1 and all non-linear responding pixels are set to 0
    badpix_left = np.ones(index_array_left.shape)
    badpix_left[mask_left] = 0

    badpix_right = np.ones(index_array_right.shape)
    badpix_right[mask_right] = 0
    
    # Final left and right side combined bad pixel mask is created
    frame_bpm_flat = np.ones(comparison.shape)
    frame_bpm_flat[15:1024, 36:933] = badpix_left
    frame_bpm_flat[5:1018, 1062:1958] = badpix_right

    # It is ensured that there are no NaN values in the final bad pixel mask
    frame_bpm_flat = np.nan_to_num(frame_bpm_flat)

    return frame_bpm_flat

###############################################################################
# process_dark_flat_frames
###############################################################################

def process_dark_flat_frames(path_dark_files, path_flat_files, indices_to_remove_dark, indices_to_remove_flat):
    '''
    Process DARK(,BACKGROUND)- and FLAT-files to create a master flat frame and 
    a bad pix map. The number of dark and flat frames provided must be the same, 
    and they must have matching exposure times. Generally a sequence of darks and 
    flats with exposure times 1, 2, 3, 4, 5 s or 2, 4, 6, 8, 10 s is used. The bad 
    pixel mask contains flags for all pixels that have a strong bias offset or that 
    respond non-linearly.
    
    Input:
        path_dark_files: list of paths to raw DARK(,BACKGROUND)-files
        path_flat_files: list of paths to raw FLAT-files
        indices_to_remove_dark: list of 1-D arrays with indices of frames to 
            be removed for each DARK(,BACKGROUND)-file. If no frames are to be 
            removed the array is empty. 
        indices_to_remove_flat: list of 1-D arrays with indices of frames to 
            be removed for each FLAT-file. If no frames are to be removed the 
            array is empty.         
        
    Output:
        frame_master_flat: master flat frame
        frame_master_bpm: master bad pixel map (1 indicates good pixel; 
                                                0 indicates bad pixel)
    
    File written by Christian Ginski; adapted by Rob van Holstein
    Function status: verified
    '''
  
    # Read dark and flat frames, compute mean and prepare flats
    list_frame_dark = []       
    list_frame_flat = []
    list_frame_flat_norm = []
    list_exptime_flat = []

    for path_dark, path_flat, indices_dark, indices_flat in zip(path_dark_files, path_flat_files, indices_to_remove_dark, indices_to_remove_flat):
        # Read data from files
        cube_dark = read_fits_files(path=path_dark, silent=True)[0]
        cube_flat, header_flat = read_fits_files(path=path_flat, silent=True)

        # Remove frames based on list of indices provided by the user
        cube_dark = np.delete(cube_dark, indices_dark, axis=0)
        cube_flat = np.delete(cube_flat, indices_flat, axis=0)

        # Compute mean of cubes
        frame_dark = np.mean(cube_dark, axis=0)
        frame_flat = np.mean(cube_flat, axis=0)
        
        # Dark-subtract flat frame
        frame_flat_dark_subtracted = frame_flat - frame_dark

        # Filter dark-subtracted flat for zeros and NaN's
        frame_flat_dark_subtracted = np.nan_to_num(frame_flat_dark_subtracted)
        frame_flat_dark_subtracted[frame_flat_dark_subtracted <= 0] = 1
        
        # Determine exposure time of each flat file and normalize flat with it
        exptime = header_flat['EXPTIME']
        frame_flat_norm = frame_flat_dark_subtracted / exptime

        # Append frames and exposure time to lists
        list_frame_dark.append(frame_dark)
        list_frame_flat.append(frame_flat_dark_subtracted)      
        list_frame_flat_norm.append(frame_flat_norm)
        list_exptime_flat.append(exptime)

    # Median combine all normalized flats
    frame_flat = np.median(np.array(list_frame_flat_norm), axis=0)
    
    # Select the left and right detector area that actualy receives signal
    frame_flat_left = frame_flat[15:1024, 36:933] 
    frame_flat_right = frame_flat[5:1018, 1062:1958]

    # Normalize left and right side of the flat with the respective median values
    frame_flat_left_norm = frame_flat_left / np.median(frame_flat_left)
    frame_flat_right_norm = frame_flat_right / np.median(frame_flat_right)
    
    # Create a baseline flat image that encompases left and right detector side
    # All values are 1, i.e. if applied this baseline flat does not influence the reduction
    # This is mainly to prevent later edge areas from containing blown-up data values
    frame_master_flat = np.ones(frame_flat.shape)
    
    # Construct the final full detector flat
    frame_master_flat[15:1024, 36:933] = frame_flat_left_norm 
    frame_master_flat[5:1018, 1062:1958] = frame_flat_right_norm    

    # Create a bad pixel map based on dark and flat frames 
    # Theoretically all detector artefacts should be identified this way
    # Create a list of bad pixel maps with pixels with a bias offset flagged as bad
    frame_bpm_dark = create_bpm_darks(list_frame_dark)

    # Createa a bad pixel map with non-linear responding pixels flagged as bad
    frame_bpm_flat = create_bpm_flats(list_frame_flat, list_exptime_flat)

    # Create the master bad pixel map by combing the one from the darks and flats    
    frame_master_bpm = frame_bpm_dark * frame_bpm_flat

    printandlog('\nCreated master flat frame and bad pixel map out of ' + str(len(path_flat_files)) + 
                ' pairs of raw FLAT- and DARK(,BACKGROUND)-file(s).')

    return frame_master_flat, frame_master_bpm

###############################################################################
# compute_fwhm_separation
###############################################################################

def compute_fwhm_separation(filter_used):
    '''
    Compute theoretical full width half maximum (FWHM) of point spread function
    of IRDIS and separation of satellite spots in CENTER-files
    
    Input:
        filter_used: string of filter from header 'ESO INS1 FILT ID'
        
    Output:
        fwhm: full width half maximum (pixels)
        separation: separation of satellite spots (pixels)
    
    File written by Rob van Holstein; based on a function by Julien Milli
    Function status: verified
    '''

    # Define diameter of telescope primary mirror M1 (m)
    D_telescope = 8.2
    
    # Define approximate separation of satellite spots (pixels) and filter central wavelength and bandwidth (m)
    if filter_used in ['FILT_BBF_Y', 'FILT_NBF_HeI']:
        separation = 31.0
        lambda_c = 1043e-9
        delta_lambda = 140e-9
    elif filter_used in ['FILT_BBF_J', 'FILT_NBF_ContJ', 'FILT_NBF_PaB']:
        separation = 37.4
        lambda_c = 1245e-9
        delta_lambda = 240e-9
    elif filter_used in ['FILT_BBF_H', 'FILT_NBF_ContH', 'FILT_NBF_FeII']:
        separation = 48.5
        lambda_c = 1625e-9
        delta_lambda = 290e-9
    elif filter_used in ['FILT_BBF_Ks', 'FILT_NBF_ContK1', 'FILT_NBF_H2', 'FILT_NBF_BrG', 'FILT_NBF_CntK2', 'FILT_NBF_CO']:
        separation = 64.5
        lambda_c = 2182e-9
        delta_lambda = 300e-9

    # Compute theoretical full width half maximum
    fwhm = np.rad2deg((lambda_c + 0.5*delta_lambda) / D_telescope) * 3600 / pixel_scale

    return fwhm, separation

###############################################################################
# process_center_frames
###############################################################################

def process_center_frames(path_center_files, 
                          indices_to_remove_center, 
                          path_object_center_files, 
                          indices_to_remove_object_center, 
                          frame_master_flat, 
                          frame_master_bpm, 
                          frame_master_sky, 
                          center_subtract_object=True, 
                          center_coordinates=(477, 521, 1503, 511), 
                          sigma_filtering=True):
    '''
    Process the CENTER frames by subtracting the background, flat-fielding, 
    removing bad pixels and computing the mean over the NDIT's
    
    Input:
        path_center_files: list of paths to raw CENTER-files
        indices_to_remove_center: list of arrays with indices of frames to remove
            for each CENTER-file
        path_object_center_files: list of paths to raw OBJECT-files to be 
            subtracted from raw CENTER-files
        indices_to_remove_object_center: list of arrays with indices of frames to
            remove for each OBJECT-file to be subtracted from raw CENTER-files           
        frame_master_flat: master flat frame
        frame_master_bpm: frame indicating location of bad pixels with 0's and good
            pixels with 1's
        frame_master_sky: master sky frame for OBJECT- (or CENTER)-files
        center_subtract_object: if True subtract the OBJECT-file taken
            closest in time from the CENTER-file (default = True). This generally
            results in a more accurate determination of the center coordinates
            as the background and any other celestial objects in the field of view
            are suppressed. When the difference between the image orientations of 
            the CENTER- and OBJECT-frames is large (i.e. difference in derotator 
            position angle for field-tracking and difference in parallactic angle 
            for pupil-tracking), the OBJECT-frame will be rotated around the initial 
            guess of the centers as defined by center_coordinates before 
            subtracting it from the CENTER-file.
        center_coordinates: length-4-tuple with initial guess of center coordinates
            x_left: initial guess of x-coordinate of center of left frame half
            y_left: initial guess of y-coordinate of center of left frame half
            x_right: initial guess of x-coordinate of center of right frame half
            y_right: initial guess of y-coordinate of center of right frame half
            Note that the center coordinates are defined in the complete frame, 
            i.e. with both detector halves (pixels; 0-based). The default value 
            is (477, 521, 1503, 511). 
        sigma_filtering: if True remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True)
       
    Output:
        list_frame_center_processed: list of processed center frames
        header: list of headers of center files
                
    File written by Rob van Holstein
    Function status: verified
    '''

    # Print if and how background is subtracted
    if center_subtract_object == True:
        printandlog('\nSubtracting OBJECT-file closest in time from CENTER-file(s).')
    else:
        if np.array_equal(frame_master_sky, np.zeros((1024, 2048))):
            printandlog('\nNot subtracting background from CENTER-file(s).')
        else:
            printandlog('\nSubtracting master sky from CENTER-file(s).')
    
    # Determine tracking mode used
    tracking_mode_used = pyfits.getheader(path_center_files[0])['ESO INS4 COMB ROT']

    # Determine filter and tracking mode used
    filter_used = pyfits.getheader(path_center_files[0])['ESO INS1 FILT ID']
                 
    # Create empty lists to store processed images and headers in
    list_frame_center_processed = []
    header = []
    
    for i, (path_center_sel, indices_center_sel, path_object_sel, indices_object_sel) in enumerate(
            zip(path_center_files, indices_to_remove_center, path_object_center_files, indices_to_remove_object_center)):
        # Read data and header from file
        cube_center, header_center = read_fits_files(path=path_center_sel, silent=True)
    
        # Remove frames based on list of indices provided by the user
        cube_center = np.delete(cube_center, indices_center_sel, axis=0)
   
        if center_subtract_object == True:
            # Use OBJECT-file closest in time to CENTER-file as backround
            cube_object, header_object = read_fits_files(path=path_object_sel, silent=True)
            
            # Remove frames based on list of indices provided by the user
            cube_object = np.delete(cube_object, indices_object_sel, axis=0)
            
            # Compute mean of frames and shift for dithering
            frame_background = np.mean(cube_object, axis=0)
            shift_x_dith = -header_object['ESO INS1 DITH POSX']
            shift_y_dith = -header_object['ESO INS1 DITH POSY']
            frame_background = ndimage.shift(frame_background, [shift_y_dith, shift_x_dith], order=3, mode='constant', cval=0.0, prefilter=True)   

            # Compute theoretical full width half maximum and separation of satellite spots
            fwhm, separation = compute_fwhm_separation(filter_used)
            
            # Compute maximum allowable image rotation before rotating the object-frame (deg)
            max_image_rotation = np.rad2deg(0.3*fwhm / separation)
            
            if tracking_mode_used == 'FIELD':
                # Determine derotator image position angle of center- and object-frame
                rot_angle_center = np.deg2rad(header_center['ESO INS4 DROT2 POSANG'])
                rot_angle_object = np.deg2rad(header_object['ESO INS4 DROT2 POSANG'])
            
            elif tracking_mode_used == 'PUPIL':
                # Determine parallactic angle of center- and object-frame
                rot_angle_center = np.deg2rad(compute_mean_angle([header_center['ESO TEL PARANG START'], header_center['ESO TEL PARANG END']]))           
                rot_angle_object = np.deg2rad(compute_mean_angle([header_object['ESO TEL PARANG START'], header_object['ESO TEL PARANG END']]))
            
            # Determine difference between image rotation angles of center- and object-frame (deg)
            delta_rot_angle = np.rad2deg(np.arctan2(np.sin(rot_angle_center - rot_angle_object), np.cos(rot_angle_center - rot_angle_object)))
            
            if np.abs(delta_rot_angle) > max_image_rotation:
                # Remove bad pixels from background frame and cut into two halves
                center_coordinates_print = tuple(x + 1 for x in center_coordinates)
                printandlog('\nRotating OBJECT-frame around the initial guess of center before subtracting it from the CENTER-file.')
                printandlog('center_coordinates = ' + str(center_coordinates_print))
                printandlog('')
                frame_background = remove_bad_pixels(cube=frame_background, frame_master_bpm=frame_master_bpm, sigma_filtering=sigma_filtering)
                frame_background_left = frame_background[:, :1024]
                frame_background_right = frame_background[:, 1024:]
                
                # Rotate frame halves about estimated centers and append them together again
                frame_background_left = rotateskimage(frame_background_left, delta_rot_angle, center=(center_coordinates[0], center_coordinates[1]), order=3)
                frame_background_right = rotateskimage(frame_background_right, delta_rot_angle, center=(center_coordinates[2] - 1024, center_coordinates[3]), order=3)
                frame_background = np.append(frame_background_left, frame_background_right, axis=1)

        else:
            # Use the master sky frame as background
            frame_background = frame_master_sky 

        # Subtract background and divide by master flat
        cube_bgsubtr_flatfielded = (cube_center - frame_background) / frame_master_flat
        
        # Remove bad pixels of each frame
        cube_badpixel_filtered = remove_bad_pixels(cube=cube_bgsubtr_flatfielded, frame_master_bpm=frame_master_bpm, sigma_filtering=sigma_filtering)
                
        # Compute mean over NDIT frames
        frame_mean = np.mean(cube_badpixel_filtered, axis=0)

        # Append reduced images and headers
        list_frame_center_processed.append(frame_mean)
        header.append(header_center)
        
        # Print which file has been processed
        if i == 0:
            printandlog('')
        printandlog('Processed file ' + str(i + 1) + '/' + str(len(path_center_files)) + ': {0:s}'.format(os.path.basename(path_center_sel)))

    return list_frame_center_processed, header

###############################################################################
# fit_2d_gaussian
###############################################################################
    
def fit_2d_gaussian(frame, x0=None, y0=None, x_stddev=1.0, y_stddev=1.0, theta=0.0, 
                    crop_radius=None, sigfactor=None, saturation_level=None):
    '''
    Fit a 2D-Gaussian to an image frame and return the center coordinates
    
    Input:
        frame: image frame to fit Gaussian to
        x0: center x-coordinate of sub-image used to fit Gaussian (pixels). Must be
            integer. If None, the center of the frame is used. Ignored when
            crop_radius is None (default = None).
        y0: center y-coordinate of sub-image used to fit Gaussian (pixels). Must be
            integer. If None, the center of the frame is used. Ignored when
            crop_radius is None (default = None).
        x_stddev: standard deviation in x-direction of Gaussian before rotating
            by theta (default = 1.0)
        y_stddev: standard deviation in y-direction of Gaussian before rotating
            by theta (default = 1.0)
        theta: rotation angle of Gaussian in rad, positive clockwise
            (default = 0.0)
        crop_radius: half the length of side of square cropped sub-images used 
            to fit Gaussian to (pixels). Must be integer. If None, the complete 
            frame is used for the fitting and the values of x0 and y0 are 
            ignored (default = None).
        sigfactor: all sub-image pixels with values smaller than 
            sigfactor*standard deviation are replaced by random Gaussian noise 
            to mask them for fitting the 2D Gaussian. If None, no pixels are
            replaced by Gaussian noise (default = None).
        saturation_level: all pixels within the smallest circle encompassing 
            the pixels with a value equal to or higher than saturation_level 
            are ignored when fitting the 2D Gaussian. We use a circle because
            strongly saturated pixels in the peak of the PSF often have values 
            lower than the saturation level. If None, no pixels are ignored 
            (default = None).
  
    Output:
        x_fit: x-coordinate of center of fitted Gaussian in complete frame (pixels)
        y_fit: y-coordinate of center of fitted Gaussian in complete frame (pixels)
        x_fit_sub_image: x-coordinate of center of fitted Gaussian in sub-image frame (pixels)
        y_fit_sub_image: y-coordinate of center of fitted Gaussian in sub-image frame (pixels)
        sub_image: cropped sub-image frame used to fit Gaussian to
                
    File written by Rob van Holstein; based on a function by Julien Milli
    Function status: verified
    '''

    # Check input and raise errors
    if frame.ndim != 2:
        raise ValueError('\n\nframe should be a 2-dimensional array.')
    if type(crop_radius) != int and crop_radius is not None:
        raise TypeError('\n\ncrop_radius must be integer or None.') 
    if type(x0) not in [int, np.int32, np.int64] and y0 is not None:
        raise TypeError('\n\nx0 must be integer or None.')
    if type(y0) not in [int, np.int32, np.int64] and x0 is not None:
        raise TypeError('\n\ny0 must be integer or None.')

    # Set approximate center coordinates to center of frame if not specified
    if x0 is None:
        x0 = frame.shape[-1] // 2
    if y0 is None:
        y0 = frame.shape[-2] // 2
    
    # Set saturation level to infinity if not specified
    if saturation_level is None:
        saturation_level = np.inf

    if crop_radius is None:
        # Use complete frame
        sub_image = np.copy(frame)
    else:
        # Cut out sub-image around center
        sub_image = np.copy(frame[y0 - crop_radius:y0 + crop_radius, \
                          x0 - crop_radius:x0 + crop_radius])
                
    # Create grid of x- and y-coordinates  
    y, x = np.mgrid[:sub_image.shape[-2], :sub_image.shape[-1]]
    
    # Determine position of maximum value in sub-image
    max_sub_image = np.max(sub_image)
    y_max, x_max = np.where(sub_image == max_sub_image)

    # If the maximum value is found on multiple pixels, take the position of first pixel
    if len(x_max) > 1:
        x_max = x_max[0]
        y_max = y_max[0]
  
    if sigfactor is not None:
        # Replace all values smaller than sigfactor*standard deviation by random Gaussian noise
        clipmed, clipstd = sigma_clipped_stats(sub_image, sigma=2)[1:]
        indices_background = np.where(sub_image <= clipmed + sigfactor*clipstd)
        subimnoise = np.random.randn(sub_image.shape[0], sub_image.shape[1])*clipstd
        sub_image[indices_background] = subimnoise[indices_background]

    # Determine accurate coordinates of center by fitting a 2D Gaussian              
    p_init1 = models.Gaussian2D(amplitude=max_sub_image, x_mean=int(x_max), y_mean=int(y_max), \
                                x_stddev=x_stddev, y_stddev=y_stddev, theta=theta)
    fit_p = fitting.LevMarLSQFitter()
    p1 = fit_p(p_init1, x, y, sub_image, maxiter=1000, acc=1e-08)
    x_fit_sub_image = p1.x_mean[0]
    y_fit_sub_image = p1.y_mean[0]

    # Compute radial distances from center of saturated pixels
    radius = np.sqrt((x - x_fit_sub_image)**2 + (y - y_fit_sub_image)**2)
    radius_saturation = radius[sub_image >= saturation_level]
    
    if any(radius_saturation) == True:
        # Mask all pixels within the smallest centered circle encompassing all saturated pixels
        mask_saturated = radius <= np.max(radius_saturation)
    
        # Remove saturated pixels
        indices_not_saturated = np.where(~mask_saturated)
        x = x[indices_not_saturated]
        y = y[indices_not_saturated]
        sub_image_masked = sub_image[indices_not_saturated]
    
        # Set saturated pixels to NaN in sub-image for plotting
        sub_image[mask_saturated] = np.nan
    
        # Determine accurate coordinates of center by fitting a 2D Gaussian to the sub-image without saturated pixels              
        p_init2 = models.Gaussian2D(amplitude=max_sub_image, x_mean=x_fit_sub_image, y_mean=y_fit_sub_image, \
                 x_stddev=x_stddev, y_stddev=y_stddev, theta=theta)
        p2 = fit_p(p_init2, x, y, sub_image_masked, maxiter=1000, acc=1e-08)
        x_fit_sub_image = p2.x_mean[0]
        y_fit_sub_image = p2.y_mean[0]
   
    # Compute spot coordinates in coordinates of complete image
    if crop_radius is None:      
        x_fit = x_fit_sub_image
        y_fit = y_fit_sub_image
    else:
        x_fit = x0 - crop_radius + x_fit_sub_image
        y_fit = y0 - crop_radius + y_fit_sub_image

    return x_fit, y_fit, x_fit_sub_image, y_fit_sub_image, sub_image

###############################################################################
# find_center_coordinates
###############################################################################
    
def find_center_coordinates(list_frame_center_processed, 
                            path_processed_center_files, 
                            center_coordinates=(477, 521, 1503, 511), 
                            param_centering=(12, None, 30000)):
    
    '''
    Find coordinates of star center from processed CENTER frames. The function
    shows the fitted coordinates of the satellite spots in an image and writes
    a REG-file that indicates the fitted centers and that can be loaded as a 
    region in the FITS-files of the processed center frames.
    
    Input:
        list_frame_center_processed: list of processed center frames
        path_processed_center_files: list of paths to processed CENTER files
        center_coordinates: length-4-tuple with initial guess of center coordinates
            x_left: initial guess of x-coordinate of center of left frame half
            y_left: initial guess of y-coordinate of center of left frame half
            x_right: initial guess of x-coordinate of center of right frame half
            y_right: initial guess of y-coordinate of center of right frame half
            Note that the center coordinates are defined in the complete frame, 
            i.e. with both detector halves (pixels; 0-based). The default value 
            is (477, 521, 1503, 511). 
        param_centering: length-3-tuple with parameters for 2D 
            Gaussian fitting the centers of the satellite spots in the CENTER-frames:
            crop_radius: half the length of side of square cropped sub-images used 
                to fit 2D Gaussian to (pixels). Must be integer. If None, the complete 
                frame is used for the fitting and center_coordinates is ignored.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of param_centering is (12, None, 30000).
                
    Output:
        center_coordinates: a length-4-tuple with the determined center 
        coordinates (pixels): (x_left, y_left, x_right, y_right)
                
    File written by Rob van Holstein; based on functions by Julien Milli and 
    Christian Ginski
    Function status: verified
    '''
    
    # Read centering parameters
    x_center_0 = np.array([center_coordinates[0], center_coordinates[2]])
    y_center_0 = np.array([center_coordinates[1], center_coordinates[3]])
    crop_radius, sigfactor, saturation_level = param_centering
    
    # Subtract 1024 from the right x-coordinate to make it valid for a frame half
    x_center_0[1] -= 1024
   
    # Read headers
    header = [pyfits.getheader(x) for x in path_processed_center_files]
    
    # Determine filter used
    filter_used = header[0]['ESO INS1 FILT ID'] 
                
    # Compute theoretical full width half maximum and separation of satellite spots
    fwhm, separation = compute_fwhm_separation(filter_used)
    
    # Set NaN-vales in color map to gray 
    cmap = mpl.cm.viridis
    cmap.set_bad(color='gray')

    # Create zero arrays to save fitted center coordinates in
    x_center_fit = np.zeros((len(list_frame_center_processed), 2))
    y_center_fit = np.copy(x_center_fit)  
            
    # For each frame
    for i, (frame_sel, header_sel, path_sel) in enumerate(zip(list_frame_center_processed, header, path_processed_center_files)):       
        # Cut frame into a left and right half
        frame_halves = [frame_sel[:, :1024], frame_sel[:, 1024:]]

        # Determine waffle pattern used ('x' or '+')
        waffle_pattern = header_sel['ESO OCS WAFFLE ORIENT']
        
        # Define approximate position angles of spots
        if waffle_pattern == 'x':
            position_angle = np.deg2rad(np.array([45., 135., 225., 315.]))
            titles_sub_image = ['upper right', 'upper left', 'lower left', 'lower right']
        else:
            printandlog('\nWARNING, waffle pattern is \'+\'. Correct functioning of centering not tested.')
            position_angle = np.deg2rad(np.array([0., 90., 180., 270.]))
            titles_sub_image = ['right', 'top', 'left', 'bottom']
        
        # Initiate plot to show sub-images used to fit coordinates of satellite spots
        path_plot_sub_images = os.path.splitext(path_sel)[0] + '.png'
        printandlog('\nCreating plot ' + path_plot_sub_images + ' showing sub-images of satellite spots with fitted coordinates.')
        fig, axs = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True, subplot_kw={'xticks': [], 'yticks': []}, figsize=(8, 4.3)) 
        center_coordinates_print = tuple(x + 1 for x in center_coordinates)
        title_main = 'object_center_coordinates = %s' % (center_coordinates_print,)
        if center_coordinates == (477, 521, 1503, 511):
            title_main += ' (default)'
        title_main += '\ncenter_param_centering = %s' % (param_centering,)         
        if param_centering == (12, None, 30000):
            title_main += ' (default)'
        fig.suptitle(title_main, horizontalalignment='center')
        fig.subplots_adjust(top=0.7)
        sub_image_position = [[(0, 1), (0, 0), (1, 0), (1, 1)], [(0, 3), (0, 2), (1, 2), (1, 3)]]
        
        # Create zero arrays to save fitted spot coordinates in 
        x_fit = np.zeros((2, 4))
        y_fit = np.zeros((2, 4))
                
        # For each frame half
        for j, (frame_half, x_center_0_sel, y_center_0_sel) in enumerate(zip(frame_halves, x_center_0, y_center_0)): 
            # For each satellite spot
            for k, angle in enumerate(zip(position_angle)):
                # Determine approximate coordinates of satellite spots         
                x_spot_0 = int(np.round(x_center_0_sel + np.cos(angle) * separation))
                y_spot_0 = int(np.round(y_center_0_sel + np.sin(angle) * separation))

                # Determine accurate coordinates of satellite spots by fitting a 2D Gaussian           
                x_fit[j, k], y_fit[j, k], x_fit_sub_image, y_fit_sub_image, sub_image = fit_2d_gaussian(frame=frame_half, x0=x_spot_0, y0=y_spot_0, \
                x_stddev=1.1*fwhm, y_stddev=0.63*fwhm, theta=angle, crop_radius=crop_radius, sigfactor=sigfactor, saturation_level=saturation_level)
                
                # Plot sub-images used
                axs[sub_image_position[j][k]].imshow(sub_image, cmap=cmap, origin='lower',interpolation='nearest')
                axs[sub_image_position[j][k]].plot(x_fit_sub_image, y_fit_sub_image, 'rx', markersize=60/crop_radius)
                axs[sub_image_position[j][k]].set_title(titles_sub_image[k])

        # Add 1024 to the x-coordinates of the satellite spots on the right frame half to make the values valid for complete frame    
        x_fit[1, :] += 1024         
        
        # Compute gradients and constants of lines between opposite spots
        a13 = (y_fit[:, 0] - y_fit[:, 2]) / (x_fit[:, 0] - x_fit[:, 2])          
        b13 = y_fit[:, 0] - a13*x_fit[:, 0]
        a24 = (y_fit[:, 1] - y_fit[:, 3]) / (x_fit[:, 1] - x_fit[:, 3])
        b24 = y_fit[:, 1] - a24*x_fit[:, 1]

        # Compute center coordinates from intersection of lines
        x_center_fit[i, :] = (b24 - b13) / (a13 - a24)
        y_center_fit[i, :] = a13*x_center_fit[i, :] + b13
                        
        # Make titles for left and right half and show sub-images of satellite spots
        ext = []
        for j in range(4):
            ext.append([axs[0, j].get_window_extent().x0, axs[0, j].get_window_extent().width])
        inv = fig.transFigure.inverted()
        width_left = ext[0][0] + (ext[1][0] + ext[1][1] - ext[0][0]) / 2
        left_center = inv.transform((width_left, 1))
        width_right = ext[2][0] + (ext[3][0] + ext[3][1] - ext[2][0]) / 2
        right_center = inv.transform((width_right, 1))
        plt.figtext(left_center[0], 0.81, 'Left frame half', va='center', ha='center', size=12)
        plt.figtext(right_center[0], 0.81, 'Right frame half', va='center', ha='center', size=12)
        plt.savefig(path_plot_sub_images, dpi = 300, bbox_inches = 'tight')
        plt.close(fig)
               
        # Write REG-file with fitted lines and circles indicating centers to use with FITS-files of processed center frames
        path_reg_file = os.path.splitext(path_sel)[0] + '.reg'
        f = open(path_reg_file, 'w')
        f.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=0 delete=1 include=1 source=1\n")
        f.write("physical\n")
        f.write('circle({0:8.2f}, {1:8.2f}, 2.5)\n'.format(x_center_fit[i, 0] + 1, y_center_fit[i, 0] + 1))
        f.write('line({0:8.2f}, {1:8.2f}, {2:8.2f}, {3:8.2f})\n'.format(x_fit[0, 0] + 1, y_fit[0, 0] + 1, x_fit[0, 2] + 1, y_fit[0, 2] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 1.0)\n'.format(x_fit[0, 0] + 1, y_fit[0, 0] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 1.0)\n'.format(x_fit[0, 2] + 1, y_fit[0, 2] + 1))
        f.write('line({0:8.2f}, {1:8.2f}, {2:8.2f}, {3:8.2f})\n'.format(x_fit[0, 1] + 1, y_fit[0, 1] + 1, x_fit[0, 3] + 1, y_fit[0, 3] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 1.0)\n'.format(x_fit[0, 1] + 1, y_fit[0, 1] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 1.0)\n'.format(x_fit[0, 3] + 1, y_fit[0, 3] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 2.5)\n'.format(x_center_fit[i, 1] + 1, y_center_fit[i, 1] + 1))
        f.write('line({0:8.2f}, {1:8.2f}, {2:8.2f}, {3:8.2f})\n'.format(x_fit[1, 0] + 1, y_fit[1, 0] + 1, x_fit[1, 2] + 1, y_fit[1, 2] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 1.0)\n'.format(x_fit[1, 0] + 1, y_fit[1, 0] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 1.0)\n'.format(x_fit[1, 2] + 1, y_fit[1, 2] + 1))
        f.write('line({0:8.2f}, {1:8.2f}, {2:8.2f}, {3:8.2f})\n'.format(x_fit[1, 1] + 1, y_fit[1, 1] + 1, x_fit[1, 3] + 1, y_fit[1, 3] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 1.0)\n'.format(x_fit[1, 1] + 1, y_fit[1, 1] + 1))
        f.write('circle({0:8.2f}, {1:8.2f}, 1.0)\n'.format(x_fit[1, 3] + 1, y_fit[1, 3] + 1))        
        f.close()  
        printandlog('\nWrote file ' + path_reg_file + ' which can be loaded as a region in DS9 (in the top bar: Region --> Load Regions) and shows the fitted coordinates of the satellite spots and the center.')

    # Print center coordinates found
    max_path_length = max([len(os.path.basename(x)) for x in path_processed_center_files])
    separator_top = ' '*(max_path_length - 9)
    printandlog('\nMeasured center coordinates per file (pixels):')
    printandlog('File name' + separator_top + '    x_left    y_left    x_right    y_right', wrap=False)
    for i, path_sel in enumerate(path_processed_center_files):
        separator = ' '*(max_path_length - len(os.path.basename(path_sel)))
        printandlog(os.path.basename(path_sel) + separator + '    %.2f    %.2f    %.2f    %.2f' 
                    % (x_center_fit[i, 0] + 1, y_center_fit[i, 0] + 1, x_center_fit[i, 1] + 1, y_center_fit[i, 1] + 1), wrap=False)    
    
    # Compute mean of fitted center coordinates    
    x_center = np.mean(x_center_fit, axis=0)
    y_center = np.mean(y_center_fit, axis=0)

    if len(list_frame_center_processed) > 1:
        # Compute standard deviation of fitted center coordinates
        x_center_std = np.std(x_center_fit, ddof=1, axis=0)
        y_center_std = np.std(y_center_fit, ddof=1, axis=0)

        # Print mean center coordinates with error
        separator = ' '*max([len('%.2f' % x) for x in np.append(x_center_std, y_center_std)])
        printandlog('\nFinal center coordinates (pixels):')
        printandlog('x_left' + separator + '         y_left' + separator + '         x_right' + separator + '         y_right', wrap=False)
        printandlog('%.2f +/- %.2f    %.2f +/- %.2f    %.2f +/- %.2f    %.2f +/- %.2f' 
                    % (x_center[0] + 1, x_center_std[0], y_center[0] + 1, y_center_std[0], x_center[1] + 1, x_center_std[1], y_center[1] + 1, y_center_std[1]), wrap=False)    
        
        # Print warning if there is significant deviation among the center coordinates found
        if any(np.append(x_center_std, y_center_std) > 0.5):
            printandlog('\nWARNING, for at least one of the fitted center coordinates the standard deviation is larger than 0.5 pixels.')
    else:
        # Print mean center coordinates without error
        printandlog('\nFinal center coordinates (pixels):')
        printandlog('x_left    y_left    x_right    y_right')
        printandlog('%.2f    %.2f    %.2f    %.2f' % (x_center[0] + 1, y_center[0] + 1, x_center[1] + 1, y_center[1] + 1))    

    #TODO: plot fitted center coordinates vs time
    #TODO: please rename the function
    # Define function to plot determined center coordinates
    def plot_center_coordinates(data, x_y, left_right):
        font_size = 10
        width_figure = 4.7 / 20 * len(data)
        if width_figure < 4.7:
            width_figure = 4.7   
        elif width_figure > 30:
            width_figure = 40                  
        plot_name = name_file_root + 'center_coordinates_' + x_y + '_' + left_right + '.png'            
        path_plot = os.path.join(path_center_dir, plot_name)
        printandlog(path_plot, wrap=False)
        image_number = np.arange(1, len(data)+1)
        plt.figure(figsize = (width_figure, 3.0))
        plt.plot(image_number, data, '-ok')
        ax = plt.gca()
        ax.set_xlabel(r'Image', fontsize=font_size)
        ax.tick_params(axis = 'x', labelsize=font_size)
        ax.set_xlim([0, max(image_number) + 1])
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        ax.set_ylabel(x_y + ' ' + left_right + ' (pixels)', fontsize=font_size)
        ax.tick_params(axis='y', labelsize=font_size)  
        ax.ticklabel_format(useOffset=False, axis='y')
        ax.grid()
        plt.tight_layout()
        plt.savefig(path_plot, dpi=300, bbox_inches='tight')
        plt.close()

    # Assemble output
    center_coordinates = (x_center[0], y_center[0], x_center[1], y_center[1])
    
    return center_coordinates

###############################################################################
# create_sub_image
###############################################################################
    
def create_sub_image(frame, x0, y0, crop_radius):
    '''
    Create a square sub-image by cropping it from an image frame
    
    Input:
        frame: image frame to be cropped
        x0: x-coordinate of center of sub-image
        y0: y-coordinate of center of sub-image
        crop_radius: half the length of the sides of the square sub-image. If
            None the complete frame is returned as the sub-image.
       
    Output:
        sub_image: sub-image
                
    File written by Rob van Holstein
    Function status: verified
    '''
    
    if crop_radius is None:
        # Use complete frame
        sub_image = np.copy(frame)
    else:
        # Cut out sub-image around center
        sub_image = np.copy(frame[y0 - crop_radius:y0 + crop_radius, \
                          x0 - crop_radius:x0 + crop_radius])
    return sub_image

###############################################################################
# process_object_frames
###############################################################################

def process_object_frames(path_object_files, 
                          file_index_object, 
                          indices_to_remove_object, 
                          frame_master_flat, 
                          frame_master_bpm, 
                          frame_master_sky, 
                          sigma_filtering=True, 
                          centering_method='center frames', 
                          center_coordinates=(477, 521, 1503, 511), 
                          param_centering=(60, None, 30000), 
                          collapse_ndit=False, 
                          show_images_center_coordinates=True):
    '''
    Process the OBJECT frames by subtracting the background, flat-fielding, 
    removing bad pixels, centering, computing the mean over the NDIT's and
    computing the single sum and difference images. 
    
    Input:
        path_object_files: list of paths to raw OBJECT-files
        file_index_object: list of file indices of OBJECT-files (0-based)
        indices_to_remove_object: list of arrays with indices of frames to remove for each OBJECT-file
        frame_master_flat: master flat frame
        frame_master_bpm: frame indicating location of bad pixels with 0's and 
            good pixels with 1's
        frame_master_sky: master sky frame for OBJECT-files
        sigma_filtering: if True remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True)
        centering_method: method to center the images. If 'center frames' or
            'manual', use fixed coordinates as provided by center_coordinates. 
            If 'gaussian', fit a 2D Gaussian to each frame. If 
            'cross-correlation', fit a 2D Gaussian to the first frame and then
            use cross-correlation to align (register) the other frames onto the 
            centered first  frame. For 'gaussian' and 'cross-correlation' 
            center_coordinates is used as initial guess of the center 
            coordinates and the determined center coordinates are plotted for
            each image (default = 'center frames').
        center_coordinates: length-4-tuple with center coordinates
            x_left: x-coordinate of center of left frame half
            y_left: y-coordinate of center of left frame half
            x_right: x-coordinate of center of right frame half
            y_right: y-coordinate of center of right frame half
            Note that the center coordinates are defined in the complete frame, 
            i.e. with both detector halves (pixels; 0-based). The default value 
            is (477, 521, 1503, 511).         
        param_centering: length-3-tuple with parameters for centering by fitting 
            a 2D Gaussian or using cross-correlation:          
            crop_radius: half the length of the sides of the square cropped 
                sub-images used to fit the 2D Gaussian to and used for 
                cross-correlating the images (pixels). Must be integer. The 
                sub-image is centered on the coordinates as provided by
                center_coordinates. If None, the complete frame is used for the 
                fitting and center_coordinates is ignored. The value of
                crop_radius is also used to create the sub-images when
                show_images_center_coordinates = True.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of param_centering is (60, None, 30000). param_centering
            is only used when centering_method is 'gaussian' or 'cross-correlation'.
        collapse_ndit: If True, compute the mean over the (NDIT) frames of a
            file before subtracting the background, flat-fielding, bad pixel
            removal and centering. If False, perform the above steps for each
            frame and after that compute the mean over the frames 
            (default = False).
        show_images_center_coordinates: If True, plot the sub-images showing the 
            center coordinates for each frame. The plots allow for checking
            whether the centering is correct and to scan the data for frames 
            with bad quality (default = True). 
        
    Output:
        cube_left_frames: cube of pre-processed left frames 
        cube_right_frames: ccube of pre-processed right frames 
        header: list of FITS-headers of raw science frames   
                
    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: verified
    '''
    
    # Print centering method selected
    center_coordinates_print = tuple(x + 1 for x in center_coordinates)
    if centering_method == 'center frames':
        printandlog('\nCentering frames with center coordinates found from CENTER-file(s):')
        printandlog('(%.2f, %.2f, %.2f, %.2f)' % center_coordinates_print)
    elif centering_method == 'gaussian':
        printandlog('\nCentering frames by fitting a 2D Gaussian.')
        printandlog('center_coordinates = ' + str(center_coordinates_print))
        printandlog('param_centering = ' + str(param_centering)) 
    elif centering_method == 'cross-correlation':
        printandlog('\nCentering frames using cross-correlation.')
        printandlog('center_coordinates = ' + str(center_coordinates_print))
        printandlog('param_centering = ' + str(param_centering)) 
    elif centering_method == 'manual':
        printandlog('\nCentering frames with user-provided center coordinates:')
        printandlog('(%.2f, %.2f, %.2f, %.2f)' % center_coordinates_print)

    # Assemble center coordinates in two arrays and subtract 1024 from the right x-coordinate to make it valid for a frame half
    x_center_0 = np.array([center_coordinates[0], center_coordinates[2] - 1024])
    y_center_0 = np.array([center_coordinates[1], center_coordinates[3]])
 
    # Read centering parameters
    crop_radius, sigfactor, saturation_level = param_centering

    # Create zero arrays and empty lists
    list_frame_index = []
    list_file_index = []
    list_shift_x = [[], []]
    list_shift_y = [[], []]
    list_left_frames = []
    list_right_frames = []
    header = []

    if centering_method in ['gaussian', 'cross-correlation']:
        # Determine filter used and compute theoretical full width half maximum
        filter_used = pyfits.getheader(path_object_files[0])['ESO INS1 FILT ID']
        fwhm = compute_fwhm_separation(filter_used)[0]
    
    if show_images_center_coordinates == True:
        # Create empty lists
        list_x_fit_sub_image = [[], []]
        list_y_fit_sub_image = [[], []]
        list_sub_image = [[], []]
                        
    if centering_method == 'cross-correlation':
        # Create zero arrays
        x_fit_template = np.zeros(2)
        y_fit_template = np.zeros(2)
        x_fit_sub_image_template = np.zeros(2)
        y_fit_sub_image_template = np.zeros(2)
        list_sub_image_template = []
        
    printandlog('')
    
    for i, (path_sel, indices_sel) in enumerate(zip(path_object_files, indices_to_remove_object)):
        # Read data and header from file
        cube_sel, header_sel = read_fits_files(path=path_sel, silent=True)
        
        # Create list of indices of frames
        frame_index = [x for x in range(0, cube_sel.shape[0])]
          
        # Remove frames and frame indices based on list of indices provided by the user
        cube_sel = np.delete(cube_sel, indices_sel, axis=0)
        frame_index = [x for x in frame_index if x not in indices_sel]  

        if collapse_ndit == True:
            # Compute mean over NDIT frames
            cube_sel = np.mean(cube_sel, axis=0, keepdims=True)  
            
            if len(frame_index) > 1:
                # Change frame index to say 'collapsed' so that this can be printed when plotting sub-images
                frame_index = ['collapsed']
            
        # Subtract background and divide by master flat
        cube_bgsubtr_flatfielded = (cube_sel - frame_master_sky) / frame_master_flat

        # Remove bad pixels of each frame
        cube_badpixel_filtered = remove_bad_pixels(cube=cube_bgsubtr_flatfielded, frame_master_bpm=frame_master_bpm, sigma_filtering=sigma_filtering)

        # Create zero arrays to save centered images in
        cube_left_centered = np.zeros((cube_badpixel_filtered.shape[-3], cube_badpixel_filtered.shape[-2], int(cube_badpixel_filtered.shape[-1] / 2)))
        cube_right_centered = np.copy(cube_left_centered)

        # Retrieve dithering shifts in x- and y-direction from header
        x_dith = int(header_sel['ESO INS1 DITH POSX'])
        y_dith = int(header_sel['ESO INS1 DITH POSY'])

        # For each frame in the cube
        for j, frame_sel in enumerate(cube_badpixel_filtered):        
            # Separate left and right part of image
            frame_left = frame_sel[:, :1024]
            frame_right = frame_sel[:, 1024:]

            # Save file number and frame number of cube 
            list_file_index.append(file_index_object[i])
            list_frame_index.append(frame_index[j])
            
            # Assembly the frame halves in a list
            frame_halves = [frame_left, frame_right]           

            # For each frame half
            for k, (frame_half, x_center_0_sel, y_center_0_sel) in enumerate(zip(frame_halves, x_center_0, y_center_0)):
                if centering_method in ['center frames', 'manual']:           
                    # Compute x- and y-shifts for left and right images
                    list_shift_x[k].append(511.5 - x_center_0_sel - x_dith)
                    list_shift_y[k].append(511.5 - y_center_0_sel - y_dith)

                    if show_images_center_coordinates == True:
                        # Create sub-image to show center coordinates in
                        crop_radius = 12
                        x_center_0_rounded = int(np.round(x_center_0_sel))
                        y_center_0_rounded = int(np.round(y_center_0_sel))
                        sub_image = create_sub_image(frame=frame_half, x0=x_center_0_rounded + x_dith, y0=y_center_0_rounded + y_dith, crop_radius=crop_radius)
    
                        # Compute center position in coordinates of sub-image 
                        x_fit_sub_image = x_center_0_sel - x_center_0_rounded + crop_radius
                        y_fit_sub_image = y_center_0_sel - y_center_0_rounded + crop_radius

                elif centering_method == 'gaussian':
                    # Determine accurate coordinates of star position by fitting a Gaussian           
                    x_fit, y_fit, x_fit_sub_image, y_fit_sub_image, sub_image = \
                    fit_2d_gaussian(frame=frame_half, x0=x_center_0_sel + x_dith, y0=y_center_0_sel + y_dith, x_stddev=fwhm, y_stddev=fwhm, \
                                    theta=0.0, crop_radius=crop_radius, sigfactor=sigfactor, saturation_level=saturation_level)
    
                    # Compute shift in x- and y-directions
                    list_shift_x[k].append(511.5 - x_fit)
                    list_shift_y[k].append(511.5 - y_fit)

                elif centering_method == 'cross-correlation':
                    if i == 0 and j == 0:
                        # Center frame halves of first frame by fitting a Gaussian
                        if k == 0:
                            printandlog('Creating templates for left and right images from first image frame.\n')
                        x_fit_template[k], y_fit_template[k], x_fit_sub_image_template[k], y_fit_sub_image_template[k] = \
                        fit_2d_gaussian(frame=frame_half, x0=x_center_0_sel + x_dith, y0=y_center_0_sel + y_dith, x_stddev=fwhm, y_stddev=fwhm, \
                                        theta=0.0, crop_radius=crop_radius, sigfactor=sigfactor, saturation_level=saturation_level)[:4]
                        
                        # Subtract dithering shifts from fitted template coordinates
                        x_fit_template[k] -= x_dith
                        y_fit_template[k] -= y_dith
                        
                        # Create template sub-image
                        list_sub_image_template.append(create_sub_image(frame=frame_half, x0=x_center_0_sel + x_dith, y0=y_center_0_sel + y_dith, crop_radius=crop_radius))
                
                    # Create sub-image to cross-correlate with template sub-image
                    sub_image = create_sub_image(frame=frame_half, x0=x_center_0_sel + x_dith, y0=y_center_0_sel + y_dith, crop_radius=crop_radius)
                    
                    # Determine required shift of image by cross-correlation with template
                    x_shift_fit, y_shift_fit = register_translation(list_sub_image_template[k], sub_image, upsample_factor=10)[0]
                       
                    # Compute shift in x- and y-directions
                    list_shift_x[k].append(511.5 - x_fit_template[k] + x_shift_fit - x_dith)
                    list_shift_y[k].append(511.5 - y_fit_template[k] + y_shift_fit - y_dith)   

                    if show_images_center_coordinates == True:
                        # Compute fit position in coordinates of sub-image                   
                        x_fit_sub_image = x_fit_sub_image_template[k] - x_shift_fit
                        y_fit_sub_image = y_fit_sub_image_template[k] - y_shift_fit
                        
                if show_images_center_coordinates == True:
                    # Append sub-image and its fitted coordinates to lists
                    list_x_fit_sub_image[k].append(x_fit_sub_image)
                    list_y_fit_sub_image[k].append(y_fit_sub_image)
                    list_sub_image[k].append(sub_image)

            # Shift left and right images to center
            cube_left_centered[j, :, :] = np.expand_dims(ndimage.shift(frame_left, [list_shift_y[0][-1], list_shift_x[0][-1]], 
                                                         order=3, mode='constant', cval=0.0, prefilter=True), axis=0)   
            cube_right_centered[j, :, :] = np.expand_dims(ndimage.shift(frame_right, [list_shift_y[1][-1], list_shift_x[1][-1]], 
                                                         order=3, mode='constant', cval=0.0, prefilter=True), axis=0)  
       
        # Compute mean images of left and right image cubes
        frame_left_centered = np.mean(cube_left_centered, axis=0)
        frame_right_centered = np.mean(cube_right_centered, axis=0)

        # Append single sum and difference images and header
        list_left_frames.append(frame_left_centered)
        list_right_frames.append(frame_right_centered)
        header.append(header_sel)
        
        # Print which file has been processed
        printandlog('Processed file ' + str(i + 1) + '/' + str(len(path_object_files)) + ': {0:s}'.format(os.path.basename(path_sel)))

    # Convert lists of single sum and difference images to image cubes
    cube_left_frames = np.stack(list_left_frames)
    cube_right_frames = np.stack(list_right_frames)    
        
    # Determine type of observations and assign directory to save plots in
    observation_type = header_sel['ESO DPR TYPE']
    if observation_type == 'OBJECT':
        path_plots_dir = path_preprocessed_dir
    elif observation_type == 'OBJECT,FLUX':
         path_plots_dir = path_flux_dir
       
    if centering_method in ['gaussian', 'cross-correlation']:
        # Convert lists of shifts to arrays of center coordinates
        x_center = 511.5 - np.array(list_shift_x).T
        y_center = 511.5 - np.array(list_shift_y).T
       
        # Add 1024 to the x-coordinates of the right frame half to make the values valid for complete frame
        x_center[:, 1] += 1024 

        # Compute mean of determined center coordinates    
        x_center_mean = np.mean(x_center, axis=0)
        y_center_mean = np.mean(y_center, axis=0)

        if len(x_center[:, 0]) > 1:
            # Define function to plot determined center coordinates
            def plot_center_coordinates(data, x_y, left_right):
                font_size = 10
                width_figure = 4.7 / 20 * len(data)
                if width_figure < 4.7:
                    width_figure = 4.7   
                elif width_figure > 30:
                    width_figure = 40                  
                plot_name = name_file_root + 'center_coordinates_' + x_y + '_' + left_right + '.png'            
                path_plot = os.path.join(path_plots_dir, plot_name)
                printandlog(path_plot, wrap=False)
                image_number = np.arange(1, len(data)+1)
                plt.figure(figsize = (width_figure, 3.0))
                plt.plot(image_number, data, '-ok')
                ax = plt.gca()
                ax.set_xlabel(r'Image', fontsize=font_size)
                ax.tick_params(axis = 'x', labelsize=font_size)
                ax.set_xlim([0, max(image_number) + 1])
                ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
                ax.set_ylabel(x_y + ' ' + left_right + ' (pixels)', fontsize=font_size)
                ax.tick_params(axis='y', labelsize=font_size)  
                ax.ticklabel_format(useOffset=False, axis='y')
                ax.grid()
                plt.tight_layout()
                plt.savefig(path_plot, dpi=300, bbox_inches='tight')
                plt.close()

            # Plot center coordinates in x- and y-direction of left and right frame halves
            printandlog('\nCreating plots showing the fitted center coordinates of each image:')
            plot_center_coordinates(data=x_center[:, 0] + 1, x_y='x', left_right='left')
            plot_center_coordinates(data=y_center[:, 0] + 1, x_y='y', left_right='left')
            plot_center_coordinates(data=x_center[:, 1] + 1, x_y='x', left_right='right')
            plot_center_coordinates(data=y_center[:, 1] + 1, x_y='y', left_right='right') 
        
            # Compute standard deviation of determined center coordinates
            x_center_std = np.std(x_center, ddof=1, axis=0)
            y_center_std = np.std(y_center, ddof=1, axis=0)
    
            # Print mean center coordinates with error
            separator = ' '*max([len('%.2f' % x) for x in np.append(x_center_std, y_center_std)])
            printandlog('\nMean center coordinates (pixels):')
            printandlog('x_left' + separator + '         y_left' + separator + '         x_right' + separator + '         y_right')
            printandlog('%.2f +/- %.2f    %.2f +/- %.2f    %.2f +/- %.2f    %.2f +/- %.2f' 
                        % (x_center_mean[0] + 1, x_center_std[0], y_center_mean[0] + 1, y_center_std[0], x_center_mean[1] + 1, x_center_std[1], y_center_mean[1] + 1, y_center_std[1]))    
            
            # Print warning if there is significant deviation among the center coordinates found
            if any(np.append(x_center_std, y_center_std) > 0.5):
                printandlog('\nWARNING, for at least one of the fitted center coordinates the standard deviation is larger than 0.5 pixels.')
        else:
            # Print mean center coordinates without error
            printandlog('\nCenter coordinates (pixels):')
            printandlog('x_left    y_left    x_right    y_right')
            printandlog('%.2f    %.2f    %.2f    %.2f' % (x_center_mean[0] + 1, y_center_mean[0] + 1, x_center_mean[1] + 1, y_center_mean[1] + 1))   

    if show_images_center_coordinates == True:
        # Compute number of figures, rows and figure size
        number_frames_max = 10
        height_frame = 3.05
        width_figure = 4.7
        number_frames = len(list_sub_image[0])
        number_figures = int(np.ceil(number_frames / number_frames_max))
        number_rows_nominal = int(np.ceil(number_frames / number_figures))
        number_rows = np.ones(number_figures, dtype=int)*number_rows_nominal
        number_rows[-1] = int(number_frames - (number_figures - 1)*number_rows_nominal)
        height_figure = number_rows*height_frame
       
        # Plot sub-images used to determine center coordinates
        cmap = mpl.cm.viridis
        cmap.set_bad(color='gray')
        
        printandlog('\nCreating plot(s) showing the sub-images of each frame and the center coordinates fitted:')
        for i in range(number_figures):
            if number_figures == 1:
                plot_name = name_file_root + 'centering_sub_images.png'            
            else:
                plot_name = name_file_root + 'centering_sub_images_' + str(i + 1) + '.png'            
            path_plot = os.path.join(path_plots_dir, plot_name)
            printandlog(path_plot, wrap=False)
            fig, axs = plt.subplots(nrows=number_rows[i], ncols=2, sharex=True, sharey=True, \
                                    subplot_kw={'xticks': [], 'yticks': []}, figsize=(width_figure, height_figure[i]))
            fig.subplots_adjust(top=0.7)
            if len(list_shift_x[0]) > 1:
                axs[0, 0].set_title('Left frame halves')
                axs[0, 1].set_title('Right frame halves')
                for j in range(number_rows[i]):
                    j2 = j + number_rows[i-1]*i
                    if j2 < number_frames:
                        if list_frame_index[j2] == 'collapsed':
                            axs[j, 0].set_ylabel(str(j2 + 1) + ': file ' + str(list_file_index[j2] + 1) + ' (collapsed)')
                        else:
                            axs[j, 0].set_ylabel(str(j2 + 1) + ': file ' + str(list_file_index[j2] + 1) + ' frame ' + str(list_frame_index[j2] + 1))
                        for k in range(2):
                            axs[j, k].imshow(list_sub_image[k][j2], cmap=cmap, origin='lower',interpolation='nearest')
                            axs[j, k].plot(list_x_fit_sub_image[k][j2], list_y_fit_sub_image[k][j2], 'rx', markersize=60/crop_radius)
            else:
                axs[0].set_title('Left frame halves')
                axs[1].set_title('Right frame halves')
                if list_frame_index[0] == 'collapsed':
                    axs[0].set_ylabel(str(1) + ': file ' + str(list_file_index[0] + 1) + ' (collapsed)')
                else:
                    axs[0].set_ylabel(str(1) + ': file ' + str(list_file_index[0] + 1) + ' frame ' + str(list_frame_index[0] + 1))
                for k in range(2):
                    axs[k].imshow(list_sub_image[k][0], cmap=cmap, origin='lower',interpolation='nearest')
                    axs[k].plot(list_x_fit_sub_image[k][0], list_y_fit_sub_image[k][0], 'rx', markersize=60/crop_radius)
            plt.savefig(path_plot, dpi=300, bbox_inches='tight')
            plt.close(fig)

    return cube_left_frames, cube_right_frames, header

###############################################################################
# compute_annulus_values
###############################################################################    
    
def compute_annulus_values(cube, param):
    '''
    Obtain values of an image cube in an annulus (cube can have more than 
    3 dimensions; the cut is made in the last 2 dimensions)
    
    Input:
        cube: image cube or frame to obtain values from
        param: (list of) length-6-tuple(s) with parameters to generate annulus coordinates:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            
    Output:
        values_annulus: values of the image cube in the annulus
        frame_annulus: frame showing the annulus used to obtain the values
        
    File written by Rob van Holstein
    Function status: verified
    '''

    # Create meshgrid coordinates to construct annulus on
    x = np.arange(0, cube.shape[-1])
    y = np.arange(0, cube.shape[-2])
    xm, ym = np.meshgrid(x, y)

    # Create empty arrays for x- and y-coordinates
    coord_x_tot = np.array([])
    coord_y_tot = np.array([])
    
    # Make sure param is a list of tuples
    if type(param) == tuple:
        param = [param]

    for param_sel in param:
        coord_center_x, coord_center_y, inner_radius, outer_radius, start_angle, end_angle = param_sel

        start_angle = np.mod(start_angle, 360)
        end_angle = np.mod(end_angle, 360)
        
        # Of each pixel calculate radius and angle in range [0, 360)
        radius = np.sqrt((xm - coord_center_x)**2 + (ym - coord_center_y)**2)
        angle = np.mod(np.rad2deg(np.arctan2(ym - coord_center_y, xm - coord_center_x)), 360)
           
        # Select pixels that satisfy provided requirements
        if start_angle < end_angle:
            coord_y, coord_x = np.nonzero(np.logical_and(np.logical_and(radius >= inner_radius, radius < outer_radius), 
                                                         np.logical_and(angle >= start_angle, angle < end_angle)))
        else:
            coord_y1, coord_x1 = np.nonzero(np.logical_and(np.logical_and(radius >= inner_radius, radius < outer_radius), 
                                                           np.logical_and(angle >= start_angle, angle < 360)))
            coord_y2, coord_x2 = np.nonzero(np.logical_and(np.logical_and(radius >= inner_radius, radius < outer_radius), 
                                                           np.logical_and(angle >= 0, angle < end_angle)))
            coord_y, coord_x = np.hstack([coord_y1, coord_y2]), np.hstack([coord_x1, coord_x2])
        
        # Append coordinates to final coordinate arrays
        coord_x_tot = np.append(coord_x_tot, coord_x).astype(np.int)
        coord_y_tot = np.append(coord_y_tot, coord_y).astype(np.int)
        
    # Determine values 
    values_annulus = cube[..., coord_y_tot, coord_x_tot]
 
    # Create map with annulus coordinates
    frame_annulus = np.zeros(cube.shape[-2:], dtype = np.float32)
    frame_annulus[coord_y_tot, coord_x_tot] = 1.0

    return values_annulus, frame_annulus

###############################################################################
# subtract_background
###############################################################################

def subtract_background(cube, annulus_background):   
    '''
    Subtract background from cube or frame
     
    Input:
        cube: image cube or frame to subtract background from
        annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)

    Output:
        cube_background_subtracted: image cube or frame with background subtracted
        background: float or array of background values for each image
   
    File written by Rob van Holstein
    Function status: verified   
    '''

    # Determine background in frame or cube and subtract it
    if cube.ndim == 2:
        background = np.median(compute_annulus_values(cube=cube, param=annulus_background)[0])
        cube_background_subtracted = cube - background
    elif cube.ndim == 3:
        background = np.median(compute_annulus_values(cube=cube, param=annulus_background)[0], axis=1)
        cube_background_subtracted = cube - background[:, np.newaxis, np.newaxis] 
   
    return cube_background_subtracted, background

###############################################################################
# process_flux_frames
###############################################################################

def process_flux_frames(path_flux_files, 
                        file_index_flux, 
                        indices_to_remove_flux, 
                        frame_master_flat, frame_master_bpm, 
                        frame_master_sky_flux, 
                        annulus_background, 
                        sigma_filtering=True, 
                        centering_method='gaussian', 
                        center_coordinates=(477, 521, 1503, 511), 
                        param_centering=(60, None, 30000), 
                        collapse_ndit=False, 
                        show_images_center_coordinates=True):
    '''
    Create a master flux-frame from the FLUX-files. Function performs the same
    steps as process_object_frames, but only uses the single-sum (total 
    intensity) image and additionally subtracts the background
    
    Input:
        path_flux_files: list of paths to raw FLUX-files
        file_index_flux: list of file indices of FLUX-files (0-based)
        indices_to_remove_flux: list of arrays with indices of frames to remove
            for each FLUX-file        
        frame_master_flat: master flat frame
        frame_master_bpm: frame indicating location of bad pixels with 0's and 
            good pixels with 1's
        frame_master_sky_flux: master sky frame for FLUX-files
        annulus_background: (list of) length-6-tuple(s) with parameters 
            to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right
                and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and
                positive rotation counterclockwise)
        sigma_filtering: if True remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True)
        centering_method: method to center the images. If 'manual', use fixed 
            coordinates as provided by center_coordinates. If 'gaussian', fit 
            a 2D Gaussian to each frame. For 'gaussian' center_coordinates is 
            used as initial guess of the center coordinates and the determined 
            center coordinates are plotted for each image (default = 'gaussian').
        center_coordinates: length-4-tuple with center coordinates
            x_left: x-coordinate of center of left frame half
            y_left: y-coordinate of center of left frame half
            x_right: x-coordinate of center of right frame half
            y_right: y-coordinate of center of right frame half
            Note that the center coordinates are defined in the complete frame, 
            i.e. with both detector halves (pixels; 0-based). The default value 
            is (477, 521, 1503, 511).         
        param_centering: length-3-tuple with parameters for centering by fitting
            a 2D Gaussian:
            crop_radius: half the length of the sides of the square cropped 
                sub-images used to fit the 2D Gaussian to (pixels). Must be 
                integer. The sub-image is centered on the coordinates as 
                provided by center_coordinates. If None, the complete frame is 
                used for the fitting and center_coordinates is ignored. The
                value of crop_radius is also used to create the sub-images when
                show_images_center_coordinates = True.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of param_centering is (60, None, 30000). param_centering
            is only used when centering_method is 'gaussian'.    
        collapse_ndit: If True, compute the mean over the (NDIT) frames of a
            file before subtracting the background, flat-fielding, bad pixel
            removal and centering. If False, perform the above steps for each
            frame and after that compute the mean over the frames 
            (default = False).
        show_images_center_coordinates: If True, plot the sub-images showing the 
            center coordinates for each frame. The plots allow for checking
            whether the centering is correct and to scan the data for frames 
            with bad quality (default = True). 
        
    Output:
        frame_master_flux: master flux frame 
        frame_annulus_background: frame showing annulus used to determine background
            
    File written by Rob van Holstein
    Function status: verified
    '''

    # Perform dark-subtraction, flat-fielding, bad pixel removal and centering
    cube_left_frames, cube_right_frames = process_object_frames(path_object_files=path_flux_files, 
                                                                file_index_object=file_index_flux, 
                                                                indices_to_remove_object=indices_to_remove_flux, 
                                                                frame_master_flat=frame_master_flat, 
                                                                frame_master_bpm=frame_master_bpm, 
                                                                frame_master_sky=frame_master_sky_flux, 
                                                                sigma_filtering=sigma_filtering, 
                                                                centering_method=centering_method, 
                                                                center_coordinates=center_coordinates, 
                                                                param_centering=param_centering, 
                                                                collapse_ndit=collapse_ndit, 
                                                                show_images_center_coordinates=show_images_center_coordinates)[:2]

    # Compute flux frame as mean of single sum cube
    cube_single_sum = cube_left_frames + cube_right_frames
    frame_flux = np.mean(cube_single_sum, axis=0)

    # Determine background and subtract it
    frame_master_flux, background = subtract_background(cube=frame_flux, 
                                                        annulus_background=annulus_background)

    printandlog('\nSubtracted background in master flux image = %.3f' % background)
      
    # Create frame showing annulus used to determine background
    frame_annulus_background = compute_annulus_values(cube=frame_flux, param=annulus_background)[1]
    
    # Print number of frames used to create master flux frame
    number_frames_total = sum([pyfits.getheader(x)['ESO DET NDIT'] for x in path_flux_files])
    number_frames_removed = sum([len(x) for x in indices_to_remove_flux])
    number_frames_used = number_frames_total - number_frames_removed
    
    printandlog('\nCreated master flux frame out of ' + str(len(path_flux_files)) + 
            ' raw FLUX-file(s) comprising a total of ' + str(number_frames_used) + ' frame(s).')

    return frame_master_flux, frame_annulus_background

###############################################################################
# annulus_1_to_0_based
###############################################################################

def annulus_1_to_0_based(annulus):
    '''
    Convert annulus parameters from 1- to 0-based indexing
    
    Input:
        annulus: (list of) length-6-tuple(s) with parameters to generate annulus:
            coord_center_x: x-coordinate of center (pixels; 1-based)
            coord_center_y: y-coordinate of center (pixels; 1-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            
    Output:
        annulus: same as for input, but with coord_center_x and coord_center_y
            0-based
    
    File written by Rob van Holstein
    Function status: verified
    ''' 

    if type(annulus) is tuple:
        annulus = (annulus[0] - 1,) + (annulus[1] - 1,) + annulus[2:]
    elif type(annulus) is list:
        for i,x in enumerate(annulus):
            x = (x[0] - 1,) + (x[1] - 1,) + x[2:]
            annulus[i] = x
    
    return annulus        
            
###############################################################################
# annulus_0_to_1_based
###############################################################################

def annulus_0_to_1_based(annulus):
    '''
    Convert annulus parameters from 0- to 1-based indexing
    
    Input:
        annulus: (list of) length-6-tuple(s) with parameters to generate annulus:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            
    Output:
        annulus: same as for input, but with coord_center_x and coord_center_y
            1-based
    
    File written by Rob van Holstein
    Function status: verified
    ''' 

    if type(annulus) is tuple:
        annulus = (annulus[0] + 1,) + (annulus[1] + 1,) + annulus[2:]
    elif type(annulus) is list:
        for i,x in enumerate(annulus):
            x = (x[0] + 1,) + (x[1] + 1,) + x[2:]
            annulus[i] = x
    
    return annulus

###############################################################################
# perform_preprocessing
###############################################################################

def perform_preprocessing(frames_to_remove=[], 
                          sigma_filtering=True, 
                          object_collapse_ndit=False, 
                          show_images_center_coordinates=True, 
                          object_centering_method='automatic', 
                          center_subtract_object=True, 
                          object_center_coordinates='automatic', 
                          center_param_centering=(12, None, 30000),
                          object_param_centering=(60, None, 30000), 
                          flux_centering_method='gaussian', 
                          flux_center_coordinates=(477, 521, 1503, 511), 
                          flux_param_centering=(60, None, 30000), 
                          flux_annulus_background='large annulus',
                          save_preprocessed_data=True, 
                          combination_method_polarization='least squares',\
                          flux_annulus_star='automatic',\
                          convert_in_contrast_per_arcsec2=False):
    '''
    Perform pre-processing of OBJECT, CENTER, SKY and FLUX-files, i.e. sorting data, 
    background subtraction, flat-fielding, bad pixel removal, centering and compution
    of the single-sum and -difference images used in the post-processing

    Input:
        frames_to_remove: list of integers and length-2-tuples of integers 
            indicating which files and frames to remove (1-based). A complete
            file can be removed by specifying its integer index, while a frame
            of specific file can be removed by specifying a tuple
            (file_index, frame_index). If no files or frames should be removed, 
            use an empty list [] (default = []). The files are sorted in
            chronological order from oldest to newest.
        sigma_filtering: if True, remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True). Applies
            to all file-types (OBJECT, CENTER, SKY and FLUX).
        object_collapse_ndit: If True, compute the mean over the (NDIT) frames of
            the OBJECT-files before subtracting the background, flat-fielding, bad 
            pixel removal and centering to speed up the preprocessing. If False, 
            perform the above steps for each frame and after that compute the 
            mean over the frames (default = False).
        show_images_center_coordinates: If True, plot the sub-images showing the 
            center coordinates for each frame of the OBJECT- and FLUX-files. 
            The plots allow for checking whether the centering is correct and 
            to scan the data for frames with bad quality (default = True).    
        object_centering_method: method to center the OBJECT-frames. If 
            'center frames' determine the center coordinates from the 
            CENTER-frames. If 'gaussian', fit a 2D Gaussian to each frame. 
            If 'cross-correlation', fit a 2D Gaussian to the first frame and then
            use cross-correlation to align (register) the other frames onto the 
            centered first  frame. For 'gaussian' and 'cross-correlation' 
            object_center_coordinates is used as initial guess of the center 
            coordinates and the determined center coordinates are plotted for
            each image. If 'manual', use fixed coordinates as provided by 
            object_center_coordinates. If 'automatic', object_centering_method is 
            set to 'center frames' if there are CENTER-files, and is set to 'gaussian' 
            if there are no CENTER-files (default = 'automatic').
        center_subtract_object: if True subtract the OBJECT-file(s) taken
            closest in time from the CENTER-file(s) (default = True). This generally
            results in a more accurate determination of the center coordinates
            as the background and any other celestial objects in the field of view
            are suppressed. When the difference between the image orientations of 
            the CENTER- and OBJECT-frames is large (i.e. difference in derotator 
            position angle for field-tracking and difference in parallactic angle 
            for pupil-tracking), the OBJECT-frame will be rotated around the initial 
            guess of the centers as defined by object_center_coordinates before 
            subtracting it from the CENTER-file. 
        object_center_coordinates: length-4-tuple with center coordinates of OBJECT-frames:
            x_left: x-coordinate of center of left frame half
            y_left: y-coordinate of center of left frame half
            x_right: x-coordinate of center of right frame half
            y_right: y-coordinate of center of right frame half
            Note that the center coordinates are defined in the complete frame, 
            i.e. with both detector halves (pixels; 0-based). If string 
            'automatic', parameter is set to (477, 534, 1503, 523) when the 
            coronagraph used is N_ALC_Ks and otherwise to (477, 521, 1503, 511)
            (default = 'automatic').
        center_param_centering: length-3-tuple with parameters for 
            2D Gaussian fitting the centers of the satellite spots in the CENTER-frames:
            crop_radius: half the length of side of square cropped sub-images used 
                to fit 2D Gaussian to (pixels). Must be integer. If None, the complete 
                frame is used for the fitting and center_coordinates is ignored.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of center_param_centering is (12, None, 30000).
            center_param_centering is only used when object_centering_method
            is 'center frames'.
        object_param_centering: length-3-tuple with parameters for centering of 
            OBJECT-frames by fitting a 2D Gaussian or using cross-correlation:
            crop_radius: half the length of the sides of the square cropped 
                sub-images used to fit the 2D Gaussian to and used for 
                cross-correlating the images (pixels). Must be integer. The 
                sub-image is centered on the coordinates as provided by
                center_coordinates. If None, the complete frame is used for the 
                fitting and center_coordinates is ignored. The value of
                crop_radius is also used to create the sub-images when
                show_images_center_coordinates = True.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of object_param_centering is (12, 7, 30000). 
            object_param_centering is only used when object_centering_method is 
            'gaussian' or 'cross-correlation'.
        flux_centering_method: method to center the FLUX-frames. If 'manual', use fixed 
            coordinates as provided by flux_center_coordinates. If 'gaussian', fit 
            a 2D Gaussian to each frame. For 'gaussian' flux_center_coordinates is 
            used as initial guess of the center coordinates and the determined 
            center coordinates are plotted for each image (default = 'gaussian').
        flux_center_coordinates: length-4-tuple with center coordinates of FLUX-frames:
            x_left: x-coordinate of center of left frame half
            y_left: y-coordinate of center of left frame half
            x_right: x-coordinate of center of right frame half
            y_right: y-coordinate of center of right frame half
            Note that the center coordinates are defined in the complete frame, 
            i.e. with both detector halves (pixels; 0-based). The default value 
            is (477, 521, 1503, 511).         
        flux_param_centering: length-3-tuple with parameters for centering of 
            FLUX-frames by fitting a 2D Gaussian:
            crop_radius: half the length of the sides of the square cropped 
                sub-images used to fit the 2D Gaussian to (pixels). Must be 
                integer. The sub-image is centered on the coordinates as 
                provided by center_coordinates. If None, the complete frame is 
                used for the fitting and center_coordinates is ignored. The
                value of crop_radius is also used to create the sub-images when
                show_images_center_coordinates = True.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of flux_param_centering is (60, None, 30000).
            flux_param_centering is only used when flux_centering_method is 'gaussian'.  
        flux_annulus_background: (list of) length-6-tuple(s) with parameters 
            to generate annulus to measure and subtract background in master flux frame:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right
                and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and
                positive rotation counterclockwise)
            If string 'large annulus' the annulus will be star-centered and 
            located far away from the star with an inner radius of 320 pixels
            and an out radius of 380 pixels (default = 'large annulus').
        save_preprocessed_data: If True, save preprocessed cubes of single-sum 
            and single-difference images in the 'preprocessed' folder so that 
            the preprocessing can be skipped when re-running the pipeline
            (default = True). 
        combination_method_polarization: method to be used to produce the 
            incident Q- and U-images, 'least squares', 'trimmed mean' or 'median' 
            (default = 'least squares')
        flux_annulus_star: (list of) length-6-tuple(s) with parameters 
            to generate an annulus (normally a full circular aperture but we let 
            the option to exclude an inner circle or a wedge in special circumstances like 
            a binary system) to measure the star total flux in master flux frame:
            coord_center_x: x-coordinate of center (pixels; 0-based; by default 255.5)
            coord_center_y: y-coordinate of center (pixels; 0-based; by default 255.5)
            inner_radius: inner radius (pixels; by default 0)
            outer_radius: outer_radius (pixels; by default 140)
            start_angle: start angle of annulus sector (deg; 0 deg to the right
                and positive rotation counterclockwise; by default 0)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and
                positive rotation counterclockwise; by default 360)        
        convert_in_contrast_per_arcsec2: if True, performs aperture photometry 
            on the master flux frame, and convert all the output frames in contrast 
            per arcsec^2.
        
    Output:
        cube_left_frames: cube of pre-processed left frames 
        cube_right_frames: cube of pre-processed right frames 
        header: list of FITS-headers of raw science frames 
        file_index_object: list of file indices of OBJECT-files (0-based)
        combination_method_polarization: method to be used to produce the 
            incident Q- and U-images, 'least squares', 'trimmed mean' or 'median'
        star_flux: float with the total star flux as computed using the parameters 
            defined by the flux_annulus_star tuple, in ADU, corrected by the DIT ratio
            between FLUX and OBJECT frames and the potential use of ND filters
            
    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: verified
    '''

    ###############################################################################
    # Checking and sorting data and creating directories 
    ###############################################################################
    
    # Check and sort data and create directories
    printandlog('\n###############################################################################')
    printandlog('# Checking and sorting data and creating directories')
    printandlog('###############################################################################') 

    # Check and sort data, and create directories   
    path_dark_files, path_flat_files, path_object_files, path_sky_files, path_center_files, \
    path_object_center_files, path_flux_files, path_sky_flux_files, \
    indices_to_remove_dark, indices_to_remove_flat, indices_to_remove_object, \
    indices_to_remove_sky, indices_to_remove_center, indices_to_remove_object_center, \
    indices_to_remove_flux, indices_to_remove_sky_flux, file_index_object, \
    file_index_flux, object_centering_method, combination_method_polarization \
    = check_sort_data_create_directories(frames_to_remove=frames_to_remove, 
                                         combination_method_polarization=combination_method_polarization, 
                                         object_centering_method=object_centering_method, 
                                         save_preprocessed_data=save_preprocessed_data, 
                                         show_images_center_coordinates=show_images_center_coordinates)

    ###############################################################################
    # Computing master flat and bad pixel map or reading static ones
    ###############################################################################
    
    if any(path_flat_files):
        # Process the dark and flat files to create a master flat and bad pixel map
        printandlog('\n###############################################################################')
        printandlog('# Processing DARK(,BACKGROUND)- and FLAT-files')
        printandlog('###############################################################################') 
        frame_master_flat, frame_master_bpm = process_dark_flat_frames(path_dark_files=path_dark_files, 
                                                                       path_flat_files=path_flat_files,
                                                                       indices_to_remove_dark=indices_to_remove_dark, 
                                                                       indices_to_remove_flat=indices_to_remove_flat)

        printandlog('')
        write_fits_files(data=frame_master_flat, path=os.path.join(path_flat_dir, name_file_root + 'master_flat.fits'), header=False, silent=False)
        write_fits_files(data=frame_master_bpm, path=os.path.join(path_bpm_dir, name_file_root + 'master_badpix.fits'), header=False, silent=False)
        
    else:
        # Determine filter used
        printandlog('\n###############################################################################')
        printandlog('# Reading static master flat and bad pixel map')
        printandlog('###############################################################################') 
        filter_used = pyfits.getheader(path_object_files[0])['ESO INS1 FILT ID']
    
        # Read static master flat    
        if filter_used in ['FILT_BBF_Y', 'FILT_NBF_HeI']:
            path_static_flat = os.path.join(path_static_calib_dir, 'master_flat_Y.fits')
            filter_static_flat = 'FILT_BBF_Y'
        elif filter_used in ['FILT_BBF_J', 'FILT_NBF_ContJ', 'FILT_NBF_PaB']:
            path_static_flat = os.path.join(path_static_calib_dir, 'master_flat_J.fits')
            filter_static_flat = 'FILT_BBF_J'
        elif filter_used in ['FILT_BBF_H', 'FILT_NBF_ContH', 'FILT_NBF_FeII']:
            path_static_flat = os.path.join(path_static_calib_dir, 'master_flat_H.fits')
            filter_static_flat = 'FILT_BBF_H'
        elif filter_used in ['FILT_BBF_Ks', 'FILT_NBF_ContK1', 'FILT_NBF_H2', 'FILT_NBF_BrG', 'FILT_NBF_CntK2', 'FILT_NBF_CO']:
            path_static_flat = os.path.join(path_static_calib_dir, 'master_flat_Ks.fits')
            filter_static_flat = 'FILT_BBF_Ks'
            
        frame_master_flat = np.squeeze(read_fits_files(path=path_static_flat, silent=True)[0])
       
        # Read static bad pixel map
        frame_master_bpm = np.squeeze(read_fits_files(path=os.path.join(path_static_calib_dir, 'master_badpix.fits'), silent=True)[0])
        printandlog('\nRead static bad pixel map and static master flat in ' + filter_static_flat + '.')
        
    ###############################################################################
    # Computing master sky for object images
    ###############################################################################
    
    if any(path_sky_files):
        # Process the sky files for the object files
        printandlog('\n###############################################################################')
        printandlog('# Processing SKY-files for OBJECT-files')
        printandlog('###############################################################################') 

        frame_master_sky = process_sky_frames(path_sky_files=path_sky_files, 
                                              indices_to_remove_sky=indices_to_remove_sky)
        
        # Write master sky-frame
        printandlog('')
        write_fits_files(data=frame_master_sky, path=os.path.join(path_sky_dir, name_file_root + 'master_sky.fits'), header=False, silent=False)
    
    else:
        # Create a master sky frame with only zeros
        frame_master_sky = np.zeros((1024, 2048))

    ###############################################################################
    # Processing center files and extracting center coordinates
    ###############################################################################
    
    if object_centering_method == 'center frames':
        # Print that we process the center files
        printandlog('\n###############################################################################')
        printandlog('# Processing CENTER-files')
        printandlog('###############################################################################')       

        if object_center_coordinates == 'automatic':
            # Determine coronagraph used and set center coordinates
            coronagraph_used = pyfits.getheader(path_object_files[0])['ESO INS COMB ICOR']
    
            if coronagraph_used == 'N_ALC_Ks':
                object_center_coordinates = (477, 534, 1503, 523)
                printandlog('\nobject_center_coordinates is \'automatic\': setting it to ' + str(tuple(x + 1 for x in object_center_coordinates)) + ' because the coronagraph used is N_ALC_Ks.')
            else:
                object_center_coordinates = (477, 521, 1503, 511)
                printandlog('\nobject_center_coordinates is \'automatic\': setting it to ' + str(tuple(x + 1 for x in object_center_coordinates)) + ' because the coronagraph used is not N_ALC_Ks.')
                    
        # Process the center files            
        list_frame_center_processed, header_center = process_center_frames(path_center_files=path_center_files, 
                                                                           indices_to_remove_center=indices_to_remove_center, 
                                                                           path_object_center_files=path_object_center_files, 
                                                                           indices_to_remove_object_center=indices_to_remove_object_center, 
                                                                           frame_master_flat=frame_master_flat, 
                                                                           frame_master_bpm=frame_master_bpm, 
                                                                           frame_master_sky=frame_master_sky, 
                                                                           center_subtract_object=center_subtract_object, 
                                                                           center_coordinates=object_center_coordinates, 
                                                                           sigma_filtering=sigma_filtering)
    
        # Write processed center frames
        path_processed_center_files = [os.path.join(path_center_dir, os.path.splitext(os.path.basename(x))[0] + '_processed.fits') for x in path_center_files]
        printandlog('')
        write_fits_files(data=list_frame_center_processed, path=path_processed_center_files, header=header_center, silent=False)
        
        # Find center coordinates and replace values of center_coordinates
        object_center_coordinates = find_center_coordinates(list_frame_center_processed=list_frame_center_processed, 
                                                            path_processed_center_files=path_processed_center_files, 
                                                            center_coordinates=object_center_coordinates, 
                                                            param_centering=center_param_centering)

    ###############################################################################
    # Creating processed and centered single-sum and -difference images
    ###############################################################################
    
    # Create reduced and centerd single-sum and -difference images
    printandlog('\n###############################################################################')
    printandlog('# Processing OBJECT-files')
    printandlog('###############################################################################') 

    if object_center_coordinates == 'automatic':
        # Determine coronagraph used and set center coordinates
        coronagraph_used = pyfits.getheader(path_object_files[0])['ESO INS COMB ICOR']

        if coronagraph_used == 'N_ALC_Ks':
            object_center_coordinates = (477, 534, 1503, 523)
            printandlog('\nobject_center_coordinates is \'automatic\': setting it to ' + str(tuple(x + 1 for x in object_center_coordinates)) + ' because the coronagraph used is N_ALC_Ks.')
        else:
            object_center_coordinates = (477, 521, 1503, 511)
            printandlog('\nobject_center_coordinates is \'automatic\': setting it to ' + str(tuple(x + 1 for x in object_center_coordinates)) + ' because the coronagraph used is not N_ALC_Ks.')
        
    cube_left_frames, cube_right_frames, header = process_object_frames(path_object_files=path_object_files, 
                                                                        file_index_object=file_index_object, 
                                                                        indices_to_remove_object=indices_to_remove_object, 
                                                                        frame_master_flat=frame_master_flat, 
                                                                        frame_master_bpm=frame_master_bpm, 
                                                                        frame_master_sky=frame_master_sky, 
                                                                        sigma_filtering=sigma_filtering, 
                                                                        centering_method=object_centering_method, 
                                                                        center_coordinates=object_center_coordinates, 
                                                                        param_centering=object_param_centering, 
                                                                        collapse_ndit=object_collapse_ndit, 
                                                                        show_images_center_coordinates=show_images_center_coordinates)
    
    if save_preprocessed_data == True:
        # Write preprocessed cubes of single-sum and single-difference images 
        printandlog('\nSaving pre-processed data so that pre-processing can be skipped the next time.')
        printandlog('')
        write_fits_files(data=cube_left_frames, path=os.path.join(path_preprocessed_dir, 'cube_left_frames.fits'), header=False, silent=False)
        write_fits_files(data=cube_right_frames, path=os.path.join(path_preprocessed_dir, 'cube_right_frames.fits'), header=False, silent=False)

        # Write path of object files to a .txt-file to be able to read headers
        with open(os.path.join(path_preprocessed_dir, 'path_object_files.txt'), 'w') as fh:
            for path_sel in path_object_files:
                fh.write('%s\n' % path_sel)
        printandlog('Wrote file ' + os.path.join(path_preprocessed_dir, 'path_object_files.txt') + '.', wrap=False)

        # Write indices of OBJECT-files to a .txt-file to be able to read them later for plot of header angles
        with open(os.path.join(path_preprocessed_dir, 'file_index_object.txt'), 'w') as fh:
            fh.write('%s' % file_index_object)
        printandlog('Wrote file ' + os.path.join(path_preprocessed_dir, 'file_index_object.txt') + '.', wrap=False)

    ###############################################################################
    # Computing master sky for flux images
    ###############################################################################
    
    if any(path_sky_flux_files):
        # Process the sky files for the flux files
        printandlog('\n###############################################################################')
        printandlog('# Processing SKY-files for FLUX-files')
        printandlog('###############################################################################') 

        frame_master_sky_flux = process_sky_frames(path_sky_files=path_sky_flux_files, 
                                                   indices_to_remove_sky=indices_to_remove_sky_flux)
        
        # Write master sky-frame
        printandlog('')
        write_fits_files(data=frame_master_sky_flux, path=os.path.join(path_sky_flux_dir, name_file_root + 'master_sky_flux.fits'), header=False, silent=False)

    else:
        # Create a master sky frame with only zeros
        frame_master_sky_flux = np.zeros((1024, 2048))
    
    ###############################################################################
    # Computing master flux image
    ###############################################################################
            
    if any(path_flux_files):
        # Print that we process the flux files
        printandlog('\n###############################################################################')
        printandlog('# Processing FLUX-files')
        printandlog('###############################################################################') 

        # Define and print annulus to determine the background from
        if type(flux_annulus_background) == tuple or type(flux_annulus_background) == list:
            printandlog('\nThe background will be determined with a user-defined annulus or several user-defined annuli:')
            if type(flux_annulus_background) == tuple:
                printandlog(annulus_0_to_1_based(flux_annulus_background))
            elif type(flux_annulus_background) == list:
                for x in flux_annulus_background:
                    printandlog(annulus_0_to_1_based(x))
        elif flux_annulus_background == 'large annulus':
            flux_annulus_background = (511.5, 511.5, 320, 380, 0, 360)
            printandlog('\nThe background will be determined with a star-centered annulus located far away from the star:')
            printandlog(annulus_0_to_1_based(flux_annulus_background))

        # Processthe flux files
        frame_master_flux, frame_annulus_background_flux = process_flux_frames(path_flux_files=path_flux_files, 
                                                                               file_index_flux=file_index_flux, 
                                                                               indices_to_remove_flux=indices_to_remove_flux, 
                                                                               frame_master_flat=frame_master_flat, 
                                                                               frame_master_bpm=frame_master_bpm, 
                                                                               frame_master_sky_flux=frame_master_sky_flux, 
                                                                               annulus_background=flux_annulus_background, 
                                                                               sigma_filtering=sigma_filtering, 
                                                                               centering_method=flux_centering_method, 
                                                                               center_coordinates=flux_center_coordinates, 
                                                                               param_centering=flux_param_centering, 
                                                                               collapse_ndit=False, 
                                                                               show_images_center_coordinates=show_images_center_coordinates)

        # Write master flux-frame and frame showing annulus used to determine background
        printandlog('')
        write_fits_files(data=frame_master_flux, path=os.path.join(path_flux_dir, name_file_root + 'master_flux.fits'), header=False, silent=False)
        write_fits_files(data=frame_annulus_background_flux, path=os.path.join(path_flux_dir, name_file_root + 'annulus_background_flux.fits'), header=False)    

        if flux_annulus_star == 'automatic':
            flux_annulus_star = (511.5, 511.5, 0, 140, 0, 360)
            printandlog('\nflux_annulus_star is \'automatic\': setting it to ' + str(tuple(x + 1 for x in flux_annulus_star)))
        printandlog('\nThe star total flux will be determined from the master flux image with an aperture located on:')
        printandlog(annulus_0_to_1_based(flux_annulus_star))
                        
        star_flux, dit_ratio,nd_ratio = determine_star_flux(\
                                            frame_master_flux,path_flux_files, \
                                            path_object_files, flux_annulus_star)
        reference_flux = star_flux*nd_ratio*dit_ratio
        table_star_flux = pd.DataFrame({'aperture parameters':[flux_annulus_star],\
                                        'encircled_flux (ADU)':[star_flux],\
                                        'DIT ratio (OBJECT/FLUX)':[dit_ratio],\
                                        'ND filter ratio (OBJECT/FLUX)':[nd_ratio],\
                                        'reference flux (ADU)':[reference_flux]})
        table_star_flux.to_csv(os.path.join(path_flux_dir, name_file_root + 'reference_flux.csv'),\
                               index=False)

    return cube_left_frames, cube_right_frames, header, file_index_object, \
            combination_method_polarization,reference_flux

################################################################################
## determine_star_flux
################################################################################
def determine_star_flux(frame_master_flux,path_flux_files,path_object_files, flux_annulus_star):   
    '''
    Determine flux of star in the master flux frame using aperture photometry 
     
    Input:
        frame_master_flux: master flux frame 
        path_flux_files: list of paths to raw FLUX-files
        path_object_files: list of paths to raw OBJECT-files
        flux_annulus_star: (list of) length-6-tuple(s) with parameters 
            to generate an annulus (normally a full circular aperture but we let 
            the option to exclude an inner circle or a wedge in special circumstances like 
            a binary system) to measure the star total flux in master flux frame:
            coord_center_x: x-coordinate of center (pixels; 0-based; by default 255.5)
            coord_center_y: y-coordinate of center (pixels; 0-based; by default 255.5)
            inner_radius: inner radius (pixels; by default 0)
            outer_radius: outer_radius (pixels; by default 140)
            start_angle: start angle of annulus sector (deg; 0 deg to the right
                and positive rotation counterclockwise; by default 0)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and
                positive rotation counterclockwise; by default 360)        
    Output:
        star_total_flux: flux of the star in ADU encircled in the region defined
                        by the parameters flux_annulus_star.
                        The reference flux to use to convert from ADU to contrast 
                        is the product star_total_flux *  dit_ratio * nd_ratio
        dit_ratio: ratio between the DIT of OBJECT frames over FLUX frames
        nd_ratio: ratio between the transmission of the neutral density used for
                    OBJECT frames over FLUX frames
        corrected by the transmission
                    of the neutral density filter used with the FLUX frame.
        
    File written by Julien Milli
    Function status: 
    '''
    # We first perform aperture photometry on the master FLUX frame    
    star_total_flux = np.sum(compute_annulus_values(cube=frame_master_flux, param=flux_annulus_star)[0])
    printandlog('The star reference flux encircled in the annulus ({0:s}) is {1:.1f} ADU'.format(\
                ', '.join(['{0:3.1f}'.format(f) for f in flux_annulus_star]),star_total_flux))

    # Determine filter used
    filter_used = pyfits.getheader(path_object_files[0])['ESO INS1 FILT ID']

    # Determine the DIT of the OBJECT
    dit_object = pyfits.getheader(path_object_files[0])['ESO DET SEQ1 DIT']

    # Determine the DIT of the FLUX
    dit_flux = pyfits.getheader(path_flux_files[0])['ESO DET SEQ1 DIT']

    # Determine the ND filter used for the OBJECT files
    nd_object_used = pyfits.getheader(path_object_files[0])['ESO INS4 FILT2 NAME']

    # Determine the ND filter used for the FLUX files
    nd_flux_used = pyfits.getheader(path_flux_files[0])['ESO INS4 FILT2 NAME']

    if 'K' in filter_used:
        BB_filter='B_Ks'
    elif 'H' in filter_used:
        BB_filter='B_H'
    elif 'Y' in filter_used:
        BB_filter='B_Y'
    elif 'J' in filter_used:
        BB_filter='B_J'
    
    # utility function to retriev the neutral density number 0.0, 1.0, 2.0 or 3.5
    def get_ND_value(ND_filter_name):
        if 'OPEN' in ND_filter_name:
            return 0
        else:
            return float(ND_filter_name[-3:])
    
    tr_object = sphere_transmission(BB_filter=BB_filter, NDset=get_ND_value(nd_object_used))
    tr_flux = sphere_transmission(BB_filter=BB_filter, NDset=get_ND_value(nd_flux_used))

    nd_ratio = tr_object/tr_flux
    dit_ratio = dit_object/dit_flux  

    return star_total_flux, dit_ratio,nd_ratio     

def sphere_transmission(BB_filter='B_H', NDset=0.):
    '''
    
    Input:
        - BB_filter: name of the broad band filter (among 'B_H',
                    'B_Hcmp2', 'B_J', 'B_Ks', 'B_ND-H', 'B_Y'). By default, assumes 'B_H'
        - NDset: ND filter (float) among 0., 1., 2. or 3.5. By default 0.
    Output:
        - transmission of the neutral density filter (between 0 and 1)
    File written by Julien Milli based on transmission curves by Arthur Vigan
    Function status: 
    '''
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
    
    return t_final_nd / t_final

###############################################################################
# compute_double_sum_double_difference
###############################################################################

def compute_double_sum_double_difference(cube_single_sum, cube_single_difference, header, double_difference_type='conventional'):
    '''  
    Compute double-sum I_Q- and I_U-images and double-difference Q- and U-images
    
    Input:
        cube_single_sum: cube of single-sum I_Q^+, I_Q^-, I_U^+ and I_U^- intensity images
        cube_single_difference: cube of single-difference Q^+, Q^-, U^+ and U^- images
        header: list of FITS-headers of raw science frames    
        double_difference_type: type of double difference to be computed, either
        'conventional' or 'normalized' (see van Holstein et al. 2019; default = 'conventional')
               
    Output:
        cube_I_Q_double_sum: cube of double-sum I_Q-images
        cube_I_U_double_sum: cube of double-sum I_U-images
        cube_Q_double_difference: cube of double-difference Q-images
        cube_U_double_difference: cube of double-difference U-images

    File written by Rob van Holstein
    Function status: verified
    ''' 

    # Filter for zeros and NaN's
    cube_single_sum = np.nan_to_num(cube_single_sum)
    cube_single_sum[cube_single_sum == 0] = 1
    cube_single_difference = np.nan_to_num(cube_single_difference)
    cube_single_difference[cube_single_difference == 0] = 1
    
    # Determine indices of files with Qplus, Qminus, Uplus and Uminus   
    stokes_parameter = np.array([x['ESO OCS DPI H2RT STOKES'] for x in header])              
    indices_Qplus = np.nonzero(stokes_parameter == 'Qplus')[0]
    indices_Qminus = np.nonzero(stokes_parameter == 'Qminus')[0]
    indices_Uplus = np.nonzero(stokes_parameter == 'Uplus')[0]
    indices_Uminus = np.nonzero(stokes_parameter == 'Uminus')[0] 
    
    # Compute double sum I_Q- and I_U-cubes
    cube_I_Q_double_sum = 0.5*(cube_single_sum[indices_Qplus, :, :] + cube_single_sum[indices_Qminus, :, :])
    cube_I_U_double_sum = 0.5*(cube_single_sum[indices_Uplus, :, :] + cube_single_sum[indices_Uminus, :, :])
       
    if double_difference_type == 'conventional':
        # Compute Q- and U-cubes using standard double difference
        printandlog('\nUsing the standard double difference to compute the Q- and U-images.')
        cube_Q_double_difference = 0.5*(cube_single_difference[indices_Qplus, :, :] - cube_single_difference[indices_Qminus, :, :])
        cube_U_double_difference = 0.5*(cube_single_difference[indices_Uplus, :, :] - cube_single_difference[indices_Uminus, :, :])
    
    elif double_difference_type == 'normalized':
        # Compute Q- and U-cubes using normalized double difference
        printandlog('\nUsing the normalized double difference to compute the Q- and U-images')   
        cube_qplus = cube_single_difference[indices_Qplus, :, :] / cube_single_sum[indices_Qplus, :, :]
        cube_qminus = cube_single_difference[indices_Qminus, :, :] / cube_single_sum[indices_Qminus, :, :]
        cube_uplus = cube_single_difference[indices_Uplus, :, :] / cube_single_sum[indices_Uplus, :, :]
        cube_uminus = cube_single_difference[indices_Uminus, :, :] / cube_single_sum[indices_Uminus, :, :]
    
        cube_q = np.nan_to_num(0.5*(cube_qplus - cube_qminus))
        cube_u = np.nan_to_num(0.5*(cube_uplus - cube_uminus))
    
        cube_Q_double_difference = cube_q*cube_I_Q_double_sum
        cube_U_double_difference = cube_u*cube_I_U_double_sum
        
    return cube_I_Q_double_sum, cube_I_U_double_sum, cube_Q_double_difference, cube_U_double_difference
    
###############################################################################
# remove_detector_artefact
###############################################################################

def remove_detector_artefact(cube, number_pixels):
    '''
    Remove vertical band detector artefact seen in Stokes Q- and U-images
    
    Input:
        cube: image cube or frame to be corrected
        number_pixels: number of pixels on the top and bottom of each pixel column
                       used to compute the median to be subtracted
    
    Output:
        cube_artefact_removed: image cube or frame with the artefact corrected
        
    File written by Rob van Holstein
    Function status: verified
    '''
    
    # If the input is a frame turn it into a cube
    cube_ndim = cube.ndim
    if cube_ndim == 2:
        cube = np.expand_dims(cube, axis=0)
      
    # Remove the vertical band artefact
    indices = np.r_[0:number_pixels, (cube.shape[-2] - number_pixels):cube.shape[-2]]
    median_columns = np.median(cube[:, indices, :], axis = 1)
    cube_artefact_removed = cube - median_columns[:, np.newaxis, :]

    # If the input is a frame turn the resulting cube back into a frame
    if cube_ndim == 2:
        cube_artefact_removed = np.squeeze(cube_artefact_removed)

    return cube_artefact_removed

###############################################################################
# determine_star_polarization
###############################################################################

def determine_star_polarization(cube_I_Q, cube_I_U, cube_Q, cube_U, annulus_star, annulus_background):   
    '''
    Determine polarization of star in annulus
     
    Input:
        cube_I_Q: cube of I_Q-images
        cube_I_U: cube of I_U-images
        cube_Q: cube of Q-images
        cube_U: cube of U-images
        annulus_star: (list of) length-6-tuple(s) with parameters to generate annulus to measure polarization of star:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
        annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
              
    Output:
        q: normalized Stokes q measured in annulus
        u: normalized Stokes u measured in annulus
   
    File written by Rob van Holstein
    Function status: verified   
    '''

    # Set axis parameter of mean and median depending on whether the cube is a 3D or 2D object
    if cube_I_Q.ndim == 2:
        mean_median_axis = None
    elif cube_I_Q.ndim == 3:
        mean_median_axis = 1
           
    # Compute flux in I_Q, I_U, Q and U in an annulus minus the background in an annulus
    I_Q = np.mean(compute_annulus_values(cube=cube_I_Q, param=annulus_star)[0], axis=mean_median_axis) - \
                  np.median(compute_annulus_values(cube=cube_I_Q, param=annulus_background)[0], axis=mean_median_axis)       
    I_U = np.mean(compute_annulus_values(cube=cube_I_U, param=annulus_star)[0], axis=mean_median_axis) - \
                  np.median(compute_annulus_values(cube=cube_I_U, param=annulus_background)[0], axis=mean_median_axis)      
    Q = np.mean(compute_annulus_values(cube=cube_Q, param=annulus_star)[0], axis=mean_median_axis) - \
                np.median(compute_annulus_values(cube=cube_Q, param=annulus_background)[0], axis=mean_median_axis)       
    U = np.mean(compute_annulus_values(cube=cube_U, param=annulus_star)[0], axis=mean_median_axis) - \
                np.median(compute_annulus_values(cube=cube_U, param=annulus_background)[0], axis=mean_median_axis)    
    
    # Compute normalized Stokes q and u
    q = Q / I_Q
    u = U / I_U
    
    return q, u

###############################################################################
# compute_mean_angle
###############################################################################
  
def compute_mean_angle(angles, degree_radian='degree', axis=None):
    '''
    Calculate mean of angles using mean of circular quantities
    
    Input:
        angles: list or array of angles
        degree_radian: if 'degree' input and output are in degree; else in radians
        axis: axis to compute mean over; None or int or tuple of ints (default = None)
    
    Output:
        mean_angle: scalar of mean angle

    File written by Rob van Holstein
    Function status: verified  
    '''
    
    # Convert angle to rad if specified in degree
    if degree_radian == 'degree':
        angles = np.deg2rad(angles)
    
    # Compute mean angle
    y = np.mean(np.sin(angles), axis = axis)
    x = np.mean(np.cos(angles), axis = axis)
    mean_angle = np.arctan2(y, x)

    # Convert mean angle to degree if input angle specified in degree
    if degree_radian == 'degree':
        mean_angle = np.rad2deg(mean_angle)
    
    return mean_angle

###############################################################################
# compute_rotation_mueller_matrix
###############################################################################

def compute_rotation_mueller_matrix(alpha):
    '''
    Calculate the Mueller matrix describing the rotation of the reference frame 

    Input:
        alpha: angle of rotation (deg)

    Output: 
        T: rotation Mueller matrix
    
    File written by Rob van Holstein
    Function status: verified
    ''' 

    # Convert angle to radians
    alpha = np.deg2rad(alpha)
    
    # Construct matrix
    T = np.zeros((4,4))
    T[0,0] = T[3,3] = 1
    T[1,1] = T[2,2] = np.cos(2*alpha)
    T[1,2] = np.sin(2*alpha)
    T[2,1] = -1*T[1,2] 
    
    return T
           
###############################################################################
# compute_reflection_mueller_matrix
###############################################################################

def compute_reflection_mueller_matrix(epsilon, Delta):
    '''
    Calculate the Mueller matrix describing a linear diattenuator and retarder 
    (e.g. reflection off a mirror)

    Input:
        epsilon: linear diattenuation (-1 <= epsilon <= 1, ideally epsilon = 0)
        Delta: linear retardance (deg; 0 deg <= Delta < 360 deg, ideally Delta = 180 deg)

    Output:     
        M: Mueller matrix of linear diattenuator and retarder
    
    File written by Rob van Holstein
    Function status: verified
    ''' 

    # Convert angle to radians   
    Delta = np.deg2rad(Delta)

    # Construct matrix 
    M = np.zeros((4,4))
    M[0,0] = M[1,1] = 1
    M[0,1] = M[1,0] = epsilon
    M[2,2] = M[3,3] = np.sqrt(1 - epsilon**2)*np.cos(Delta)
    M[2,3] = np.sqrt(1 - epsilon**2)*np.sin(Delta)
    M[3,2] = -1*M[2,3]
    
    return M

###############################################################################
# compute_polarizer_mueller_matrix
###############################################################################

def compute_polarizer_mueller_matrix(d):
    '''
    Calculate the Mueller matrix of a polarizer

    Input:
        d: diattenuation of polarizer (-1 <= d <= 1)

    Output:      
        M: Mueller matrix of a polarizer
        
    File written by Rob van Holstein
    Function status: verified
    '''  

    # Construct matrix     
    M = np.zeros((4,4))
    M[0,0] = M[1,1] = 0.5
    M[0,1] = M[1,0] = 0.5*d
        
    return M

###############################################################################
# compute_irdis_model_coefficient_matrix
###############################################################################

def compute_irdis_model_coefficient_matrix(p1, p2, a1, a2, theta_hwp1, theta_hwp2, theta_der1, theta_der2, dates, filter_used):
    '''
    Calculate coefficient matrix for the IRDIS polarimetry model correction

    Input:
        p1: array or list of mean parallactic angles of the first measurements used to compute the double sum and double difference (p^+)
        p2: array or list of mean parallactic angles of the second measurements used to compute the double sum and double difference (p^-)
        a1: array or list of mean altitude angles of the first measurements used to compute the double sum and double difference (a^+)
        a2: array or list of mean altitude angles of the second measurements used to compute the double sum and double difference (a^-)
        theta_hwp1: array or list of mean HWP angles of the first measurements used to compute the double sum and double difference (theta_hwp^+)
        theta_hwp2: array or list of mean HWP angles of the second measurements used to compute the double sum and double difference (theta_hwp^-)
        theta_der1: array or list of mean derotator angles of the first measurements used to compute the double sum and double difference (theta_der^+)
        theta_der2: array or list of mean derotator angles of the second measurements used to compute the double sum and double difference (theta_der^-)
        dates: array or list of dates of observations (format 'YYYY-MM-DD', e.g.: '2018-09-30')
        filter_used: string specifying filter used for observations ('FILT_BBF_Y', 'FILT_BBF_J', 'FILT_BBF_H' or 'FILT_BBF_Ks')

    Output:      
        X: coefficient matrix containing the double difference row vectors describing the measurements
        
    File written by Rob van Holstein
    Function status: verified
    '''   
                
    # Define model parameters that do not depend on the date of the observations
    delta_hwp = -0.613158589269
    delta_der = 0.500072483779
            
    if filter_used in ['FILT_BBF_Y', 'FILT_NBF_HeI']:
        Delta_UT = 171.891576898
        Delta_M4 = 171.891576898
        epsilon_hwp = -0.00021492286258857182
        Delta_hwp = 184.24489040687035
        epsilon_der = -0.0009423930026144313                                     
        Delta_der = 126.11957538036766
        d_CI = 0.9801695109615369
        filter_model = 'FILT_BBF_Y'
                                    
    elif filter_used in ['FILT_BBF_J', 'FILT_NBF_ContJ', 'FILT_NBF_PaB']:
        Delta_UT = 173.414169049                
        Delta_M4 = 173.414169049
        epsilon_hwp = -0.00043278581895049085
        Delta_hwp = 177.52027378388257
        epsilon_der = -0.008303978181252019                                      
        Delta_der = 156.0584333408133
        d_CI = 0.9894796343284551
        filter_model = 'FILT_BBF_J'
        
    elif filter_used in ['FILT_BBF_H', 'FILT_NBF_ContH', 'FILT_NBF_FeII']:
        Delta_UT = 174.998748608
        Delta_M4 = 174.998748608
        epsilon_hwp = -0.00029657803108325395
        Delta_hwp = 170.67214967707864
        epsilon_der = -0.002260131403393225                                       
        Delta_der = 99.32313652084311
        d_CI = 0.9955313968849352
        filter_model = 'FILT_BBF_H'
                                                                           
    elif filter_used in ['FILT_BBF_Ks', 'FILT_NBF_ContK1', 'FILT_NBF_H2', 'FILT_NBF_BrG', 'FILT_NBF_CntK2', 'FILT_NBF_CO']:
        Delta_UT = 176.302288996           
        Delta_M4 = 176.302288996
        epsilon_hwp = -0.00041456866069250524
        Delta_hwp = 177.61874393785442
        epsilon_der = 0.0035517563420643166                                        
        Delta_der = 84.13439892002613
        d_CI = 0.9841908773870153
        filter_model = 'FILT_BBF_Ks'

    printandlog('\nUsing model parameters corresponding to filter ' + filter_model[5:] + '.')
    
    # Define Mueller matrices of HWP, derotator and common path and IRDIS
    M_hwp = compute_reflection_mueller_matrix(epsilon_hwp, Delta_hwp)
    M_der = compute_reflection_mueller_matrix(epsilon_der, Delta_der)
    M_CId_left = compute_polarizer_mueller_matrix(d_CI)
    M_CId_right = compute_polarizer_mueller_matrix(-d_CI)

    # Define date of recoating of M1 and M3
    date_recoating = time.strptime('2017-04-16', "%Y-%m-%d")       

    # Print whether observations are taken before or after recoating or both
    dates_unique = [time.strptime(x, '%Y-%m-%d') for x in list(set(dates))]
    dates_before_recoating = [x < date_recoating for x in dates_unique]

    if all(dates_before_recoating):
        printandlog('\nUsing epsilon_UT and epsilon_M4 from before the recoating of M1 and M3 that took place on April 16 2017.')
    elif not any(dates_before_recoating):
        printandlog('\nUsing epsilon_UT and epsilon_M4 from after the recoating of M1 and M3 that took place on April 16 2017.')
    else:
        printandlog('\nUsing epsilon_UT and epsilon_M4 from both before and after the recoating of M1 and M3 that took place on April 16 2017.')

    # Convert observation dates to allow comparison to date of recoating
    dates = [time.strptime(x, '%Y-%m-%d') for x in dates]

    # Compute coefficient matrix describing the measurements
    n = len(p1) 
    X = np.zeros((n, 4))

    for i in range(n):
        
        # Define model parameters that depend on the date of the observations
        if filter_used in ['FILT_BBF_Y', 'FILT_NBF_HeI']:
            if dates[i] < date_recoating:
                epsilon_UT = 0.023607413903534567 
                epsilon_M4 = 0.018211735858186456
            else:
                epsilon_UT = 0.01745394681183012
                epsilon_M4 = 0.018194769704367342
            
        elif filter_used in ['FILT_BBF_J', 'FILT_NBF_ContJ', 'FILT_NBF_PaB']:
            if dates[i] < date_recoating:
                epsilon_UT = 0.016685701811847004
                epsilon_M4 = 0.012844478639635984
            else:
                epsilon_UT = 0.01213513552053676
                epsilon_M4 = 0.013046513475544473
   
        elif filter_used in ['FILT_BBF_H', 'FILT_NBF_ContH', 'FILT_NBF_FeII']:
            if dates[i] < date_recoating:
                epsilon_UT = 0.012930082215499676
                epsilon_M4 = 0.009845229155837451 
            else:
                epsilon_UT = 0.009032205030412622
                epsilon_M4 = 0.009220985704954044   
                
        elif filter_used in ['FILT_BBF_Ks', 'FILT_NBF_ContK1', 'FILT_NBF_H2', 'FILT_NBF_BrG', 'FILT_NBF_CntK2', 'FILT_NBF_CO']:
            if dates[i] < date_recoating:
                epsilon_UT = 0.010575041739439168
                epsilon_M4 = 0.007766430613637994 
            else:
                epsilon_UT = 0.0074593218421309315
                epsilon_M4 = 0.008078285701524354
  
        # Define Mueller matrices of UT and M4
        M_UT = compute_reflection_mueller_matrix(epsilon_UT, Delta_UT)
        M_M4 = compute_reflection_mueller_matrix(epsilon_M4, Delta_M4) 
    
        # Compute rotation Mueller matrices
        T_p1 = compute_rotation_mueller_matrix(p1[i])
        T_p2 = compute_rotation_mueller_matrix(p2[i])
        T_a1 = compute_rotation_mueller_matrix(a1[i])
        T_a2 = compute_rotation_mueller_matrix(a2[i])
        T_hwp1 = compute_rotation_mueller_matrix(theta_hwp1[i] + delta_hwp)
        T_hwp1min = compute_rotation_mueller_matrix(-(theta_hwp1[i] + delta_hwp))
        T_hwp2 = compute_rotation_mueller_matrix(theta_hwp2[i] + delta_hwp)
        T_hwp2min = compute_rotation_mueller_matrix(-(theta_hwp2[i] + delta_hwp))
        T_der1 = compute_rotation_mueller_matrix(theta_der1[i] + delta_der)
        T_der1min = compute_rotation_mueller_matrix(-(theta_der1[i] + delta_der))
        T_der2 = compute_rotation_mueller_matrix(theta_der2[i] + delta_der)
        T_der2min = compute_rotation_mueller_matrix(-(theta_der2[i] + delta_der))
        
        # Compute Mueller matrices of two measurements
        M_UTM41 = M_M4.dot(T_a1.dot(M_UT.dot(T_p1)))
        M_CPI1 = T_der1min.dot(M_der.dot(T_der1.dot(T_hwp1min.dot(M_hwp.dot(T_hwp1)))))
        M1_left = M_CId_left.dot(M_CPI1.dot(M_UTM41))
        M1_right = M_CId_right.dot(M_CPI1.dot(M_UTM41))
        
        M_UTM42 = M_M4.dot(T_a2.dot(M_UT.dot(T_p2)))
        M_CPI2 = T_der2min.dot(M_der.dot(T_der2.dot(T_hwp2min.dot(M_hwp.dot(T_hwp2)))))
        M2_left = M_CId_left.dot(M_CPI2.dot(M_UTM42))
        M2_right = M_CId_right.dot(M_CPI2.dot(M_UTM42))

        # Compute single differences row vectors
        SD1 = M1_left[0, :] - M1_right[0, :]
        SD2 = M2_left[0, :] - M2_right[0, :]
        
        # Compute double difference row vector
        DD = 0.5*(SD1 - SD2)

        # Assembly double difference row vectors in a coefficient matrix 
        X[i, :] = DD
        
    return X   

###############################################################################
# correct_instrumental_polarization_effects
###############################################################################

def correct_instrumental_polarization_effects(cube_I_Q_double_sum, 
                                              cube_I_U_double_sum, 
                                              cube_Q_double_difference, 
                                              cube_U_double_difference, 
                                              header, 
                                              file_index_object, 
                                              annulus_star, 
                                              annulus_background, 
                                              combination_method_polarization='least squares', 
                                              trimmed_mean_prop_to_cut_polar=0.1, 
                                              combination_method_intensity='mean', 
                                              trimmed_mean_prop_to_cut_intens=0.1, 
                                              single_posang_north_up=True):
    '''
    Calculate incident I_Q-, I_U-, Q- and U-images by correcting for the instrumental polarization effects of IRDIS using the polarimetric instrument model

    Input:
        cube_I_Q_double_sum: cube of double-sum intensity I_Q-images in order of HWP cycles
        cube_I_U_double_sum: cube of double-sum intensity I_U-images in order of HWP cycles  
        cube_Q_double_difference: cube of double-difference Stokes Q-images in order of HWP cycles
        cube_U_double_difference: cube of double-difference Stokes U-images in order of HWP cycles 
        header: list of FITS-headers of OBJECT-files
        file_index_object: list of file indices of OBJECT-files (0-based)
        annulus_star: (list of) length-6-tuple(s) with parameters to generate annulus to measure polarization of star:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
        annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
        combination_method_polarization: method to be used to produce the incident Q- and U-images, 
            'least squares', 'trimmed mean' or 'median' (default = 'least squares')
        trimmed_mean_prop_to_cut_polar: fraction to cut off of both tails of the distribution if 
            combination_method_polarization = 'trimmed mean' (default = 0.1) 
        combination_method_intensity: method to be used to produce the incident I_Q- and I_U-images, 
            'mean', 'trimmed mean' or 'median' (default = 'mean')
        trimmed_mean_prop_to_cut_intens: fraction to cut off of both tails of the distribution if 
            trimmed_mean_prop_to_cut_intens = 'trimmed mean' (default = 0.1) 
        single_posang_north_up: if True the images produced are oriented with North up; if False the images have the image orientation of the
            raw frames (default = True); only valid for observations taken in field-tracking mode with a single derotator 
            position angle; parameter is ignored for pupil-tracking observations or field-tracking observations with more 
            than one derotator position angle, because in these cases the final images produced always have North up

    Output:
        frame_I_Q_incident: polarimetric model-corrected incident I_Q-image
        frame_I_U_incident: polarimetric model-corrected incident I_U-image
        frame_Q_incident: polarimetric model-corrected incident Q-image
        frame_U_incident: polarimetric model-corrected incident U-image
        cube_I_Q_incident: cube of polarimetric model-corrected incident I_Q-images for each HWP cycle
        cube_I_U_incident: cube of polarimetric model-corrected incident I_U-images for each HWP cycle
        cube_Q_incident: cube of polarimetric model-corrected incident Q-images for each HWP cycle
        cube_U_incident: cube of polarimetric model-corrected incident U-images for each HWP cycle
           
    File written by Rob van Holstein
    Function status: verified
    '''    

    ###############################################################################
    # Read headers and extract (component) angles and date
    ###############################################################################

    # Determine filter used    
    filter_used = header[0]['ESO INS1 FILT ID']   
    
    # Determine tracking mode used
    tracking_mode_used = header[0]['ESO INS4 COMB ROT']
  
    if tracking_mode_used == 'FIELD':
        printandlog('\nComputing the derotator angle for field-tracking mode.')
    elif tracking_mode_used == 'PUPIL':
        printandlog('\nComputing the derotator angle for pupil-tracking mode.')
    
    # Read header values and compute component angles
    p = np.zeros(len(header))
    a_start = np.copy(p)
    theta_hwp = np.copy(p)
    theta_der = np.copy(p)
    mjd = np.copy(p)
    exposure_time = np.copy(p)
    NDIT = np.copy(p)
    derotator_position_angle = np.copy(p)   
    dates = np.empty(p.shape, dtype='object')
    
    for i, header_sel in enumerate(header):
        p[i] = compute_mean_angle([header_sel['ESO TEL PARANG START'], header_sel['ESO TEL PARANG END']])
        a_start[i] = header_sel['ESO TEL ALT']
        theta_hwp[i] = np.mod(compute_mean_angle([header_sel['ESO INS4 DROT3 BEGIN'], header_sel['ESO INS4 DROT3 END']]) - 152.15, 360)
        if tracking_mode_used == 'FIELD':
            theta_der[i] = np.mod(compute_mean_angle([header_sel['ESO INS4 DROT2 BEGIN'], header_sel['ESO INS4 DROT2 END']]), 360)
        elif tracking_mode_used == 'PUPIL':
            theta_der[i] = np.mod(compute_mean_angle([header_sel['ESO INS4 DROT2 BEGIN'], header_sel['ESO INS4 DROT2 END']]) + 0.5*pupil_offset, 360)    
        mjd[i] = header_sel['MJD-OBS']
        exposure_time[i] = header_sel['EXPTIME']
        NDIT[i] = header_sel['ESO DET NDIT']
        derotator_position_angle[i] = header_sel['ESO INS4 DROT2 POSANG']
        dates[i] = header_sel['DATE'][:10]
    
    # Define a Julian date with respect to noon at Paranal and round it to night number 
    # (Paranal is 70 deg in longitude or a fraction 0.2 away from Greenwich and mjd is
    # defined with respect to midnight at Greenwich)
    jd_paranal_noon = mjd - 0.7
    night_number = np.floor(jd_paranal_noon)
    
    # Determine indices to separate observations taken at different nights
    indices_cut = np.where(np.append(0.0, np.diff(night_number)) != 0.0)[0]
    indices_cut = np.append(indices_cut, len(mjd))
    
    # Print number of splines used for the interpolation
    number_of_nights = len(indices_cut)
    if number_of_nights == 1:
        printandlog('\nInterpolating the altitude angle using a single spline for this single night of observations.' )
    else:
        printandlog('\nInterpolating the altitude angles using a separate spline for each of the ' + str(number_of_nights) + ' nights.')
    
    # Calculate mean Julian date halfway each exposure
    file_execution_time = NDIT * (0.938 + exposure_time) + 2.4
    mjd_half = mjd + 0.5 * file_execution_time / msd
    
    # Compute mean altitude angle by interpolation (and extrapolation) using a separate spline for each night
    index_start = 0
    a = np.zeros(a_start.shape)
    
    for i, index_end in enumerate(indices_cut):   
        altitude_spline = interpolate.InterpolatedUnivariateSpline(mjd[index_start:index_end], a_start[index_start:index_end], k=3)
        a[index_start:index_end] = altitude_spline(mjd_half[index_start:index_end])
 
    # Split header values of first and second measurements used to compute the double sum and double difference
    stokes_parameter = np.array([x['ESO OCS DPI H2RT STOKES'] for x in header])     
    indices_QUplus = np.nonzero(np.logical_or(stokes_parameter == 'Qplus', stokes_parameter == 'Uplus'))[0]
    indices_QUminus = np.nonzero(np.logical_or(stokes_parameter == 'Qminus', stokes_parameter == 'Uminus'))[0]

    p1 = p[indices_QUplus]
    p2 = p[indices_QUminus]
    a1 = a[indices_QUplus]
    a2 = a[indices_QUminus]
    theta_hwp1 = theta_hwp[indices_QUplus]
    theta_hwp2 = theta_hwp[indices_QUminus]
    theta_der1 = theta_der[indices_QUplus]
    theta_der2 = theta_der[indices_QUminus]  
    file_index_object1 = np.array(file_index_object)[indices_QUplus]
    file_index_object2 = np.array(file_index_object)[indices_QUminus]
    
    # Plot header values as a function of time
    plot_name = name_file_root + 'header_angles.png'
    printandlog('\nCreating plot ' + plot_name + ' showing the parallactic, altitude, HWP and derotator angles of the observations.')
    angles = np.hstack([p1, p2, a1, a2, theta_hwp1, theta_hwp2, theta_der1, theta_der2])
    font_size = 10
    plt.figure(figsize = (12, 8))
    plt.plot(file_index_object1 + 1, p1, 'ob', label = 'Parallactic angle 1')    
    plt.plot(file_index_object2 + 1, p2, 'sb', label = 'Parallactic angle 2')    
    plt.plot(file_index_object1 + 1, a1, 'or', label = 'Altitude angle 1')       
    plt.plot(file_index_object2 + 1, a2, 'sr', label = 'Altitude angle 2')       
    plt.plot(file_index_object1 + 1, theta_hwp1, 'og', label = 'HWP angle 1')       
    plt.plot(file_index_object2 + 1, theta_hwp2, 'sg', label = 'HWP angle 2')    
    plt.plot(file_index_object1 + 1, theta_der1, 'ok', label = 'Derotator angle 1') 
    plt.plot(file_index_object2 + 1, theta_der2, 'sk', label = 'Derotator angle 2')
    ax = plt.gca()
    ax.set_xlabel(r'File number', fontsize = font_size)
    ax.tick_params(axis = 'x', labelsize = font_size)
    ax.set_ylim([np.min(angles) - 10, np.max(angles) + 10])
    ax.set_ylabel(r'Angle ($^\circ$)', fontsize = font_size)
    ax.tick_params(axis = 'y', labelsize = font_size)  
    ax.grid()
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.savefig(os.path.join(path_reduced_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.close()

    ###############################################################################
    # Compute model coefficient matrix, IP and cross-talk elements
    ###############################################################################

    # Construct coefficient matrix describing the measurements
    X = compute_irdis_model_coefficient_matrix(p1=p1, p2=p2, a1=a1, a2=a2, theta_hwp1=theta_hwp1, theta_hwp2=theta_hwp2,
                                               theta_der1=theta_der1, theta_der2=theta_der2, dates=dates, filter_used=filter_used) 
 
    # Obtain elements describing IP and cross-talk/polarized transmission/rotation of double difference and assume V_in = 0
    IP = X[:, 0]
    X_QU = X[:, 1:3]

    # Obtain indices of Q and U double differences
    indices_Q = np.nonzero(stokes_parameter[np.logical_or(stokes_parameter == 'Qplus', stokes_parameter == 'Uplus')] == 'Qplus')[0]
    indices_U = np.nonzero(stokes_parameter[np.logical_or(stokes_parameter == 'Qplus', stokes_parameter == 'Uplus')] == 'Uplus')[0]
    
    # Extract intrumental polarization of Q- and U-measurements and polarized transmission and cross-talk/rotation elements
    IP_Q = IP[indices_Q]
    IP_U = IP[indices_U]
    QQ_Q = X_QU[indices_Q, 0]
    UQ_Q = X_QU[indices_Q, 1]
    QQ_U = X_QU[indices_U, 0]
    UQ_U = X_QU[indices_U, 1]

    ###############################################################################
    # Compute rotation angles of images
    ###############################################################################
    
    # Determine number of unique derotator position angles
    number_derotator_position_angles = len(np.unique(derotator_position_angle))
       
    # Compute indices of Qplus, Qminus, Uplus and Uminus
    indices_Qplus = np.nonzero(stokes_parameter == 'Qplus')[0]
    indices_Qminus = np.nonzero(stokes_parameter == 'Qminus')[0]
    indices_Uplus = np.nonzero(stokes_parameter == 'Uplus')[0]
    indices_Uminus = np.nonzero(stokes_parameter == 'Uminus')[0] 

    if tracking_mode_used == 'FIELD' and number_derotator_position_angles == 1 and single_posang_north_up == False:
        # Set rotation angles of images equal to zero
        printandlog('\nSetting the image rotation angles for field-tracking mode equal to zero so that North is not up in the final images.')
        rotation_angles_Q = np.zeros(len(indices_Qplus))
        rotation_angles_U = np.zeros(len(indices_Uplus))

    elif tracking_mode_used == 'FIELD':
        # Compute rotation angle of images
        printandlog('\nComputing the image rotation angles for field-tracking mode to orient the images with North up. The true North correction used equals %.2f deg.' % true_north_correction)
        rotation_angles_Q = -derotator_position_angle[indices_Qplus] - true_north_correction
        rotation_angles_U = -derotator_position_angle[indices_Uplus] - true_north_correction
        
    elif tracking_mode_used == 'PUPIL':    
        # Compute mean parallactic angles of first and second measurements used to compute I_Q and Q, and I_U and U
        printandlog('\nComputing the image rotation angles for pupil-tracking mode to orient the images with North up. The true North correction used equals %.2f deg.' % true_north_correction)
        p_Q = compute_mean_angle(np.vstack([p[indices_Qplus], p[indices_Qminus]]), axis = 0)
        p_U = compute_mean_angle(np.vstack([p[indices_Uplus], p[indices_Uminus]]), axis = 0)
        
        # Compute rotation angle of images
        rotation_angles_Q = -p_Q - pupil_offset - true_north_correction
        rotation_angles_U = -p_U - pupil_offset - true_north_correction

    ###############################################################################
    # Measure polarization signal in annulus on star
    ###############################################################################

    # If the annuli for the star and background are centered on the star and rotationally symmetric
    if np.all(np.array(annulus_star, ndmin=2)[:, np.array([0, 1, 4, 5])] == np.array([511.5, 511.5, 0, 360])) and \
       np.all(np.array(annulus_background, ndmin=2)[:, np.array([0, 1, 4, 5])] == np.array([511.5, 511.5, 0, 360])):
        # Do not rotate the cubes of double-sum and double-difference images and compute normalized Stokes q and u in an annulus   
        printandlog('\nNot rotating the images used to determine the polarization signal in an annulus because the annuli used for the star and the background are centered on the star and rotationally symmetric.')

        q_annulus, u_annulus = determine_star_polarization(cube_I_Q=cube_I_Q_double_sum, 
                                                           cube_I_U=cube_I_U_double_sum, 
                                                           cube_Q=cube_Q_double_difference, 
                                                           cube_U=cube_U_double_difference, 
                                                           annulus_star=annulus_star, 
                                                           annulus_background=annulus_background)

    else:
        # Rotate the cubes of double-sum and double-difference images so that the annuli are at the right position
        printandlog('\nRotating the images used to determine the polarization signal in annulus because the annuli used for the star and the background are not centered on the star and/or not rotationally symmetric.')
        cube_I_Q_annulus = np.zeros(cube_I_Q_double_sum.shape)
        cube_Q_annulus = np.zeros(cube_Q_double_difference.shape)
        cube_I_U_annulus = np.zeros(cube_I_U_double_sum.shape)
        cube_U_annulus = np.zeros(cube_U_double_difference.shape)
        
        for i, (frame_I_Q, frame_Q, rotation_angle_Q) in enumerate(zip(cube_I_Q_double_sum, cube_Q_double_difference, rotation_angles_Q)):
            cube_I_Q_annulus[i, :, :] = rotate(frame_I_Q, rotation_angle_Q, reshape=False)
            cube_Q_annulus[i, :, :] = rotate(frame_Q, rotation_angle_Q, reshape=False)
            
        for i, (frame_I_U, frame_U, rotation_angle_U) in enumerate(zip(cube_I_U_double_sum, cube_U_double_difference, rotation_angles_U)):
            cube_I_U_annulus[i, :, :] = rotate(frame_I_U, rotation_angle_U, reshape=False)
            cube_U_annulus[i, :, :] = rotate(frame_U, rotation_angle_U, reshape=False)
         
        # Compute normalized Stokes q and u in an annulus   
        q_annulus, u_annulus = determine_star_polarization(cube_I_Q=cube_I_Q_annulus, 
                                                           cube_I_U=cube_I_U_annulus, 
                                                           cube_Q=cube_Q_annulus, 
                                                           cube_U=cube_U_annulus, 
                                                           annulus_star=annulus_star, 
                                                           annulus_background=annulus_background)
    
    ###############################################################################
    # Fit polarization of star from measurements using model coefficient matrix
    ###############################################################################
         
    # Create y-data array                               
    ydata = np.zeros(len(np.hstack([q_annulus, u_annulus])))
    ydata[indices_Q] = q_annulus
    ydata[indices_U] = u_annulus
  
    # Define function that computes the sum of squared residuals of fit
    def residuals(parameters, ydata, X):
        return np.sum((ydata - X.dot(np.array([1, parameters[0], parameters[1], 0])))**2) 

    # Fit data to model
    fit_result = optimize.minimize(residuals, np.array([0, 0]), args=(ydata, X), method='SLSQP', 
                                   bounds=[(-1,1), (-1,1)], tol=1e-14, options={'maxiter':100, 'disp':False})
 
    # Retrieve normalized Stokes q and u and compute DoLP and AoLP of star
    q_star, u_star = fit_result.x
    DoLP_star = np.sqrt(q_star**2 + u_star**2)
    AoLP_star = np.mod(np.rad2deg(0.5 * np.arctan2(u_star, q_star)), 180)
    
    # Print resulting star polarization
    printandlog('\nFitted star polarization:')
    printandlog('q_star =    %.4f %%' % (100*q_star))
    printandlog('u_star =    %.4f %%' % (100*u_star))
    printandlog('DoLP_star = %.4f %%' % (100*DoLP_star))
    printandlog('AoLP_star = %.2f deg' % AoLP_star)
           
    # Obtain fitted q and u and determine residuals
    yfit = X.dot(np.array([1, q_star, u_star, 0])) 
    q_fit = yfit[indices_Q]
    u_fit = yfit[indices_U]
    res_Q = q_annulus - q_fit
    res_U = u_annulus - u_fit
    
    # Compute and print coefficient of determination (R squared)
    rsq = 1 - np.sum((ydata - yfit)**2) / (np.var(ydata)*len(ydata))  
    printandlog('R squared of fit = %.4f' % rsq)
    
    ###############################################################################
    # Plot model-predicted IP, star polarization and cross-talk elements vs HWP cycle number
    ###############################################################################

    # Plot model, measured and fitted normalized Stokes parameter vs HWP cycle number 
    plot_name = name_file_root + 'model_ip_star_pol.png'
    printandlog('\nCreating plot ' + plot_name + ' showing the model-predicted IP, measured polarization signal and model-predicted IP + fitted star polarization vs. HWP cycle number.')
    font_size = 10  
    marker_size = 4
    x_max = max([len(IP_Q), len(IP_U)])
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[5, 2]}, sharex=True, figsize = (9.4, 10.34))  
    ax1.plot([0, x_max + 1],[0, 0], '-k', zorder=1) 
    ph1, = ax1.plot(np.arange(1, len(IP_Q) + 1), 100*IP_Q, '--b', markersize = marker_size, label = '$q$ model', zorder=5)  
    ph2, = ax1.plot(np.arange(1, len(IP_U) + 1), 100*IP_U, '--r', markersize = marker_size, label = '$u$ model', zorder=4)
    ph3, = ax1.plot(np.arange(1, len(q_annulus) + 1), 100*q_annulus, 'ob', markersize = marker_size, label = '$q$ meas.', zorder=9)  
    ph4, = ax1.plot(np.arange(1, len(u_annulus) + 1), 100*u_annulus, 'or', markersize = marker_size, label = '$u$ meas.', zorder=8) 
    ph5, = ax1.plot(np.arange(1, len(q_fit) + 1), 100*q_fit, '-b', markersize = marker_size, label = '$q$ model + star fit.', zorder=7)  
    ph6, = ax1.plot(np.arange(1, len(u_fit) + 1), 100*u_fit, '-r', markersize = marker_size, label = '$u$ model + star fit.', zorder=6)  
    ax1.set_ylabel(r'Normalized Stokes parameter (%)', fontsize = font_size)
    ax1.tick_params(axis = 'y', labelsize = font_size)
    ax1.grid()
    l1 = ax1.legend([ph1, ph3, ph5], ['$q$ model', '$q$ meas.', '$q$ model + $q/u_\mathrm{in}$ fit.'], loc = 'center left', prop={'size': font_size}, handlelength = 2.5, ncol = 1)
    ax1.legend([ph2, ph4, ph6], ['$u$ model', '$u$ meas.', '$u$ model + $q/u_\mathrm{in}$ fit.'], loc = 'lower right', prop={'size': font_size}, handlelength = 2.5, ncol = 1)
    ax1.add_artist(l1)
    ax2.plot([0, x_max + 1],[0, 0], '-k', zorder=1) 
    ax2.plot(np.arange(1, len(res_Q) + 1), 100*res_Q, 'ob', markersize = marker_size, label = '$q$ model', zorder=5)  
    ax2.plot(np.arange(1, len(res_U) + 1), 100*res_U, 'or', markersize = marker_size, label = '$u$ model', zorder=4)
    ax2.set_xlabel(r'HWP cycle', fontsize = font_size)
    ax2.tick_params(axis = 'x', labelsize = font_size)
    ax2.set_xlim([0, x_max + 1])
    ax2.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    ax2.set_ylabel(r'Residuals (%)', fontsize = font_size)
    ax2.tick_params(axis = 'y', labelsize = font_size)
    ax2.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(path_reduced_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.close()

    # Plot elements Q->Q and U->Q from model as a function of HWP cycle number
    plot_name = name_file_root + 'model_crosstalk_transmission.png'
    printandlog('\nCreating plot ' + plot_name + ' showing the model-predicted polarized transmission and crosstalk/rotation elements vs. HWP cycle number.')
    font_size = 10
    x_max = max([len(QQ_Q), len(QQ_U)])
    plt.figure(figsize = (5.9, 3.8))
    plt.plot([0, x_max + 1],[0, 0], '-k')  
    plt.plot(np.arange(1, len(QQ_Q) + 1), QQ_Q, 'o-b', label = r'Qin $\rightarrow$ Q')  
    plt.plot(np.arange(1, len(UQ_Q) + 1), UQ_Q, 'o-', color = (0.5, 0, 0.5), label = r'Uin $\rightarrow$ Q')
    plt.plot(np.arange(1, len(QQ_U) + 1), QQ_U, 'o-', color = (1, 0.5, 0), label = r'Qin $\rightarrow$ U')
    plt.plot(np.arange(1, len(UQ_U) + 1), UQ_U, 'o-r', label = r'Uin $\rightarrow$ U')  
    ax = plt.gca()
    ax.set_xlabel(r'HWP cycle', fontsize = font_size)
    ax.tick_params(axis = 'x', labelsize = font_size)
    ax.set_xlim([0, x_max + 1])
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    ax.set_ylabel(r'Retardance/rotation elements from model', fontsize = font_size)
    ax.tick_params(axis = 'y', labelsize = font_size)  
    ax.grid()
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.savefig(os.path.join(path_reduced_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.close()
      
    ###############################################################################
    # Compute incident Q- and U-images by correcting for instrumental polarization effects
    ###############################################################################

    # Subtract intrumental polarization from Q- and U-images
    cube_Q_IP_subtracted = cube_Q_double_difference - IP_Q[:, np.newaxis, np.newaxis]*cube_I_Q_double_sum
    cube_U_IP_subtracted = cube_U_double_difference - IP_U[:, np.newaxis, np.newaxis]*cube_I_U_double_sum

#TODO: remove simple cADI below
#    cube_Q_IP_subtracted -= np.median(cube_Q_IP_subtracted, axis=0)
#    cube_U_IP_subtracted -= np.median(cube_U_IP_subtracted, axis=0)
    
    # Derotate IP-subtracted Q- and U-images
    cube_Q_derotated = np.zeros(cube_Q_IP_subtracted.shape)
    cube_U_derotated = np.zeros(cube_U_IP_subtracted.shape)

    for i, (frame_Q, rotation_angle_Q) in enumerate(zip(cube_Q_IP_subtracted, rotation_angles_Q)):
        cube_Q_derotated[i, :, :] = rotate(frame_Q, rotation_angle_Q, reshape=False)
    
    for i, (frame_U, rotation_angle_U) in enumerate(zip(cube_U_IP_subtracted, rotation_angles_U)):
        cube_U_derotated[i, :, :] = rotate(frame_U, rotation_angle_U, reshape=False)
         
    # Calculate Q- and U-images incident on telescope by solving system of equations per HWP cycle
    number_frames_min = min(len(cube_Q_derotated), len(cube_U_derotated))
    cube_Q_incident = np.zeros((number_frames_min,) + cube_Q_derotated.shape[-2:])
    cube_U_incident = np.copy(cube_Q_incident)

    for i, (frame_Q, index_Q, frame_U, index_U) in enumerate(zip(cube_Q_derotated, indices_Q, cube_U_derotated, indices_U)):
        # Apply model correction of Q->Q and U->Q elements by solving equations per HWP cycle
        Y = np.stack([frame_Q, frame_U])    
        Y_stretched = Y.reshape(Y.shape[0], Y.shape[1]*Y.shape[2])
        X_hwp_cycle = np.stack([X_QU[index_Q], X_QU[index_U]])
        cube_QU_incident_stretched = np.linalg.solve(X_hwp_cycle, Y_stretched)
        cube_QU_incident = cube_QU_incident_stretched.reshape(cube_QU_incident_stretched.shape[0], Y.shape[1], Y.shape[2]) 
        cube_Q_incident[i, :, :] = cube_QU_incident[0, :, :]
        cube_U_incident[i, :, :] = cube_QU_incident[1, :, :]
        
    # Create incident Q- and U-images
    if combination_method_polarization == 'least squares':
        # Obtain incident Q- and U-images from the least squares solution
        printandlog('\nComputing the incident Q- and U-images using least squares.')
          
        # Create cube Y with Q/U measurements in order of HWP cycles
        Y = np.zeros(((len(p1),) + cube_Q_derotated.shape[-2:]))
        Y[indices_Q, :, :] = cube_Q_derotated
        Y[indices_U, :, :] = cube_U_derotated
    
        # Make images 1D, compute linear least-squares solution and reshape final images to 2D again   
        Y_stretched = Y.reshape(Y.shape[0], Y.shape[1]*Y.shape[2])     
        QU_in_stretched = np.linalg.lstsq(X_QU, Y_stretched, rcond=None)[0]
        QU_in = QU_in_stretched.reshape(QU_in_stretched.shape[0], Y.shape[1], Y.shape[2])    
        frame_Q_incident = QU_in[0, :, :] 
        frame_U_incident = QU_in[1, :, :]
    
    elif combination_method_polarization == 'trimmed mean': 
        # Compute incident Q- and U-images from the trimmed mean of incident cubes
        printandlog('\nComputing the incident Q- and U-images using the trimmed mean with a proportion to cut equal to ' + str(trimmed_mean_prop_to_cut_polar) + '.')
        frame_Q_incident = trim_mean(cube_Q_incident, proportiontocut=trimmed_mean_prop_to_cut_polar, axis=0)
        frame_U_incident = trim_mean(cube_U_incident, proportiontocut=trimmed_mean_prop_to_cut_polar, axis=0)
    
    elif combination_method_polarization == 'median': 
        # Compute incident Q- and U-images from the median of incident cubes 
        printandlog('\nComputing the incident Q- and U-images using the median.')
        frame_Q_incident = np.median(cube_Q_incident, axis=0)
        frame_U_incident = np.median(cube_U_incident, axis=0)

    ###############################################################################
    # Compute incident I_Q- and I_U-images
    ###############################################################################

#TODO: remove simple cADI below
#    cube_I_Q_double_sum -= np.median(cube_I_Q_double_sum, axis=0)
#    cube_I_U_double_sum -= np.median(cube_I_U_double_sum, axis=0)
    
    # Derotate I_Q- and I_U-images
    cube_I_Q_incident = np.zeros(cube_I_Q_double_sum.shape)
    cube_I_U_incident = np.zeros(cube_I_U_double_sum.shape)
    
    for i, (frame_I_Q, rotation_angle_Q) in enumerate(zip(cube_I_Q_double_sum, rotation_angles_Q)):
        cube_I_Q_incident[i, :, :] = rotate(frame_I_Q, rotation_angle_Q, reshape=False)
    
    for i, (frame_I_U, rotation_angle_U) in enumerate(zip(cube_I_U_double_sum, rotation_angles_U)):
        cube_I_U_incident[i, :, :] = rotate(frame_I_U, rotation_angle_U, reshape=False)

    # Create incident I_Q- and I_U-images
    if combination_method_intensity == 'mean':
        # Compute incident I_Q- and I_U-images from the mean 
        printandlog('\nComputing the incident I_Q- and I_U-images using the mean.')
        frame_I_Q_incident = np.mean(cube_I_Q_incident, axis = 0)
        frame_I_U_incident = np.mean(cube_I_U_incident, axis = 0)  
    
    elif combination_method_intensity == 'trimmed mean': 
        # Compute incident I_Q- and I_U-images from the trimmed mean
        printandlog('\nComputing the incident I_Q- and I_U-images using the trimmed mean with a proportion to cut equal to ' + str(trimmed_mean_prop_to_cut_intens) + '.')
        frame_I_Q_incident = trim_mean(cube_I_Q_incident, proportiontocut=trimmed_mean_prop_to_cut_intens, axis=0)
        frame_I_U_incident = trim_mean(cube_I_U_incident, proportiontocut=trimmed_mean_prop_to_cut_intens, axis=0)
    
    elif combination_method_intensity == 'median': 
        # Compute incident I_Q- and I_U-images from the median
        printandlog('\nComputing the incident I_Q- and I_U-images using the median.')
        frame_I_Q_incident = np.median(cube_I_Q_incident, axis=0)
        frame_I_U_incident = np.median(cube_I_U_incident, axis=0)

    # Make cubes of I_Q- and I_U-images the same length  
    cube_I_Q_incident = cube_I_Q_incident[:number_frames_min, :, :]
    cube_I_U_incident = cube_I_U_incident[:number_frames_min, :, :]
    
    return frame_I_Q_incident, frame_I_U_incident, frame_Q_incident, frame_U_incident, cube_I_Q_incident, cube_I_U_incident, cube_Q_incident, cube_U_incident

###############################################################################
# compute_azimuthal_stokes_parameters
###############################################################################
    
def compute_azimuthal_stokes_parameters(frame_Q, frame_U, rotation_angle=0, center_coordinates=None):
    '''  
    Compute images of azimuthal stokes parameters Q_phi and U_phi using defintions of
    van Holstein et al. (2019)
    
    Input:
        frame_Q: Stokes Q-image
        frame_U: Stokes U-image
        rotation_angle: additional image rotation (deg)
        center_coordinates: tuple containing center coordinates in x and y; if
            None, use center of frame
        
    Output:
        Q_phi: azimuthal Stokes Q_phi-image
        U_phi: azimuthal Stokes U_phi-image
        phi: frame showing values of the angle phi (deg)
        
    Note:
        Q_phi > 0 for azimuthal polarization and Q_phi < 0 for radial polarization
    
    File written by Rob van Holstein
    Function status: verified
    ''' 

    # Create grid and compute angle 
    x = np.arange(0, frame_Q.shape[-1])
    y = np.arange(0, frame_Q.shape[-2])
    xm, ym = np.meshgrid(x, y)
    
    if center_coordinates is None:
        x_center = 0.5*x[-1]
        y_center = 0.5*y[-1]
    else:
        x_center = center_coordinates[0]
        y_center = center_coordinates[1]

    phi = np.arctan2((x_center - xm), (ym - y_center)) + np.deg2rad(rotation_angle)
    
    # Compute Q_phi- and U_phi-images
    frame_Q_phi = -frame_Q*np.cos(2*phi) - frame_U*np.sin(2*phi)
    frame_U_phi = frame_Q*np.sin(2*phi) - frame_U*np.cos(2*phi)
    
    # Convert frame showing phi to degrees
    phi = np.rad2deg(phi)
    
    return frame_Q_phi, frame_U_phi, phi

###############################################################################
# compute_final_images
###############################################################################
    
def compute_final_images(frame_I_Q, frame_I_U, frame_Q, frame_U, header, single_posang_north_up=True):
    ''' 
    Compute final images of polarimetric data reduction
    
    Input:
        frame_I_Q: I_Q-image
        frame_I_U: I_U-image
        frame_Q: Q-image
        frame_U: U-image
        header: list of headers of raw science frames; if None, no Q_phi- and U_phi-images are created and they are set to None in the output     
        single_posang_north_up: if True the images produced are oriented with North up; if False the images have the image orientation of the
            raw frames (default = True); only valid for observations taken in field-tracking mode with a single derotator 
            position angle; parameter is ignored for pupil-tracking observations or field-tracking observations with more 
            than one derotator position angle, because in these cases the final images produced always have North up
        
    Output:
        frame_I_tot: total-intensity image
        frame_Q_phi: image of Q_phi (None if header is set to None)
        frame_U_phi: image of U_phi (None if header is set to None)
        frame_I_pol: polarized-intensity image
        frame_AoLP: image of angle of linear polarization computed from Q- and U-images
        frame_DoLP: image of degree of linear polarization computed from Q-, U-, I_Q- and I-U-images
        frame_q: image of normalized Stokes q
        frame_u: image of normalized Stokes u
        frame_AoLP_norm: image of angle of linear polarization computed using the normalized q- and u-images
        frame_DoLP_norm: image of degree of linear polarization computed using the normalized q- and u-images

    File written by Rob van Holstein
    Function status: verified         
    '''

    # Compute total intensity image  
    frame_I_tot = 0.5*(frame_I_Q + frame_I_U)
    
    if header != None:
        # Determine tracking mode used
        tracking_mode_used = header[0]['ESO INS4 COMB ROT']
        
        # Determine number of unique derotator position angles
        derotator_position_angle = [x['ESO INS4 DROT2 POSANG'] for x in header]
        number_derotator_position_angles = len(np.unique(derotator_position_angle))
    
        if tracking_mode_used == 'FIELD' and number_derotator_position_angles == 1 and single_posang_north_up == False:
            # Compute image position angle
            rotation_angle = -derotator_position_angle[0] - true_north_correction
        else:
            rotation_angle = 0.0
        
        frame_Q_phi, frame_U_phi, frame_azimuthal_angle = compute_azimuthal_stokes_parameters(frame_Q, frame_U, rotation_angle=-rotation_angle)   
    else:
        frame_Q_phi = None
        frame_U_phi = None

    # Compute polarized intensity image
    frame_I_pol = np.sqrt(frame_Q**2 + frame_U**2)
    
    # Compute images of angle and degree of linear polarization from Q-, U-, I_Q- and I-U-images
    frame_AoLP = np.mod(np.rad2deg(0.5 * np.arctan2(frame_U, frame_Q)), 180)
    frame_DoLP = frame_I_pol / frame_I_tot

    # Compute corrected q- and u-images
    frame_q = frame_Q / frame_I_Q
    frame_u = frame_U / frame_I_U

    # Compute images of angle and degree of linear polarization using the normalized q- and u-images
    frame_AoLP_norm = np.mod(np.rad2deg(0.5 * np.arctan2(frame_u, frame_q)), 180)
    frame_DoLP_norm = np.sqrt(frame_q**2 + frame_u**2)

    return frame_I_tot, frame_Q_phi, frame_U_phi, frame_I_pol, frame_AoLP, frame_DoLP, frame_q, frame_u, frame_AoLP_norm, frame_DoLP_norm

###############################################################################
# perform_postprocessing
###############################################################################

def perform_postprocessing(cube_left_frames, 
                           cube_right_frames, 
                           header, 
                           file_index_object, 
                           annulus_star='automatic', 
                           annulus_background='large annulus', 
                           double_difference_type='conventional', 
                           remove_vertical_band_detector_artefact=True, 
                           combination_method_polarization='least squares', 
                           trimmed_mean_prop_to_cut_polar=0.1, 
                           combination_method_intensity='mean', 
                           trimmed_mean_prop_to_cut_intens=0.1, 
                           single_posang_north_up=True, 
                           normalized_polarization_images=False, \
                           reference_flux=None):
    '''
    Perform post-processing of data, including applying the model-based correction
    for the instrumental polarization effects, and save final images to FITS-files
    
    Input:
        cube_left_frames: cube of pre-processed left frames 
        cube_right_frames: cube of pre-processed right frames 
        header: list of FITS-headers of OBJECT-files 
        file_index_object: list of file indices of OBJECT-files (0-based)
        annulus_star: (list of) length-6-tuple(s) with parameters to generate annulus to measure polarization of star:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            If string 'ao residuals' the annulus will be automatically determined and 
            be star-centered and located over the AO residuals. The inner and 
            outer radius of the annulus will depend on the filter used. If 
            'star aperture' a small aparture located at the position of
            the central star will be used. If 'automatic', annulus_star will 
            first be set to 'ao residuals' in case of coronagraphic data, and to 
            'star aperture' in case of non-coronagraphic data (default = 'automatic').
        annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            outer_radius: outer_radius (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg to the right and positive rotation counterclockwise)
            If string 'large annulus' the annulus will be star-centered and 
            located far away from the star with an inner radius of 360 pixels
            and an outer radius of 420 pixels (default = 'large annulus').
        double_difference_type: type of double difference to be computed, either
        'conventional' or 'normalized' (see van Holstein et al. 2019; default = 'conventional')
        remove_vertical_band_detector_artefact: If True remove the vertical band detector artefact seen in 
            the double-difference Q- and U-images. If False don't remove it (default = True).
        combination_method_polarization: method to be used to produce the incident Q- and U-images, 
            'least squares', 'trimmed mean' or 'median' (default = 'least squares')
        trimmed_mean_prop_to_cut_polar: fraction to cut off of both tails of the distribution if 
            combination_method_polarization = 'trimmed mean' (default = 0.1) 
        combination_method_intensity: method to be used to produce the incident I_Q- and I_U-images, 
            'mean', 'trimmed mean' or 'median' (default = 'mean')
        trimmed_mean_prop_to_cut_intens: fraction to cut off of both tails of the distribution if 
            trimmed_mean_prop_to_cut_intens = 'trimmed mean' (default = 0.1) 
        single_posang_north_up: if True the images produced are oriented with North up; if False the images have the image orientation of the
            raw frames (default = True); only valid for observations taken in field-tracking mode with a single derotator 
            position angle; parameter is ignored for pupil-tracking observations or field-tracking observations with more 
            than one derotator position angle, because in these cases the final images produced always have North up
        normalized_polarization_images: if True create final images of degree of linear polarization, normalized Stokes q and u
            and degree and angle of linear polarization computed from q and u; such images only have meaning if all flux in the images
            originates from the astrophysical source of interest (default = False)
        reference_flux: if the option convert_in_contrast_per_arcsec2 is set to True,
            then the image is converted in contrast by using the number of counts 
            for the star specified by this value.
        
    File written by Rob van Holstein
    Function status: verified      
    '''

    ###############################################################################
    # Define annulus for star polarization and background
    ###############################################################################
    
    printandlog('\n###############################################################################')
    printandlog('# Defining the annulus for star polarization and background')
    printandlog('###############################################################################') 
                
    # Determine filter and coronagraph used    
    filter_used = header[0]['ESO INS1 FILT ID']   
    coronagraph_used = header[0]['ESO INS COMB ICOR']

    # If annulus_star is 'automatic', set to 'ao residuals' for coronagraphic data and to 'star aperture' for non-coronagraphic data
    if annulus_star == 'automatic':
        if coronagraph_used == 'N_NS_CLEAR':
            printandlog('\nannulus_star is \'automatic\': changing it to \'star aperture\', because the data is non-coronagraphic.')
            annulus_star = 'star aperture'
        else:
            printandlog('\nannulus_star is \'automatic\': changing it to \'ao residuals\', because the data is coronagraphic.')
            annulus_star = 'ao residuals'

    # Define and print annulus to determine the star polarization from
    if type(annulus_star) == tuple or type(annulus_star) == list:
        printandlog('\nThe star polarization will be determined with a user-defined annulus or several user-defined annuli:')
        if type(annulus_star) == tuple:
            printandlog(annulus_0_to_1_based(annulus_star))
        elif type(annulus_star) == list:
            for x in annulus_star:
                printandlog(annulus_0_to_1_based(x))
    elif annulus_star == 'ao residuals':
        if filter_used in ['FILT_BBF_Y', 'FILT_NBF_HeI']:
            annulus_star = (511.5, 511.5, 40, 65, 0, 360)
        elif filter_used in ['FILT_BBF_J', 'FILT_NBF_ContJ', 'FILT_NBF_PaB']:
            annulus_star = (511.5, 511.5, 45, 75, 0, 360)
        elif filter_used in ['FILT_BBF_H', 'FILT_NBF_ContH', 'FILT_NBF_FeII']:
            annulus_star = (511.5, 511.5, 60, 95, 0, 360)
        elif filter_used in ['FILT_BBF_Ks', 'FILT_NBF_ContK1', 'FILT_NBF_H2', 'FILT_NBF_BrG', 'FILT_NBF_CntK2', 'FILT_NBF_CO']:
            annulus_star = (511.5, 511.5, 80, 120, 0, 360)  
        printandlog('\nThe star polarization will be determined with a star-centered annulus located over the AO residuals:')
        printandlog(annulus_0_to_1_based(annulus_star))
        if coronagraph_used == 'N_NS_CLEAR':
            printandlog('\nWARNING, the data is non-coronagraphic so there might be little flux at the AO residuals. Determining the star polarization using an aperture at the position of the central star (\'star aperture\') will probably yield better results.')
    elif annulus_star == 'star aperture':
        annulus_star = (511.5, 511.5, 0, 11, 0, 360)
        printandlog('\nThe star polarization will be determined with an aparture located at the position of the central star:')
        printandlog(annulus_0_to_1_based(annulus_star))
    
    # Define and print annulus to determine the background from
    if type(annulus_background) == tuple or type(annulus_background) == list:
        printandlog('\nThe background will be determined with a user-defined annulus or several user-defined annuli:')
        if type(annulus_background) == tuple:
            printandlog(annulus_0_to_1_based(annulus_background))
        elif type(annulus_background) == list:
            for x in annulus_background:
                printandlog(annulus_0_to_1_based(x))
    elif annulus_background == 'large annulus':
        annulus_background = (511.5, 511.5, 360, 420, 0, 360)
        printandlog('\nThe background will be determined with a star-centered annulus located far away from the central star:')
        printandlog(annulus_0_to_1_based(annulus_background))
    
    ###############################################################################
    # Compute double sum and difference
    ###############################################################################
    
    printandlog('\n###############################################################################')
    printandlog('# Computing the double sum and double difference')
    printandlog('###############################################################################') 
          
    # Compute single sum and single difference cubes   
    cube_single_sum = cube_left_frames + cube_right_frames
    cube_single_difference = cube_left_frames - cube_right_frames 
    
    # Compute double sum and double difference cubes       
    cube_I_Q_double_sum, cube_I_U_double_sum, cube_Q_double_difference, cube_U_double_difference \
    = compute_double_sum_double_difference(cube_single_sum, cube_single_difference, header, double_difference_type=double_difference_type)
    
    ###############################################################################
    # Remove vertical band detector artefact 
    ###############################################################################
    
    if remove_vertical_band_detector_artefact == True:
        # Remove IRDIS vertical band detector artefact seen in double-difference Q- and U-images
        number_pixels_artefact = 60
        printandlog('\nRemoving the vertical band detector artefact in the double-difference Q- and U-images by subtracting the median of the top and bottom ' + str(number_pixels_artefact) + ' pixels of each column.')
        cube_Q_artefact_removed = remove_detector_artefact(cube=cube_Q_double_difference, number_pixels=number_pixels_artefact)
        cube_U_artefact_removed = remove_detector_artefact(cube=cube_U_double_difference, number_pixels=number_pixels_artefact)
    elif remove_vertical_band_detector_artefact == False:
        printandlog('\nNot removing the vertical band detector artefact.')
        cube_Q_artefact_removed = cube_Q_double_difference
        cube_U_artefact_removed = cube_U_double_difference
       
    ###############################################################################
    # Correct for instrumental polarization effects using the polarimetric instrument model
    ###############################################################################
    
    printandlog('\n###############################################################################')
    printandlog('# Correcting images for instrumental polarization effects')
    printandlog('###############################################################################') 

    frame_I_Q_incident, frame_I_U_incident, frame_Q_incident, frame_U_incident, cube_I_Q_incident, cube_I_U_incident, cube_Q_incident, cube_U_incident \
    = correct_instrumental_polarization_effects(cube_I_Q_double_sum=cube_I_Q_double_sum, 
                                                cube_I_U_double_sum=cube_I_U_double_sum, 
                                                cube_Q_double_difference=cube_Q_artefact_removed, 
                                                cube_U_double_difference=cube_U_artefact_removed, 
                                                header=header, 
                                                file_index_object=file_index_object, 
                                                annulus_star=annulus_star, 
                                                annulus_background=annulus_background, 
                                                combination_method_polarization=combination_method_polarization, 
                                                trimmed_mean_prop_to_cut_polar=trimmed_mean_prop_to_cut_polar, 
                                                combination_method_intensity=combination_method_intensity, 
                                                trimmed_mean_prop_to_cut_intens=trimmed_mean_prop_to_cut_intens, 
                                                single_posang_north_up=single_posang_north_up)
    
    ###############################################################################
    # Subtract background in I_Q-, I_U-, Q- and U-frames
    ###############################################################################
 
    printandlog('\n###############################################################################')
    printandlog('# Subtracting the backgrounds and measuring and removing the star polarization')
    printandlog('###############################################################################') 

    # Determine background in corrected I_Q-, I_U-, Q- and U-frames and subtract it
    frame_I_Q_background_subtracted, background_frame_I_Q = subtract_background(cube=frame_I_Q_incident, annulus_background=annulus_background)
    frame_I_U_background_subtracted, background_frame_I_U = subtract_background(cube=frame_I_U_incident, annulus_background=annulus_background)
    frame_Q_background_subtracted, background_frame_Q = subtract_background(cube=frame_Q_incident, annulus_background=annulus_background)
    frame_U_background_subtracted, background_frame_U = subtract_background(cube=frame_U_incident, annulus_background=annulus_background)
   
    # Print resulting background values
    printandlog('\nSubtracted backgrounds in the incident I_Q-, I_U-, Q- and U-images:')
    printandlog('Background I_Q = %.3f' % background_frame_I_Q)    
    printandlog('Background I_U = %.3f' % background_frame_I_U)
    printandlog('Background Q = %.3f' % background_frame_Q)
    printandlog('Background U = %.3f' % background_frame_U)
                 
    ###############################################################################
    # Determine star polarization in frames and subtract it 
    ###############################################################################    
    
    # Compute normalized Stokes q and u, DoLP and AoLP in an annulus on the star
    q_star, u_star = determine_star_polarization(cube_I_Q=frame_I_Q_background_subtracted, 
                                                 cube_I_U=frame_I_U_background_subtracted, 
                                                 cube_Q=frame_Q_background_subtracted, 
                                                 cube_U=frame_U_background_subtracted, 
                                                 annulus_star=annulus_star, 
                                                 annulus_background=annulus_background)
    DoLP_star = np.sqrt(q_star**2 + u_star**2)
    AoLP_star = np.mod(np.rad2deg(0.5 * np.arctan2(u_star, q_star)), 180)
    
    # Print resulting star polarization
    printandlog('\nMeasured star polarization in the background-subtracted I_Q-, I_U-, Q- and U-images:')
    printandlog('q_star =    %.4f %%' % (100*q_star))
    printandlog('u_star =    %.4f %%' % (100*u_star))
    printandlog('DoLP_star = %.4f %%' % (100*DoLP_star))
    printandlog('AoLP_star = %.2f deg' % AoLP_star)
    printandlog('\nThe measured star polarization should be very similar to that fitted above when correcting the instrumental polarization effects. The difference in q, u and DoLP should be < 0.1% or << 0.1% depending on the noise in images.')

    # Subtract star polarization
    frame_Q_star_polarization_subtracted = frame_Q_background_subtracted - q_star*frame_I_Q_background_subtracted
    frame_U_star_polarization_subtracted = frame_U_background_subtracted - u_star*frame_I_U_background_subtracted
    
    # Subtract very small residual background
    frame_Q_star_polarization_subtracted, background_frame_Q_star_polarization_subtracted = subtract_background(cube=frame_Q_star_polarization_subtracted, 
                                                                                                                annulus_background=annulus_background)
    frame_U_star_polarization_subtracted, background_frame_U_star_polarization_subtracted = subtract_background(cube=frame_U_star_polarization_subtracted, 
                                                                                                                annulus_background=annulus_background)
    
    # Print resulting background values
    printandlog('\nSubtracted residual backgrounds in the star-polarization-subtracted Q- and U-images:')
    printandlog('Background Q = %.3f' % background_frame_Q_star_polarization_subtracted)
    printandlog('Background U = %.3f' % background_frame_U_star_polarization_subtracted)

    if len(cube_I_Q_incident) > 1 and len(cube_I_U_incident) > 1:
        ###############################################################################
        # Subtract background in I_Q-, I_U-, Q- and U-cubes
        ###############################################################################
        
        # Determine background in corrected I_Q-, I_U-, Q- and U-frames and subtract it
        cube_I_Q_background_subtracted, background_cube_I_Q = subtract_background(cube=cube_I_Q_incident, annulus_background=annulus_background)
        cube_I_U_background_subtracted, background_cube_I_U = subtract_background(cube=cube_I_U_incident, annulus_background=annulus_background)
        cube_Q_background_subtracted, background_cube_Q = subtract_background(cube=cube_Q_incident, annulus_background=annulus_background)
        cube_U_background_subtracted, background_cube_U = subtract_background(cube=cube_U_incident, annulus_background=annulus_background)
       
        # Print resulting background values
        printandlog('\nSubtracted mean backgrounds in the incident I_Q-, I_U-, Q- and U-image cubes:')
        printandlog('Mean background cube I_Q = %.3f' % np.mean(background_cube_I_Q))
        printandlog('Mean background cube I_U = %.3f' % np.mean(background_cube_I_U))
        printandlog('Mean background cube Q = %.3f' % np.mean(background_cube_Q))
        printandlog('Mean background cube U = %.3f' % np.mean(background_cube_U))
        
        ###############################################################################
        # Determine star polarization in cubes and plot it as function of HWP cycle number
        ###############################################################################    
    
        # Compute normalized Stokes q and u, DoLP and AoLP in an annulus on the star as a function of HWP cycle number
        q_star_HWP_cycle, u_star_HWP_cycle = determine_star_polarization(cube_I_Q=cube_I_Q_background_subtracted, 
                                                                         cube_I_U=cube_I_U_background_subtracted, 
                                                                         cube_Q=cube_Q_background_subtracted, 
                                                                         cube_U=cube_U_background_subtracted, 
                                                                         annulus_star=annulus_star, 
                                                                         annulus_background=annulus_background)
        DoLP_star_HWP_cycle = np.sqrt(q_star_HWP_cycle**2 + u_star_HWP_cycle**2)
        AoLP_star_HWP_cycle = np.mod(np.rad2deg(0.5 * np.arctan2(u_star_HWP_cycle, q_star_HWP_cycle)), 180)
        
        # Compute spread of q_star, u_star, DoLP_star and AoLP_star
        sigma_q_star = np.std(q_star_HWP_cycle, ddof=1)
        sigma_u_star = np.std(u_star_HWP_cycle, ddof=1)
        sigma_DoLP_star = np.std(DoLP_star_HWP_cycle, ddof=1)
        sigma_AoLP_star = np.std(AoLP_star_HWP_cycle, ddof=1)
       
        # Print resulting spread of star polarization
        printandlog('\nMeasured spread (standard deviation) of the star polarization with HWP cycle number:')
        printandlog('sigma_q_star =    %.4f %%' % (100*sigma_q_star))
        printandlog('sigma_u_star =    %.4f %%' % (100*sigma_u_star))
        printandlog('sigma_DoLP_star = %.4f %%' % (100*sigma_DoLP_star))
        printandlog('sigma_AoLP_star = %.2f deg' % sigma_AoLP_star)
        
        # Plot q, u and DoLP from annulus as function of HWP cycle number
        plot_name_star_quDoLP = name_file_root + 'star_pol_quDoLP.png'
        printandlog('\nCreating plot ' + plot_name_star_quDoLP + ' showing the measured star polarization as a function of HWP cycle number.')          
        font_size = 10
        x_max = len(q_star_HWP_cycle) + 1
        plt.figure(figsize = (4.7, 3.0))
        plt.plot([0, x_max], [0, 0], '-k')     
        plt.plot(np.arange(1, x_max), 100*DoLP_star_HWP_cycle, 'ok', label = 'DoLP')
        plt.plot([0, x_max], [100*DoLP_star, 100*DoLP_star], '-k')
        plt.plot(np.arange(1, x_max), 100*q_star_HWP_cycle, 'ob', label = 'q')  
        plt.plot([0, x_max], [100*q_star, 100*q_star], '-b')
        plt.plot(np.arange(1, x_max), 100*u_star_HWP_cycle, 'or', label = 'u')   
        plt.plot([0, x_max], [100*u_star, 100*u_star], '-r')     
        ax = plt.gca()
        ax.set_xlabel(r'HWP cycle', fontsize = font_size)
        ax.tick_params(axis = 'x', labelsize = font_size)
        ax.set_xlim([0, x_max])
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        ax.set_ylabel(r'Meas. star polarization (%)', fontsize = font_size)
        ax.tick_params(axis = 'y', labelsize = font_size)  
        ax.grid()
        plt.legend(loc = 'best')
        plt.tight_layout()
        plt.savefig(os.path.join(path_reduced_dir, plot_name_star_quDoLP), dpi = 300, bbox_inches = 'tight')
        plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name_star_quDoLP), dpi = 300, bbox_inches = 'tight')
        plt.close()
        
        # Plot AoLP from annulus as function of HWP cycle number
        plot_name_star_AoLP = name_file_root + 'star_pol_AoLP.png'
        printandlog('\nCreating plot ' + plot_name_star_AoLP + ' showing the measured angle of linear polarization of the star as a function of HWP cycle number.')
        plt.figure(figsize = (4.7, 3.0))
        plt.plot([0, x_max],[0, 0], '-k')     
        plt.plot(np.arange(1, x_max), AoLP_star_HWP_cycle, 'ok', label = 'AoLP') 
        plt.plot([0, x_max], [AoLP_star, AoLP_star], '-k')
        ax = plt.gca()
        ax.set_xlabel(r'HWP cycle', fontsize = font_size)
        ax.tick_params(axis = 'x', labelsize = font_size)
        ax.set_xlim([0, x_max])
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
        ax.set_ylabel(r'Meas. star AoLP ($^\circ$)', fontsize = font_size)
        ax.tick_params(axis = 'y', labelsize = font_size)  
        ax.grid()
        plt.legend(loc = 'best')
        plt.tight_layout()
        plt.savefig(os.path.join(path_reduced_dir, plot_name_star_AoLP), dpi = 300, bbox_inches = 'tight')
        plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name_star_AoLP), dpi = 300, bbox_inches = 'tight')
        plt.close()   
        
        printandlog('\nHorizontal trends in the data points of the plots ' + plot_name_star_quDoLP + ' and ' + plot_name_star_AoLP + ' indicate that the instrumental polarization effects have been removed successfully. However, this is only true provided that:')
        printandlog('\n1) the observations are taken with a sufficiently large range of parallactic and altitude angles,')
        printandlog('2) the observations are taken with a sufficiently high signal-to-noise ratio, and,')
        printandlog('3) the annulus for the star is placed in a region where there is only starlight.')
        printandlog('\nA non-zero measured star polarization then indicates the star is truly polarized, which is often caused by the presence micron-sized particles in the line of sight. This star polarization can therefore indicate the presence of an unresolved (inner) circumstellar disk, starlight passing through a resolved (outer) part of a circumstellar disk or the presence of interstellar dust between the star and the Earth.')

    ###############################################################################
    # Compute final images
    ###############################################################################
    
    printandlog('\n###############################################################################')
    printandlog('# Computing the final images and writing them to FITS-files')
    printandlog('###############################################################################') 

    # Compute final images with the star polarization still present
    frame_I_tot, frame_Q_phi, frame_U_phi, frame_I_pol, frame_AoLP, frame_DoLP, frame_q, frame_u, frame_AoLP_norm, frame_DoLP_norm \
    = compute_final_images(frame_I_Q=frame_I_Q_background_subtracted, 
                           frame_I_U=frame_I_U_background_subtracted, 
                           frame_Q=frame_Q_background_subtracted, 
                           frame_U=frame_U_background_subtracted, 
                           header=header, 
                           single_posang_north_up=single_posang_north_up)

    # Compute final images with the star polarization subtracted   
    frame_Q_phi_star_polarization_subtracted, frame_U_phi_star_polarization_subtracted, frame_I_pol_star_polarization_subtracted, \
    frame_AoLP_star_polarization_subtracted, frame_DoLP_star_polarization_subtracted, frame_q_star_polarization_subtracted, \
    frame_u_star_polarization_subtracted, frame_AoLP_norm_star_polarization_subtracted, frame_DoLP_norm_star_polarization_subtracted \
    = compute_final_images(frame_I_Q=frame_I_Q_background_subtracted, 
                           frame_I_U=frame_I_U_background_subtracted, 
                           frame_Q=frame_Q_star_polarization_subtracted, 
                           frame_U=frame_U_star_polarization_subtracted, 
                           header=header, 
                           single_posang_north_up=single_posang_north_up)[1:]

    # Create frames that show annuli used to retrieve star and background signals   
    frame_annulus_star = compute_annulus_values(cube=frame_I_Q_background_subtracted, param=annulus_star)[1]
    frame_annulus_background = compute_annulus_values(cube=frame_I_Q_background_subtracted, param=annulus_background)[1]

##TODO: Convert images to mJansky/". Perhaps it is best to give this as extra output and retain the original ones in counts too. What do you think?
## If so, just give the new images new clear names and add them to the section # Write .fits-files just below here. Note that you need to convert
## the images with and without star polarization.
#    
#    ###############################################################################
#    # Optionally convert final Q-, U-, Qphi-, Uphi- and Ipol-images to mJansky/"
#    ###############################################################################
#
#    if bla == True:  
#        # Convert final Q-, U-, Qphi-, Uphi- and Ipol-images with star polarization to mJansky/"
#        dummy = 5
#        # Convert final Q-, U-, Qphi-, Uphi- and Ipol-images with star polarization to mJansky/"
        
    ###############################################################################
    # Print image orientation of final images
    ###############################################################################
    
    # Determine tracking mode used
    tracking_mode_used = header[0]['ESO INS4 COMB ROT']
    
    # Determine number of unique derotator position angles
    derotator_position_angle = [x['ESO INS4 DROT2 POSANG'] for x in header]
    number_derotator_position_angles = len(np.unique(derotator_position_angle))

    # Calculate rotation and determine whether it is clockwise or counterclockwise
    rotation_angle = -derotator_position_angle[0] - true_north_correction
    if np.sign(rotation_angle) == -1:
        rotation_direction = 'clockwise'
    elif np.sign(rotation_angle) == 1:
        rotation_direction = 'counterclockwise'
        
    # Print image orientation of final images
    if tracking_mode_used == 'FIELD' and number_derotator_position_angles == 1 and single_posang_north_up == False and np.sign(rotation_angle) != 0:
        printandlog('\nFinal images are rotated %.4f deg ' % (np.abs(rotation_angle)) + rotation_direction + ' with respect to North up.')
    else:
        printandlog('\nFinal images are oriented with North up.')
    
    ###############################################################################
    # Write .fits-files
    ###############################################################################
            
    # List files of the images with the star polarization present and define their file names
    frames_to_write = [frame_I_Q_background_subtracted, frame_I_U_background_subtracted, frame_I_tot, frame_Q_background_subtracted, 
                       frame_U_background_subtracted, frame_Q_phi, frame_U_phi, frame_I_pol, frame_AoLP]
    file_names = ['I_Q', 'I_U', 'I_tot', 'Q', 'U', 'Q_phi', 'U_phi', 'I_pol', 'AoLP']
    
    if normalized_polarization_images == True:
        # Add images of DoLP, normalized Stokes q and u and AoLP and DoLP created using q- and u-images
        printandlog('\nWARNING, the images DoLP.fits, q_norm.fits, u_norm.fits, AoLP_norm.fits and DoLP_norm.fits are only valid if all flux in the images originates from the astrophysical source of interest. This is generally the case for observations of for example solar system objects or galaxies. The images are generally not valid for observations of circumstellar disks or companions because in that case a large part of the flux in the total intensity images originates from the central star.')
        frames_to_write += [frame_DoLP, frame_q, frame_u, frame_AoLP_norm, frame_DoLP_norm]       
        file_names += ['DoLP', 'q_norm', 'u_norm', 'AoLP_norm', 'DoLP_norm']
    
    # Write files of the images with the star polarization present
    printandlog('')
    for frame, file_name in zip(frames_to_write, file_names):
        write_fits_files(data=frame, path=os.path.join(path_reduced_dir, name_file_root + file_name + '.fits'), header=False)

    # Convert in contrast per arcsec^2 and write the files of the images with the star polarization present to contrast per arcsec 
    if convert_in_contrast_per_arcsec2 == True:      
        for frame, file_name in zip(frames_to_write, file_names):
            write_fits_files(data=frame/reference_flux/pixel_scale**2,\
                             path=os.path.join(path_reduced_dir, \
                             name_file_root + file_name + '_contrast_per_arcsec2.fits'), header=False)

    # Write frames that show annuli used to retrieve star and background signals in reduced directory
    write_fits_files(data=frame_annulus_star, path=os.path.join(path_reduced_dir, name_file_root + 'annulus_star.fits'), header=False)
    write_fits_files(data=frame_annulus_background, path=os.path.join(path_reduced_dir, name_file_root + 'annulus_background.fits'), header=False)    

    # List files of the images with the star polarization subtracted and define their file names
    frames_to_write = [frame_I_Q_background_subtracted, frame_I_U_background_subtracted, frame_I_tot, frame_Q_star_polarization_subtracted, frame_U_star_polarization_subtracted, \
                       frame_Q_phi_star_polarization_subtracted, frame_U_phi_star_polarization_subtracted, frame_I_pol_star_polarization_subtracted, frame_AoLP_star_polarization_subtracted]
    file_names = ['I_Q', 'I_U', 'I_tot', 'Q_star_pol_subtr', 'U_star_pol_subtr', 'Q_phi_star_pol_subtr', 'U_phi_star_pol_subtr', 'I_pol_star_pol_subtr', 'AoLP_star_pol_subtr']

    if normalized_polarization_images == True:
        # Add images of DoLP, normalized Stokes q and u and AoLP and DoLP created using q- and u-images
        frames_to_write += [frame_DoLP_star_polarization_subtracted, frame_q_star_polarization_subtracted, frame_u_star_polarization_subtracted, \
                            frame_AoLP_norm_star_polarization_subtracted, frame_DoLP_norm_star_polarization_subtracted]       
        file_names += ['DoLP_star_pol_subtr', 'q_norm_star_pol_subtr', 'u_norm_star_pol_subtr', 'AoLP_norm_star_pol_subtr', 'DoLP_norm_star_pol_subtr']

    # Write files of the images with the star polarization subtracted
    for frame, file_name in zip(frames_to_write, file_names):
        write_fits_files(data=frame, path=os.path.join(path_reduced_star_pol_subtr_dir, name_file_root + file_name + '.fits'), header=False)

    # Write frames that show annuli used to retrieve star and background signals in reduced_star_pol_subtr directory
    write_fits_files(data=frame_annulus_star, path=os.path.join(path_reduced_star_pol_subtr_dir, name_file_root + 'annulus_star.fits'), header=False)
    write_fits_files(data=frame_annulus_background, path=os.path.join(path_reduced_star_pol_subtr_dir, name_file_root + 'annulus_background.fits'), header=False)    


    ########################################################################################
    # Optionally convert final Q-, U-, Qphi-, Uphi- and Ipol-images to contrast per arcsec^2
    ########################################################################################

    if convert_in_contrast_per_arcsec2 == True:  
        # Convert final Q-, U-, Qphi-, Uphi- and Ipol-images with star polarization subtracted to contrast per arcsec    
        for frame, file_name in zip(frames_to_write, file_names):
            write_fits_files(data=frame/reference_flux/pixel_scale**2, \
                    path=os.path.join(path_reduced_star_pol_subtr_dir, \
                    name_file_root + file_name + '_contrast_per_arcsec2.fits'), header=False)

    # Write frames that show annuli used to retrieve star and background signals in reduced_star_pol_subtr directory
    write_fits_files(data=frame_annulus_star, path=os.path.join(path_reduced_star_pol_subtr_dir, name_file_root + 'annulus_star.fits'), header=False)
    write_fits_files(data=frame_annulus_background, path=os.path.join(path_reduced_star_pol_subtr_dir, name_file_root + 'annulus_background.fits'), header=False)    





## TODO: I suggest you write the main function to perform ADI here (with sub functions above), and that gets
## the input cube_left_frames, cube_right_frames all the way at the bottom of the script, below the call
## of perform_postprocessing. That way the ADI part shares the pre-processed data, but is totally 
## separate from the polarimetry part. Should make the implementation much easier.
#
################################################################################
## perform_adi
################################################################################
#
#def perform_adi(cube_left_frames, cube_right_frames):
#    '''
#    ...
#    '''
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           
                           


###############################################################################
# run_demo
###############################################################################    

def run_demo(path_main_dir):
    '''
    Run example reduction with data from T Cha as published in Pohl et al.
    (2017) and  and used in van Holstein et al. (2019)
    
    Input:
        path_main_dir: string specifying path to main directory
    
    File written by Rob van Holstein
    Function status: verified    
    '''

    print_wrap('\nStarting example reduction using data of the circumstellar ' +
               'disk of T Cha (1 HWP cycle) as published in Pohl et al. (2017) ' +
               'and used in van Holstein et al. (2019).')
    
    # Define path of raw directory
    path_raw_dir = os.path.join(path_main_dir, 'raw')
  
    # Define names of example data
    name_example_files = ['SPHER.2016-02-20T06_59_15.857.fits',
                          'SPHER.2016-02-20T07_00_44.006.fits', 
                          'SPHER.2016-02-20T07_02_09.295.fits',
                          'SPHER.2016-02-20T07_02_49.297.fits',
                          'SPHER.2016-02-20T07_03_27.243.fits',
                          'SPHER.2016-02-20T07_04_07.281.fits',
                          'SPHER.2016-02-20T08_25_14.427.fits']
    
    # Download example data to directory with example data
    data_exists = all([os.path.exists(os.path.join(path_raw_dir, file)) for file in name_example_files])
    if not data_exists:
        # Ask if user wants to download the data
        proceed_download = input_wrap('\nTo run the demo, 56.6 MB of data needs to be downloaded. Proceed? (y/n) ')

        if proceed_download == 'y':
            if not os.path.exists(path_raw_dir):
                # Create raw directory
                os.makedirs(path_raw_dir)
                print_wrap('\nCreated raw directory {0:s}'.format(path_raw_dir))

            # Define URL of example data
            url_example_data = 'https://github.com/robvanholstein/IRDAP/raw/master/irdap/example_data/'
    
            # Download example data
            print_wrap('\nDownloading data:\n')
            for i, file in enumerate(name_example_files):
                urllib.request.urlretrieve(url_example_data + file, os.path.join(path_raw_dir, file))
                print_wrap('Downloaded file ' + str(i + 1) + '/' + str(len(name_example_files)) + ': ' + file)    
        elif proceed_download == 'n':
            print_wrap('\nNot downloading data, aborting.')
        else:
            print_wrap('\nThe provided input \'' + str(proceed_download) + '\' is not valid.') 
    else:
        print_wrap('\nExample data has been downloaded already.')

    if data_exists or proceed_download == 'y':
        # Create a configuration file in the main directory
        make_config(path_main_dir)
        
        # Run the pipeline
        run_pipeline(path_main_dir)
            
###############################################################################
# make_config
###############################################################################    

def make_config(path_main_dir):
    '''
    Create the default configuration file in the main directory
    
    Input:
        path_main_dir: string specifying path to main directory
    
    File written by Rob van Holstein
    Function status: verified    
    '''
    
    # Define paths of default configuration file and configuration file to be written
    path_default_config_file = os.path.join(os.path.dirname(__file__), 'config.conf')
    path_config_file_write = os.path.join(path_main_dir, os.path.basename(path_default_config_file))
    
    if os.path.exists(path_config_file_write):
        # Ask if user wants to overwrite the configuration file
        overwrite = input_wrap('\nThe configuration file ' + path_config_file_write + 
                               ' already exists. Do you want to overwrite it? (y/n) ')
        if overwrite == 'n':
            print_wrap('\nNo configuration file has been created.')
        elif overwrite != 'y':
            print_wrap('\nThe provided input \'' + str(overwrite) + '\' is not valid.')
    
    if not os.path.exists(path_config_file_write) or overwrite == 'y':
        # Copy default configuration file to main directory
        shutil.copyfile(path_default_config_file, path_config_file_write)
        print_wrap('\nCreated a default configuration file ' + path_config_file_write + '.')

###############################################################################
# mean_combine_images
###############################################################################  

def mean_combine_images(path_main_dir, path_read_dirs):
    '''
    Mean-combine the images of two or more reductions
    
    Input:
        path_main_dir: string specifying path to main directory
        path_read_dirs: list of strings specifying paths to directories to read data from
        
    File written by Rob van Holstein
    Function status: verified    
    '''

    for reduced_dir in ['reduced', 'reduced_star_pol_subtr']:    
        # Define path of directory to write combined images to
        path_write_dir_sel = os.path.join(path_main_dir, reduced_dir + '_combined')
        
        if not os.path.exists(path_write_dir_sel):
            # Create directory if it does not exist yet
            os.makedirs(path_write_dir_sel)
        
        # Create empty lists to store path of files in and whether or not normalized polarization images are present      
        paths_I_Q = []    
        paths_I_U = []    
        paths_Q = []    
        paths_U = []    
        paths_Q_phi = []    
        paths_U_phi = []    
        normalized_polarization_images = []  
        
        for path_sel in path_read_dirs:
            # Define path of directory containing files to be read and check if it exists
            path_reduced_dir_sel = os.path.join(path_sel, reduced_dir)
            
            if not os.path.exists(path_reduced_dir_sel):
                raise IOError('\n\nThe directory ' + path_reduced_dir_sel + ' does not exist.')
            
            # Retrieve paths of files to be read
            path_I_Q = glob.glob(os.path.join(path_reduced_dir_sel, '*_I_Q*.fits'))
            if len(path_I_Q) != 1:
                raise IOError('\n\nThere are no or multiple I_Q-images in the directory ' + path_reduced_dir_sel + '.')            
            paths_I_Q.append(path_I_Q[0])
            
            path_I_U = glob.glob(os.path.join(path_reduced_dir_sel, '*_I_U*.fits'))
            if len(path_I_U) != 1:
                raise IOError('\n\nThere are no or multiple I_U-images in the directory ' + path_reduced_dir_sel + '.')            
            paths_I_U.append(path_I_U[0])
    
            path_Q = [x for x in glob.glob(os.path.join(path_reduced_dir_sel, '*_Q*.fits')) if '_I_Q' not in x and '_Q_phi' not in x and '_q_norm' not in x]
            if len(path_Q) != 1:
                raise IOError('\n\nThere are no or multiple Q-images in the directory ' + path_reduced_dir_sel + '.')            
            paths_Q.append(path_Q[0])
            
            path_U = [x for x in glob.glob(os.path.join(path_reduced_dir_sel, '*_U*.fits')) if '_I_U' not in x and '_U_phi' not in x and '_u_norm' not in x]
            if len(path_U) != 1:
                raise IOError('\n\nThere are no or multiple U-images in the directory ' + path_reduced_dir_sel + '.')            
            paths_U.append(path_U[0])
            
            path_Q_phi = glob.glob(os.path.join(path_reduced_dir_sel, '*_Q_phi*.fits'))
            if len(path_Q_phi) != 1:
                raise IOError('\n\nThere are no or multiple Q_phi-images in the directory ' + path_reduced_dir_sel + '.')            
            paths_Q_phi.append(path_Q_phi[0])
            
            path_U_phi = glob.glob(os.path.join(path_reduced_dir_sel, '*_U_phi*.fits'))
            if len(path_U_phi) != 1:
                raise IOError('\n\nThere are no or multiple U_phi-images in the directory ' + path_reduced_dir_sel + '.')            
            paths_U_phi.append(path_U_phi[0])       
    
            # Check if normalized polarization images are present
            path_list =" ".join(glob.glob(os.path.join(os.path.join(path_sel, reduced_dir), '*.fits')))  
            normalized_polarization_images.append(all([x in path_list for x in ['DoLP', 'q_norm', 'u_norm', 'AoLP_norm', 'DoLP_norm']]))   
        
        if reduced_dir == 'reduced':
            # Define name_file_root based on names of files in each read directory
            name_file_roots = [os.path.basename(x)[:os.path.basename(x).index('_I_Q')] for x in paths_I_Q]
            target_name = [x[:-10] for x in name_file_roots][0]
            dates = [x[-10:] for x in name_file_roots]
            dates_sorted = sorted(dates, key=lambda d: tuple(map(int, d.split('-'))))
            date_first = dates_sorted[0]
            date_last = dates_sorted[-1]
            if date_first[:4] != date_last[:4]:
                date_root = date_first + '_' + date_last
            elif date_first[5:7] != date_last[5:7]:
                date_root = date_first[:5] + date_first[5:] + '_' + date_last[5:]
            elif date_first[8:] != date_last[8:]:
                date_root = date_first[:8] + date_first[8:] + '_' + date_last[8:]
            else:
                date_root = date_first
            name_file_root = target_name + date_root + '_'
        
        # Read files and mean-combine them
        frame_I_Q = np.mean(np.vstack(read_fits_files(paths_I_Q, silent=True)[0]), axis=0)      
        frame_I_U = np.mean(np.vstack(read_fits_files(paths_I_U, silent=True)[0]), axis=0)      
        frame_Q = np.mean(np.vstack(read_fits_files(paths_Q, silent=True)[0]), axis=0)      
        frame_U = np.mean(np.vstack(read_fits_files(paths_U, silent=True)[0]), axis=0)      
        frame_Q_phi = np.mean(np.vstack(read_fits_files(paths_Q_phi, silent=True)[0]), axis=0)      
        frame_U_phi = np.mean(np.vstack(read_fits_files(paths_U_phi, silent=True)[0]), axis=0)      
            
        # Compute final images
        frame_I_tot, _, _, frame_I_pol, frame_AoLP, frame_DoLP, frame_q, frame_u, frame_AoLP_norm, frame_DoLP_norm \
        = compute_final_images(frame_I_Q=frame_I_Q, 
                               frame_I_U=frame_I_U, 
                               frame_Q=frame_Q, 
                               frame_U=frame_U, 
                               header=None, 
                               single_posang_north_up=True)
       
        # List files of the images and define their (base) file names
        frames_to_write = [frame_I_Q, frame_I_U, frame_I_tot, frame_Q, frame_U, 
                           frame_Q_phi, frame_U_phi, frame_I_pol, frame_AoLP]
        file_names = ['I_Q', 'I_U', 'I_tot', 'Q', 'U', 'Q_phi', 'U_phi', 'I_pol', 'AoLP']
            
        if all(normalized_polarization_images):
            # Add images of DoLP, normalized Stokes q and u and AoLP and DoLP created using q- and u-images
            frames_to_write += [frame_DoLP, frame_q, frame_u, frame_AoLP_norm, frame_DoLP_norm]       
            file_names += ['DoLP', 'q_norm', 'u_norm', 'AoLP_norm', 'DoLP_norm']
        
        # Add substring '_star_pol_subtr' to appropriate files
        if reduced_dir == 'reduced_star_pol_subtr':
            for i in range(len(file_names)):
                if file_names[i] not in ['I_Q', 'I_U', 'I_tot']:
                    file_names[i] += '_star_pol_subtr'

        # Write files of the combined images 
        for frame, file_name in zip(frames_to_write, file_names):
            write_fits_files(data=frame, path=os.path.join(path_write_dir_sel, name_file_root + file_name + '.fits'), header=False, silent=True)
    
    # Write TXT-file with directories from which the original images come        
    f = open(os.path.join(path_main_dir, 'images_mean_combined.txt'), 'w+')
    f.write('The combined images are computed as the mean of the corresponding images in the following directories:')
    for path_sel in path_read_dirs:
        f.write('\n' + path_sel)
    f.close()
    
    print('\nSuccessfully mean-combined the images located in:')
    for path_sel in path_read_dirs: 
        print(path_sel)
       
###############################################################################
# create_overview_headers_main
###############################################################################    
    
def create_overview_headers_main(path_main_dir):
    '''
    Create an overview of relevant FITS-headers and write it to a text-file
    from the __main__.py file
    
    Input:
        path_main_dir: string specifying path to main directory
        
    File written by Rob van Holstein
    Function status: verified
    '''

    # Define path of raw directory
    path_raw_dir = os.path.join(path_main_dir, 'raw')

    # Check if raw directory exists, if not create it
    if not os.path.exists(path_raw_dir):
        os.makedirs(path_raw_dir)
        print_wrap('\nThe raw directory {0:s} did not exist. It was created but you need to put your raw FITS-files there.'.format(path_raw_dir))
    else:
        # List FITS-files in raw directory
        path_raw_files = glob.glob(os.path.join(path_raw_dir,'*.fits'))
        
        # Check if raw folder contains FITS-files
        if len(path_raw_files) == 0:
            print_wrap('\nThe raw directory {0:s} does not contain FITS-files. You need to put your raw FITS-files in this folder.'.format(path_raw_dir))
        else:
            # Define the base of the name of each file to be generated
            header_first_file = pyfits.getheader(path_raw_files[0])
            target_name = header_first_file['ESO OBS TARG NAME']
            date_obs = header_first_file['DATE-OBS']
            
            name_file_root = target_name.replace(' ', '_') + '_' + date_obs[:10].replace(' ', '_') + '_'
        
            # Define path to header overview to be created
            path_overview = os.path.join(path_main_dir, name_file_root + 'headers.txt')
          
            # Create overview of headers
            create_overview_headers(path_raw_dir, path_overview, log=False)
            print_wrap('\nCreated overview of the headers ' + path_overview + ' and ' + \
                       path_overview.replace('.txt','.csv') + '.')
            
###############################################################################
# run_pipeline
###############################################################################

def run_pipeline(path_main_dir):
    '''  
    Run the pipeline
    
    Input:
        path_main_dir: string specifying path of main directory containing the 
            configuration file 'config.conf' and the directory 'raw' with the
            raw data
       
    File written by Rob van Holstein
    Function status: verified
    ''' 

    # Start taking time
    time_start = time.time()
    
    # Turn off pyfits FITS-warnings
    warnings.filterwarnings('ignore', category=UserWarning)

    # Print that data-reduction starts    
    print_wrap('\nStarting data reduction.')
        
    ###############################################################################
    # Define global variables
    ###############################################################################
    
    # Define which variables should be global
    global pupil_offset
    global true_north_correction
    global pixel_scale
    global msd
    global path_raw_dir
    global path_calib_dir
    global path_flat_dir
    global path_bpm_dir
    global path_sky_dir
    global path_center_dir
    global path_flux_dir
    global path_sky_flux_dir
    global path_preprocessed_dir
    global path_reduced_dir
    global path_reduced_star_pol_subtr_dir
    global name_file_root
    global path_log_file
    global path_overview
    global path_static_calib_dir
    global convert_in_contrast_per_arcsec2
    
    # Define pupil-offset (deg) in pupil-tracking mode (SPHERE User Manual P99.0, 6th public release, P99 Phase 1)
    pupil_offset = 135.99
        
    # Define true North correction (deg) in pupil-tracking mode (SPHERE User Manual P99.0, 6th public release, P99 Phase 1)
    true_north_correction = -1.75

    # Define pixel scale in arcsec/px from Maire et al. 2016 (average value valid between all filters)
    pixel_scale = 0.01225
        
    # Define mean solar day (s)
    msd = 86400
    
    # Define paths of directories
    path_raw_dir = os.path.join(path_main_dir, 'raw')
    path_log_config_dir = os.path.join(path_main_dir, 'logs')
    path_calib_dir = os.path.join(path_main_dir, 'calibration')
    path_flat_dir = os.path.join(path_calib_dir, 'flat')
    path_bpm_dir = os.path.join(path_calib_dir, 'bpm')
    path_sky_dir = os.path.join(path_calib_dir, 'sky')
    path_center_dir = os.path.join(path_calib_dir, 'center')
    path_flux_dir = os.path.join(path_calib_dir, 'flux')
    path_sky_flux_dir = os.path.join(path_calib_dir, 'sky_flux')
    path_preprocessed_dir = os.path.join(path_main_dir, 'preprocessed')
    path_reduced_dir = os.path.join(path_main_dir, 'reduced')
    path_reduced_star_pol_subtr_dir = os.path.join(path_main_dir, 'reduced_star_pol_subtr')
        
    # Check if raw directory exists, if not create it
    if not os.path.exists(path_raw_dir):
        os.makedirs(path_raw_dir)
        raise IOError('\n\nThe raw directory {0:s} did not exist. It was created but you need to put your raw FITS-files there.'.format(path_raw_dir))
    
    # List FITS-files in raw directory
    path_raw_files = glob.glob(os.path.join(path_raw_dir,'*.fits'))
    
    # Check if raw folder contains FITS-files
    if len(path_raw_files) == 0:
        raise IOError('\n\nThe raw directory {0:s} does not contain FITS-files. You need to put your raw FITS-files in this folder.'.format(path_raw_dir))

    # Define the base of the name of each file to be generated
    header_target_name = [pyfits.getheader(x) for x in path_raw_files if 'ESO OBS TARG NAME' in pyfits.getheader(x)]
    if len(header_target_name) == 0:
        raise IOError('\n\nThe target name has not been found in the headers. There is likely no on-sky data in the raw directory {0:s}.'.format(path_raw_dir))
    else:
        target_name = header_target_name[0]['ESO OBS TARG NAME']
        date_obs = header_target_name[0]['DATE-OBS']
    
    name_file_root = target_name.replace(' ', '_') + '_' + date_obs[:10].replace(' ', '_') + '_'

    # Find path of log file from previous reduction
    path_log_file_old = glob.glob(os.path.join(path_main_dir,'*_log_*'))

    if len(path_log_file_old) > 1:
        raise IOError('\n\nThere should only be one log file in the directory ' + path_main_dir + '. Please remove the latest one.')        
    elif len(path_log_file_old) == 1:
        # Extract path of old log file and its number
        path_log_file_old = path_log_file_old[0]
        log_file_old_number = int(os.path.splitext(path_log_file_old)[0].split('log_')[1])
    else:
        # Set file number to zero
        log_file_old_number = 0

    # Define paths of log file, header overview and directory containing static calibrations
    path_log_file = os.path.join(path_main_dir, name_file_root + 'log_' + str(log_file_old_number + 1) + '.txt')
    path_overview = os.path.join(path_main_dir, name_file_root + 'headers.txt')
    path_static_calib_dir = os.path.join(os.path.dirname(__file__), 'static_calibs')
     
    ###############################################################################
    # Check if configuration and log file exist
    ###############################################################################   
    
    # Define the path of the configuration file and check if it exists
    path_config_file = os.path.join(path_main_dir, 'config.conf')
    if not os.path.exists(path_config_file):
        raise IOError('\n\nThere is no configuration file ' + path_config_file + '. Run \'irdap --makeconfig\' first.')
       
    # Find path to copy of configuration file from previous reduction
    path_config_file_copy_old = glob.glob(os.path.join(path_main_dir,'*_config_*'))    
    
    if len(path_config_file_copy_old) > 1:
        raise IOError('\n\nThere should only be one copy of the configuration file in the directory ' + path_main_dir + '. Please remove the latest one.') 
    elif len(path_config_file_copy_old) == 1:
        # Extract path of old copy of configuration file and its number
        path_config_file_copy_old = path_config_file_copy_old[0]
        config_file_copy_old_number = int(os.path.splitext(path_config_file_copy_old)[0].split('config_')[1])
    else:
        # Set file number to zero
        config_file_copy_old_number = 0
        
    # Raise error if there is a log file of the previous reduction but no copy of the configuration file or vice versa
    if len(path_log_file_old) != 0 and len(path_config_file_copy_old) == 0:
        raise IOError('\n\nThere is a log file of the previous reduction, but no copy of the configuration file. Please remove the log file.')
    
    if len(path_log_file_old) == 0 and len(path_config_file_copy_old) != 0:
        raise IOError('\n\nThere is a copy of the configuration file of the previous reduction, but no log file. Please remove the copy of the configuration file.')

    ###############################################################################
    # Check if there is a configuration file and start writing log file
    ###############################################################################
            
    # Create boolean saying whether log file already existed
    log_file_existed = len(path_log_file_old) != 0
    
    if log_file_existed == True:       
        # Save relevant lines from pre-processing
        log_file_lines = [x.rstrip('\n') for x in open(path_log_file_old, 'r')]
        if '# Starting post-processing' in log_file_lines:
            log_file_lines = log_file_lines[log_file_lines.index('# Starting pre-processing') - 2: \
                                            log_file_lines.index('# Starting post-processing') - 2]
        else:
            log_file_lines = None
               
    # Start writing log file
    time_now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    printandlog('\n###############################################################################')
    printandlog('# Important notice')
    printandlog('###############################################################################')
    printandlog('\nWhen publishing data reduced with IRDAP, please cite van Holstein et al.')
    printandlog('(2019): <ADS link>.')
    printandlog('For data in pupil-tracking mode please additionally cite van Holstein et al.')
    printandlog('(2017): http://adsabs.harvard.edu/abs/2017SPIE10400E..15V.')
    printandlog('\nFull documentation: https://robvanholstein.github.io/IRDAP')
    printandlog('Feedback, questions, comments: vanholstein@strw.leidenuniv.nl')
    printandlog('\nIRDAP Copyright (C) 2019 R.G. van Holstein')
    printandlog('\n###############################################################################')
    printandlog('# Starting IRDAP')
    printandlog('###############################################################################')
    printandlog('\nRunning IRDAP version ' + __version__)
    printandlog('\nDate: ' + time_now[:time_now.find(' ')])
    printandlog('Time: ' + time_now[time_now.find(' ') + 1:])   
    printandlog('\n###############################################################################')
    printandlog('# Preparing data reduction')
    printandlog('###############################################################################') 
    printandlog('\nCreated log file ' + path_log_file + '.')

    ###############################################################################
    # Retrieve and define input parameters
    ###############################################################################
    
    # Read configuration file
    printandlog('\nReading configuration file ' + path_config_file + '.')
    
    sigma_filtering, \
    object_collapse_ndit, \
    object_centering_method, \
    skip_preprocessing, \
    frames_to_remove, \
    annulus_star, \
    annulus_background, \
    normalized_polarization_images, \
    center_subtract_object, \
    center_param_centering, \
    object_center_coordinates, \
    object_param_centering, \
    flux_centering_method, \
    flux_center_coordinates, \
    flux_param_centering, \
    flux_annulus_background, \
    flux_annulus_star, \
    double_difference_type, \
    single_posang_north_up, \
    convert_in_contrast_per_arcsec2 \
    = read_config_file(path_config_file)
    
    # Define some fixed input parameters
    save_preprocessed_data = True
    show_images_center_coordinates = True
    remove_vertical_band_detector_artefact = True
    combination_method_polarization = 'least squares'
    combination_method_intensity = 'mean'
    trimmed_mean_prop_to_cut_polar = 0.1
    trimmed_mean_prop_to_cut_intens = 0.1
    
    # Check validity of input of skip_preprocessing
    if skip_preprocessing not in [True, False]:
        raise ValueError('\n\n\'skip_preprocessing\' should be either True or False. Before starting another reduction, please delete the log file ' + path_log_file + '.')   

    if skip_preprocessing == True:
        # Raise errors if there is no log file from the previous reduction or pre-processing was not finished
        if log_file_existed is False:
            raise IOError('\n\nThere is no log file from a previous reduction. Remove the log file that has been created for the current reduction attempt and then run IRDAP with \'skip_preprocessing\' = False to perform the pre-processing of the raw data and save the results.')
        elif log_file_lines is None:
            raise IOError('\n\nThe pre-processing part of the reduction is not complete in the previous log file. Remove the log file that has been created last and then run IRDAP with \'skip_preprocessing\' = False to perform the pre-processing of the raw data and save the results.')
        
    ###############################################################################
    # Make a copy of configuration file
    ###############################################################################
    
    # Read lines of configuration file
    config_file_lines = [x for x in open(path_config_file, 'r')]
    
    # Define path of new copy of configuration file
    path_config_file_copy_new = os.path.join(path_main_dir, name_file_root + \
                                             os.path.splitext(os.path.basename(path_config_file))[0] + \
                                             '_' + str(config_file_copy_old_number + 1) + os.path.splitext(path_config_file)[1])

    if skip_preprocessing == True:
        if len(path_config_file_copy_old) == 0:
            raise IOError('\n\nThere was no copy of the configuration file from a previous reduction. Run IRDAP first with \'skip_preprocessing\' = False to perform the pre-processing of the raw data and save the results.')

        # Read lines of existing copy of configuration file
        config_file_lines_old = [x for x in open(path_config_file_copy_old, 'r')]
        
        # Define indices of lines pertaining to pre-processing input parameters
        n0 = config_file_lines.index('[Basic pre-processing options]\n') + 2
        n1 = config_file_lines.index('[Basic post-processing options]\n') - 1
        n2 = config_file_lines.index('[Advanced pre-processing options]\n') + 3
        n3 = config_file_lines.index('[Advanced post-processing options]\n') - 1 
        
        # Replace pre-processing lines of configuration file by those from the previous copy 
        config_file_lines[n0:n1] = config_file_lines_old[n0:n1]
        config_file_lines[n2:n3] = config_file_lines_old[n2:n3]

    # Write lines in copy of configuration file and save to main directory
    with open(path_config_file_copy_new, 'w') as f:
        for x in config_file_lines:
            f.write(x)
    printandlog('\nCreated a copy of the used configuration file ' + path_config_file_copy_new + '.')
    if skip_preprocessing == True:
        printandlog('\nBecause the pre-processing is skipped, the pre-processing input parameters from the previous reduction have been used for this copy.')

    ###############################################################################
    # Copy previous log and configuration files to separate directory
    ###############################################################################
   
    if log_file_existed and len(path_config_file_copy_old) != 0:
        log_config_dir_existed = os.path.exists(path_log_config_dir)
        if not log_config_dir_existed:
            # Create directory to save log file and copy of configuration file of previous reduction to
            os.makedirs(path_log_config_dir)
    
        # Move log file and copy of configuration file to directory
        shutil.move(path_log_file_old, os.path.join(path_log_config_dir, os.path.basename(path_log_file_old)))
        shutil.move(path_config_file_copy_old, os.path.join(path_log_config_dir, os.path.basename(path_config_file_copy_old)))
    
        if not log_config_dir_existed:
            printandlog('\nCreated directory ' + path_log_config_dir + ' and moved log file and copy of configuration file of previous reduction there.')
        else:
            printandlog('\nMoved log file and copy of configuration file of previous reduction to directory ' + path_log_config_dir + '.')
        
    ###############################################################################
    # Check whether input values are valid (note that checks do not account for all possibilities)
    ###############################################################################
        
    printandlog('\nChecking input parameters of configuration file.')
    
    # frames_to_remove
    if type(frames_to_remove) is not list:
          raise TypeError('\n\n\'frames_to_remove\' should be an empty list or a list of integers and/or length-2-tuples of integers.')
      
    if len(frames_to_remove) != 0:
        if any([type(x) not in [int, tuple] for x in frames_to_remove]):
            raise TypeError('\n\n\'frames_to_remove\' should be an empty list or a list of integers and/or length-2-tuples of integers.')
        if any([len(x) != 2 for x in frames_to_remove if type(x) is tuple]):
            raise TypeError('\n\n\'frames_to_remove\' should be an empty list or a list of integers and/or length-2-tuples of integers.')
        if any([type(y) is not int for x in frames_to_remove if type(x) is tuple for y in x]):
            raise TypeError('\n\n\'frames_to_remove\' should be an empty list or a list of integers and/or length-2-tuples of integers.')
        if any([x < 1 for x in frames_to_remove if type(x) is int]):
            raise ValueError('\n\nThe file indices in \'frames_to_remove\' are 1-based and therefore should be larger than zero.')    
        if any([y < 1 for x in frames_to_remove if type(x) is tuple for y in x]):
            raise ValueError('\n\nThe file indices in \'frames_to_remove\' are 1-based and therefore should be larger than zero.')    
    
    # sigma_filtering 
    if sigma_filtering not in [True, False]:
        raise ValueError('\n\n\'sigma_filtering\' should be either True or False.')   
    
    # object_collapse_ndit
    if object_collapse_ndit not in [True, False]:
        raise ValueError('\n\n\'object_collapse_ndit\' should be either True or False.')   
                              
    # show_images_center_coordinates
    if show_images_center_coordinates not in [True, False]:
        raise ValueError('\n\n\'show_images_center_coordinates\' should be either True or False.')   
    
    # object_centering_method
    if object_centering_method not in ['automatic', 'center frames', 'gaussian', 'cross-correlation', 'manual']:
        raise ValueError('\n\n\'object_centering_method\' should be \'automatic\', \'center frames\', \'gaussian\', \'cross-correlation\' or \'manual\'.')
    
    # center_subtract_object
    if center_subtract_object not in [True, False]:
        raise ValueError('\n\n\'center_subtract_object\' should be either True or False.')   
       
    # object_center_coordinates
    if type(object_center_coordinates) not in [str, tuple]:
        raise TypeError('\n\n\'object_center_coordinates\' should be \'automatic\' or a tuple of length 4 containing floats or integers.')
        
    if object_center_coordinates == ():
        raise TypeError('\n\n\'object_center_coordinates\' should be \'automatic\' or a tuple of length 4 containing floats or integers.')
        
    if type(object_center_coordinates) is str and object_center_coordinates != 'automatic':
        raise ValueError('\n\n\'object_center_coordinates\' should be \'automatic\' or a tuple of length 4 containing floats or integers.')
        
    if type(object_center_coordinates) is tuple and len(object_center_coordinates) != 4:
        raise TypeError('\n\n\'object_center_coordinates\' should be \'automatic\' or a tuple of length 4 containing floats or integers.')
    
    if type(object_center_coordinates) is tuple and any([type(x) not in [int, float] for x in object_center_coordinates]):
        raise TypeError('\n\n\'object_center_coordinates\' should be \'automatic\' or a tuple of length 4 containing floats or integers.')

    # center_param_centering
    if type(center_param_centering) is not tuple or len(center_param_centering) != 3:
        raise TypeError('\n\n\'center_param_centering\' should be a tuple of length 3 containing floats, integers or None.')
    
    if type(center_param_centering[0]) is not int and center_param_centering[0] is not None:
        raise TypeError('\n\nThe first element of \'center_param_centering\' (\'crop_radius\') should be a positive integer or None.')
    
    if type(center_param_centering[0]) is int and center_param_centering[0] <= 0:
        raise ValueError('\n\nThe first element of \'center_param_centering\' (\'crop_radius\') should be a positive integer or None.')
    
    if type(center_param_centering[1]) not in [int, float] and center_param_centering[1] is not None:
        raise TypeError('\n\nThe second element of \'center_param_centering\' (\'sigfactor\') should be a positive integer, positive float or None.')
    
    if type(center_param_centering[1]) in [int, float] and center_param_centering[1] <= 0:
        raise ValueError('\n\nThe second element of \'center_param_centering\' (\'sigfactor\') should be a positive integer, positive float or None.')
    
    if type(center_param_centering[2]) not in [int, float] and center_param_centering[2] is not None:
        raise TypeError('\n\nThe third element of \'center_param_centering\'(\'saturation_level\') should be a positive integer, positive float or None.')
    
    if type(center_param_centering[2]) in [int, float] and center_param_centering[2] <= 0:
        raise ValueError('\n\nThe third element of \'center_param_centering\'(\'saturation_level\') should be a positive integer, positive float or None.')
    
    # object_param_centering
    if type(object_param_centering) is not tuple or len(object_param_centering) != 3:
        raise TypeError('\n\n\'object_param_centering\' should be a tuple of length 3 containing floats, integers or None.')
    
    if type(object_param_centering[0]) is not int and object_param_centering[0] is not None:
        raise TypeError('\n\nThe first element of \'object_param_centering\' (\'crop_radius\') should be a positive integer or None.')
    
    if type(object_param_centering[0]) is int and object_param_centering[0] <= 0:
        raise ValueError('\n\nThe first element of \'object_param_centering\' (\'crop_radius\') should be a positive integer or None.')
    
    if type(object_param_centering[1]) not in [int, float] and object_param_centering[1] is not None:
        raise TypeError('\n\nThe second element of \'object_param_centering\' (\'sigfactor\') should be a positive integer, positive float or None.')
    
    if type(object_param_centering[1]) in [int, float] and object_param_centering[1] <= 0:
        raise ValueError('\n\nThe second element of \'object_param_centering\' (\'sigfactor\') should be a positive integer, positive float or None.')
    
    if type(object_param_centering[2]) not in [int, float] and object_param_centering[2] is not None:
        raise TypeError('\n\nThe third element of \'object_param_centering\'(\'saturation_level\') should be a positive integer, positive float or None.')
    
    if type(object_param_centering[2]) in [int, float] and object_param_centering[2] <= 0:
        raise ValueError('\n\nThe third element of \'object_param_centering\'(\'saturation_level\') should be a positive integer, positive float or None.')
    
    # flux_centering_method
    if flux_centering_method not in ['gaussian', 'manual']:
        raise ValueError('\n\n\'flux_centering_method\' should be either \'gaussian\' or \'manual\'.')
    
    # flux_center_coordinates
    if type(flux_center_coordinates) is not tuple or len(flux_center_coordinates) != 4:
        raise TypeError('\n\n\'flux_center_coordinates\' should be a tuple of length 4 containing floats or integers.')
    
    if any([type(x) not in [int, float] for x in flux_center_coordinates]):
        raise TypeError('\n\n\'flux_center_coordinates\' should be a tuple of length 4 containing floats or integers.')
    
    # flux_param_centering
    if type(flux_param_centering) is not tuple or len(flux_param_centering) != 3:
        raise TypeError('\n\n\'flux_param_centering\' should be a tuple of length 3 containing floats, integers or None.')
    
    if type(flux_param_centering[0]) is not int and flux_param_centering[0] is not None:
        raise TypeError('\n\nThe first element of \'flux_param_centering\' (\'crop_radius\') should be a positive integer or None.')
    
    if type(flux_param_centering[0]) is int and flux_param_centering[0] <= 0:
        raise ValueError('\n\nThe first element of \'flux_param_centering\' (\'crop_radius\') should be a positive integer or None.')
    
    if type(flux_param_centering[1]) not in [int, float] and flux_param_centering[1] is not None:
        raise TypeError('\n\nThe second element of \'flux_param_centering\' (\'sigfactor\') should be a positive integer, positive float or None.')
    
    if type(flux_param_centering[1]) in [int, float] and flux_param_centering[1] <= 0:
        raise ValueError('\n\nThe second element of \'flux_param_centering\' (\'sigfactor\') should be a positive integer, positive float or None.')
    
    if type(flux_param_centering[2]) not in [int, float] and flux_param_centering[2] is not None:
        raise TypeError('\n\nThe third element of \'flux_param_centering\'(\'saturation_level\') should be a positive integer, positive float or None.')
    
    if type(flux_param_centering[2]) in [int, float] and flux_param_centering[2] <= 0:
        raise ValueError('\n\nThe third element of \'flux_param_centering\'(\'saturation_level\') should be a positive integer, positive float or None.')
    
    # flux_annulus_background
    if flux_annulus_background == [] or flux_annulus_background == ():
        raise TypeError('\n\n\'flux_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')

    if type(flux_annulus_background) not in [str, tuple, list]:
        raise TypeError('\n\n\'flux_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(flux_annulus_background) is str and flux_annulus_background != 'large annulus':
        raise ValueError('\n\n\'flux_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(flux_annulus_background) is tuple and len(flux_annulus_background) != 6:
        raise TypeError('\n\n\'flux_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(flux_annulus_background) is tuple and any([type(x) not in [int, float] for x in flux_annulus_background]):
        raise TypeError('\n\n\'flux_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        
    if type(flux_annulus_background) is list:
        if any([type(x) is not tuple for x in flux_annulus_background]):
            raise TypeError('\n\n\'flux_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif any([len(x) != 6 for x in flux_annulus_background]):
            raise TypeError('\n\n\'flux_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif any([type(y) not in [int, float] for x in flux_annulus_background for y in x]):
            raise TypeError('\n\n\'flux_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    # save_preprocessed_data
    if save_preprocessed_data not in [True, False]:
        raise ValueError('\n\n\'save_preprocessed_data\' should be either True or False.')   
    
    # double_difference_type
    if double_difference_type not in ['conventional', 'normalized']:
        raise ValueError('\n\n\'double_difference_type\' should be either \'conventional\' or \'normalized\'.')
    
    # remove_vertical_band_detector_artefact    
    if remove_vertical_band_detector_artefact not in [True, False]:
        raise ValueError('\n\n\'remove_vertical_band_detector_artefact\' should be either True or False.')   
    
    # annulus_star   
    if type(annulus_star) not in [str, tuple, list]:
        raise TypeError('\n\n\'annulus_star\' should be \'automatic\', \'ao residuals\', \'star aperture\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')

    if annulus_star == [] or annulus_star == ():
        raise TypeError('\n\n\'annulus_star\' should be \'automatic\', \'ao residuals\', \'star aperture\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')

    if type(annulus_star) is str and annulus_star not in ['automatic', 'ao residuals', 'star aperture']:
        raise ValueError('\n\n\'annulus_star\' should be \'automatic\', \'ao residuals\', \'star aperture\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(annulus_star) is tuple and len(annulus_star) != 6:
        raise TypeError('\n\n\'annulus_star\' should be \'automatic\', \'ao residuals\', \'star aperture\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(annulus_star) is tuple and any([type(x) not in [int, float] for x in annulus_star]):
        raise TypeError('\n\n\'annulus_star\' should be \'automatic\', \'ao residuals\', \'star aperture\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        
    if type(annulus_star) is list:        
        if any([type(x) is not tuple for x in annulus_star]):
            raise TypeError('\n\n\'annulus_star\' should be \'automatic\', \'ao residuals\', \'star aperture\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif any([len(x) != 6 for x in annulus_star]):
            raise TypeError('\n\n\'annulus_star\' should be \'automatic\', \'ao residuals\', \'star aperture\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif any([type(y) not in [int, float] for x in annulus_star for y in x]):
            raise TypeError('\n\n\'annulus_star\' should be \'automatic\', \'ao residuals\', \'star aperture\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    # annulus_background
    if type(annulus_background) not in [str, tuple, list]:
        raise TypeError('\n\n\'annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')

    if annulus_background == [] or annulus_background == ():
        raise TypeError('\n\n\'annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        
    if type(annulus_background) is str and annulus_background != 'large annulus':
        raise ValueError('\n\n\'annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(annulus_background) is tuple and len(annulus_background) != 6:
        raise TypeError('\n\n\'annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(annulus_background) is tuple and any([type(x) not in [int, float] for x in annulus_background]):
        raise TypeError('\n\n\'annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        
    if type(annulus_background) is list:
        if any([type(x) is not tuple for x in annulus_background]):
            raise TypeError('\n\n\'annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif any([len(x) != 6 for x in annulus_background]):
            raise TypeError('\n\n\'annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif any([type(y) not in [int, float] for x in annulus_background for y in x]):
            raise TypeError('\n\n\'annulus_background\' should be \'large annulus\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')

    #flux_annulus_star 
    if type(flux_annulus_star) not in [str, tuple, list]:
        raise TypeError('\n\n\'flux_annulus_star\' should be \'automatic\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')

    if flux_annulus_star == [] or flux_annulus_star == ():
        raise TypeError('\n\n\'flux_annulus_star\' should be \'automatic\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')

    if type(flux_annulus_star) is str and flux_annulus_star not in ['automatic']:
        raise ValueError('\n\n\'flux_annulus_star\' should be \'automatic\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(flux_annulus_star) is tuple and len(flux_annulus_star) != 6:
        raise TypeError('\n\n\'flux_annulus_star\' should be \'automatic\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
    
    if type(flux_annulus_star) is tuple and any([type(x) not in [int, float] for x in flux_annulus_star]):
        raise TypeError('\n\n\'flux_annulus_star\' should be \'automatic\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        
    if type(flux_annulus_star) is list:        
        if any([type(x) is not tuple for x in flux_annulus_star]):
            raise TypeError('\n\n\'flux_annulus_star\' should be \'automatic\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif any([len(x) != 6 for x in flux_annulus_star]):
            raise TypeError('\n\n\'flux_annulus_star\' should be \'automatic\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif any([type(y) not in [int, float] for x in flux_annulus_star for y in x]):
            raise TypeError('\n\n\'flux_annulus_star\' should be \'automatic\', a length-6 tuple of floats or integers or a list of length-6 tuples of floats or integers.')
        elif flux_annulus_star[2] != 0:
            printandlog('Warning, the aperture inner radius for the star photometry (3rd element of flux_annulus_star) is set to {0:.1f}. It should be set to 0 to encompass the star flux unless very specific circumstances'.format(flux_annulus_star[2]))        
        elif flux_annulus_star[3] >140:
            printandlog('Warning, the aperture outer radius for the star photometry (4th element of flux_annulus_star) is set to {0:.1f}. A cluster of bad pixels is present between 140px and 160px that might bias the photometry of the star'.format(flux_annulus_star[3]))        

    # combination_method_polarization
    if combination_method_polarization not in ['least squares', 'trimmed mean', 'median']:
        raise ValueError('\n\n\'combination_method_polarization\' should be \'least squares\', \'trimmed mean\' or \'median\'.')
    
    # trimmed_mean_prop_to_cut_polar    
    if type(trimmed_mean_prop_to_cut_polar) not in [int, float]:
        raise TypeError('\n\n\'trimmed_mean_prop_to_cut_polar\' should be of type int or float.')
    
    if not 0 <= trimmed_mean_prop_to_cut_polar <= 1:
        raise ValueError('\n\n\'trimmed_mean_prop_to_cut_polar\' should be in range 0 <= trimmed_mean_prop_to_cut_polar <= 1.')
    
    # combination_method_intensity
    if combination_method_intensity not in ['mean', 'trimmed mean', 'median']:
        raise ValueError('\n\n\'combination_method_intensity\' should be \'mean\', \'trimmed mean\' or \'median\'.')
    
    # trimmed_mean_prop_to_cut_intens
    if type(trimmed_mean_prop_to_cut_intens) not in [int, float]:
        raise TypeError('\n\n\'trimmed_mean_prop_to_cut_intens\' should be of type int or float.')
    
    if not 0 <= trimmed_mean_prop_to_cut_intens <= 1:
        raise ValueError('\n\n\'trimmed_mean_prop_to_cut_intens\' should be in range 0 <= trimmed_mean_prop_to_cut_intens <= 1.')
    
    # single_posang_north_up
    if single_posang_north_up not in [True, False]:
        raise ValueError('\n\n\'single_posang_north_up\' should be either True or False.')   
    
    # convert_in_contrast_per_arcsec2
    if convert_in_contrast_per_arcsec2 not in [True, False]:
        raise ValueError('\n\n\'convert_in_contrast_per_arcsec2\' should be either True or False.')   
    
    # normalized_polarization_images
    if normalized_polarization_images not in [True, False]:
        raise ValueError('\n\n\'normalized_polarization_images\' should be either True or False.')   

    printandlog('\nThe input parameters have passed all checks.')

    ###############################################################################
    # Convert input from 1-based to 0-based indexing
    ###############################################################################
    
    # object_center_coordinates
    if type(object_center_coordinates) is tuple:
        object_center_coordinates = tuple(x - 1 for x in object_center_coordinates)

    # flux_center_coordinates
    flux_center_coordinates = tuple(x - 1 for x in flux_center_coordinates)
    
    # flux_annulus_background
    flux_annulus_background = annulus_1_to_0_based(flux_annulus_background)
    
    # frames_to_remove
    for i,x in enumerate(frames_to_remove):
        if type(x) is tuple:
            x = tuple(y - 1 for y in x)
        else:
            x -= 1
        frames_to_remove[i] = x
    
    # annulus_star
    annulus_star = annulus_1_to_0_based(annulus_star)

    # annulus_background
    annulus_background = annulus_1_to_0_based(annulus_background)
    
    ###############################################################################
    # Run pre-processing and post-processing functions
    ###############################################################################
       
    if skip_preprocessing == False:
        # Pre-process raw data
        printandlog('\n###############################################################################')
        printandlog('# Starting pre-processing')
        printandlog('###############################################################################')
        printandlog('\nStarting pre-processing of raw data.')
        
        cube_left_frames, cube_right_frames, header, file_index_object,\
        combination_method_polarization, reference_flux \
        = perform_preprocessing(frames_to_remove=frames_to_remove, 
                                sigma_filtering=sigma_filtering, 
                                object_collapse_ndit=object_collapse_ndit, 
                                show_images_center_coordinates=show_images_center_coordinates, 
                                object_centering_method=object_centering_method, 
                                center_subtract_object=center_subtract_object, 
                                object_center_coordinates=object_center_coordinates, 
                                center_param_centering=center_param_centering,
                                object_param_centering=object_param_centering,
                                flux_centering_method=flux_centering_method, 
                                flux_center_coordinates=flux_center_coordinates, 
                                flux_param_centering=flux_param_centering, 
                                flux_annulus_background=flux_annulus_background, 
                                save_preprocessed_data=save_preprocessed_data, 
                                combination_method_polarization=combination_method_polarization,\
                                flux_annulus_star = flux_annulus_star,\
                                convert_in_contrast_per_arcsec2=convert_in_contrast_per_arcsec2)
    
        # Print that post-processing starts
        printandlog('\n###############################################################################')
        printandlog('# Starting post-processing')
        printandlog('###############################################################################') 
        printandlog('\nContinuing with the pre-processed data.')
                    
    elif skip_preprocessing == True:       
        # Write lines from pre-processing of previous log file
        for line in log_file_lines:                         
            print(line, file=open(path_log_file, 'a'))
        
        # Define paths to read pre-processed data and headers from
        path_cube_left_frames = os.path.join(path_preprocessed_dir, 'cube_left_frames.fits')
        path_cube_right_frames = os.path.join(path_preprocessed_dir, 'cube_right_frames.fits')
        path_object_files_text = os.path.join(path_preprocessed_dir, 'path_object_files.txt')
        path_file_index_object = os.path.join(path_preprocessed_dir, 'file_index_object.txt')
        
        if os.path.exists(path_cube_left_frames) and os.path.exists(path_cube_right_frames) and os.path.exists(path_object_files_text) and os.path.exists(path_file_index_object):
            # Print that post-processing starts
            printandlog('\n###############################################################################')
            printandlog('# Starting post-processing')
            printandlog('###############################################################################') 
            printandlog('\nSkipping pre-processing and reading pre-processed data and headers.')
            printandlog('')
            
            # Read pre-processed single-sum and difference- images        
            cube_left_frames = read_fits_files(path=path_cube_left_frames, silent=False)[0]
            cube_right_frames = read_fits_files(path=path_cube_right_frames, silent=False)[0] 
    
            # Read headers
            header = [pyfits.getheader(x.rstrip('\n')) for x in open(path_object_files_text, 'r')]
            printandlog('Read headers from OBJECT-files specified in ' + path_object_files_text + '.')
            
            # Read indices of OBJECT-files
            file_index_object = literal_eval([x for x in open(path_file_index_object, 'r')][0])
            printandlog('Read indices of OBJECT-files from ' + path_file_index_object + '.')

            printandlog('\nRetained part of previous log file pertaining to pre-processing.')
                        
        else:
            raise IOError('\n\nThe files ' + path_cube_left_frames + ', ' + path_cube_right_frames + ', ' + path_object_files_text + ' and/or ' + path_file_index_object + ' do not exist. Run IRDAP first with \'skip_preprocessing\' = False to perform the pre-processing of the raw data and save the results.')
    
    # Perform post-processing of data           
    perform_postprocessing(cube_left_frames=cube_left_frames, 
                           cube_right_frames=cube_right_frames, 
                           header=header, 
                           file_index_object=file_index_object, 
                           annulus_star=annulus_star, 
                           annulus_background=annulus_background, 
                           double_difference_type=double_difference_type, 
                           remove_vertical_band_detector_artefact=remove_vertical_band_detector_artefact, 
                           combination_method_polarization=combination_method_polarization, 
                           trimmed_mean_prop_to_cut_polar=trimmed_mean_prop_to_cut_polar, 
                           combination_method_intensity=combination_method_intensity, 
                           trimmed_mean_prop_to_cut_intens=trimmed_mean_prop_to_cut_intens,
                           single_posang_north_up=single_posang_north_up, 
                           normalized_polarization_images=normalized_polarization_images)

#TODO: Add main ADI function
#    perform_adi()
    
    # Print time elapsed
    time_end = time.time()
    d = datetime.datetime(1, 1, 1) + datetime.timedelta(seconds = time_end - time_start)
    printandlog('\nTime elapsed: %d h %d min %d s' % (d.hour, d.minute, d.second)) 
    
    # Print request to cite IRDAP
    printandlog('\n###############################################################################')
    printandlog('# Important notice')
    printandlog('###############################################################################')
    printandlog('\nWhen publishing data reduced with IRDAP, please cite van Holstein et al.')
    printandlog('(2019): <ADS link>.')
    printandlog('For data in pupil-tracking mode please additionally cite van Holstein et al.')
    printandlog('(2017): http://adsabs.harvard.edu/abs/2017SPIE10400E..15V.')
    printandlog('\nFull documentation: https://robvanholstein.github.io/IRDAP')
    printandlog('Feedback, questions, comments: vanholstein@strw.leidenuniv.nl')
    printandlog('\nIRDAP Copyright (C) 2019 R.G. van Holstein')