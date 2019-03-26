'''

@author: Rob the Modellist and Christian the Artist
'''

###############################################################################
# Input
###############################################################################

# Definition of target name, observation date and paths
target_name = 'TEST'
date_obs = '0000'
path_main_dir = r'C:\Users\Rob\Desktop\IP Measurement Test\GQ Lup'
path_static_flat_badpixelmap = r'C:\Users\Rob\Documents\PhD\CentralFiles\irdis_polarimetry_pipeline'

# Options for pre-processing
skip_preprocessing = False

#TODO: Write parts of readme
#TODO: Keep this or always do sigmafiltering? If keep then write part README
sigmafiltering_sky = True
sigmafiltering_object = True
sigmafiltering_flux = True
param_annulus_background_flux = 'large annulus'

save_preprocessed_data = True

# Options for post-processing
double_difference_type = 'standard'
remove_vertical_band_detector_artefact = True
param_annulus_star = 'ao residuals'
param_annulus_background = 'large annulus'
combination_method_polarization_images = 'trimmed mean'
trimmed_mean_proportiontocut_polarization_images = 0.10
combination_method_total_intensity_images = 'mean'
trimmed_mean_proportiontocut_total_intensity_images = 0.10
images_north_up = True
create_images_DoLP_AoLP_q_u_norm = True

#TODO: turn combination_method_polarization_images and trimmed_mean_proportiontocut_polarization_images
# into a single variable (also for the total intensity ones)? Instead of trimmed mean, give a number, 
# also for mean can probably give 0 (should test it though).


###############################################################################
# README
###############################################################################

'''
target_name:

String of name of target, e.g. 'T_Cha'. The target name and observation date 
(see below) are used at the beginning of the name of the reduced files 
generated. 

date_obs:
    
String of date of observations, e.g. '2018-03-23'. The target name (see below) 
and observation date are used at the beginning of the name of the reduced files
generated. 

path_main_dir:
    
String of path to location of data directory, e.g. 'C:\Data\T Cha'. The raw
data, i.e. FITS-files containing the OBJECT, CENTER, SKY and/or FLUX images, 
need to be present in a subdirectory 'Raw' within this directory, so for the 
example above in 'C:\Data\T Cha\Raw'. One can use the directory separator 
('\' above) specific to the operating system in use. The reduced calibration 
and object data will be written to subdirectories called 'Calibrations' and 
'Reduced', respectively.

path_static_flat_badpixelmap:

String of path to location of directory containing the FITS-files of the static
master flats and static bad pixel map, e.g. C:\Data\CentralFiles'. It might be 
useful to keep these static files in a central directory so that they are 
easily accessible when reducing various data sets.

save_preprocessed_data:
    
If True, save preprocessed cubes of single-sum and single-difference images in 
the 'preprocessed' folder so that the preprocessing can be skipped when 
re-running the pipeline.

double_difference_type:

Type of double difference to be computed, either 'standard' or 'normalized' 
(see van Holstein et al. 2019). In almost all cases one would use 'standard'.
When there are large variations in seeing and sky transparency among the
measurements, using 'normalized' can suppress spurious polarization signals and
improve the quality of the final images.

remove_vertical_band_detector_artefact:

If True remove the vertical band detector artefact seen in the double-
difference Q- and U-images by subtracting the median of the top and bottom 60
pixel of each pixel column in the images. If False don't remove the artefact.
In almost all cases one would use True. One would only use False in case 
there is a lot of astrophysical signal within in the region where the median 
is computed.

param_annulus_star: 
    
Parameter(s) defining with which annulus/annuli the star polarization will be
determined. If 'ao residuals' the annulus will be star-centered and located 
over the AO residuals. The inner radius and width of the annulus will depend on
the filter used. If 'star aperture' a small aparture located at the position of
the central star will be used. One can manually define the annulus/annuli by
specifying a (list of) length-6-tuple(s) of floats with parameters:
- x-coordinate of center (pixels, 0-based, center of image is at 511.5)
- y-coordinate of center (pixels, 0-based, center of image is at 511.5)
- inner radius (pixels)
- width (pixels)
- start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
- end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
For example, when using 'star aperture' the annulus used will be: 
(511.5, 511.5, 0, 11, 0, 360). The coordinates and angles of the annulus are
defined with respect to the final orientation of the image. Generally this is 
North up, except when images_north_up = False (see below) for observations in 
field-tracking mode with a single derotator POSANG. The annulus used can be
checked with annulus_star.fits that can be found in the directories with the
reduced images.

param_annulus_background: 
    
Parameter(s) defining with which annulus/annuli the background will be
determined. If 'large annulus' the annulus will be star-centered and located 
far away from the star with an inner radius of 360 pixels and a width of 60
pixels. One can manually define the annulus/annuli by specifying a (list of) 
length-6-tuple(s) of floats with parameters:
- x-coordinate of center (pixels, 0-based, center of image is at 511.5)
- y-coordinate of center (pixels, 0-based, center of image is at 511.5)
- inner radius (pixels)
- width (pixels)
- start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
- end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
For example, when using 'large annulus' the annulus used will be: 
(511.5, 511.5, 360, 60, 0, 360). The coordinates and angles of the annulus are
defined with respect to the final orientation of the image. Generally this is 
North up, except when images_north_up = False (see below) for observations in 
field-tracking mode with a single derotator POSANG. The annulus used can be
checked with annulus_background.fits that can be found in the directories with 
the reduced images.

combination_method_polarization_images:

Method to be used to produce the incident Q- and U-images, i.e. the images that
are corrected for the instrumental polarization effects. Valid values are 
'least squares', 'trimmed mean' or 'median'. The recommended option is 'trimmed
mean.' With 'least squares' the images are obtained by solving for every pixel 
the system of equations describing the measurements using linear least squares 
(see Eq. 35 of van Holstein et al. 2019). With 'trimmed mean' or 'median' the 
images are obtained by solving the system of equations for each pair of double-
difference Q- and U-images (each HWP cycle) separately, and then computing the 
trimmed mean or median over all resulting images. 'least squares' is the most 
accurate option, but any unremoved bad pixels will still be visible in the 
images. Using 'median' will remove these bad pixels, but is the least accurate 
option and also yields images with a lower signal-to-noise ratio as is clear
from images of circumstellar disks. Using 'trimmed mean' will yield images that
have essentially the same accuracy and signal-to-noise ratio images produced 
using 'least squares', but without the bad pixels. Therefore 'trimmed mean' is 
the recommended option.

trimmed_mean_proportiontocut_polarization_images:

Fraction to cut off of both tails of the distribution if 'trimmed mean' is used
for combination_method_polarization_images. Parameter is ignored in case 
'least squares' or 'median' is used. Value should be in range
0 <= trimmed_mean_proportiontocut_polarization_images <= 1. In most cases a 
value of 0.1 or 0.15 removes the bad pixels well while producing images very 
similar to those obtained with 'least squares'.
    
combination_method_total_intensity_images:

Method to be used to produce the incident I_Q- and I_U-images. These images are 
computed by combining the I_Q- or I_U-images of all HWP cycles using the 
'mean', 'trimmed mean' or 'median'. 'mean' yields the most accurate images, but
any unremoved bad pixels will still be visible in the images. Using 'median' 
will remove these bad pixels, but is the least accurate option. 'trimmed mean'
produces images similar to 'mean', but without the bad pixels. It is generally
recommended to use either 'trimmed mean' or 'mean'.

trimmed_mean_proportiontocut_total_intensity_images:

Fraction to cut off of both tails of the distribution if 'trimmed mean' is used
for combination_method_total_intensity_images. Parameter is ignored in case 
'mean' or 'median' is used. Value should be in range
0 <= trimmed_mean_proportiontocut_total_intensity_images <= 1. In most cases a 
value of 0.1 or 0.15 removes the bad pixels well while producing images similar
to those obtained with 'mean'.

images_north_up:

For observations taken in field-tracking mode with a single derotator position 
angle (INS4.DROT2.POSANG either 0 or with a fixed offset), the final images 
are oriented with North up if True, and oriented as the raw frames if False.
Parameter is ignored for field-tracking observations with more than one 
derotator position angle or observations taken in pupil-tracking mode, because 
in these cases the final images produced always have North up.

create_images_DoLP_AoLP_q_u_norm: 
    
if True create final images of degree of linear polarization, normalized Stokes
q and u and degree and angle of linear polarization computed from q and u:
- DoLP = sqrt(Q^2 + U^2) / (0.5*(I_Q + I_U))
- q = Q / I_Q
- u = U / I_U
- AoLP_norm = 0.5 * arctan(u / q)
- DoLP_norm = sqrt(q^2 + u^2)
These images are only valid if all flux in the images originates from the 
astrophysical source of interest. This is generally the case for observations 
of for example solar system objects or galaxies. The images are generally not 
valid for observations of circumstellar disks or companions because in that 
case a large part of the flux in the total intensity images originates from the 
central star. AoLP_norm and DoLP_norm are potentially more accurate than 
AoLP = 0.5 * arctan(U / Q) and DoLP, especially when there are significant 
variations in seeing and sky transparency among the measurements.


'''

###############################################################################
# Import packages
###############################################################################

# Import packages
import os
import glob
import time
import datetime
import math
import subprocess
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy import optimize
from scipy import ndimage
from scipy.stats import trim_mean
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model

#import pdb # for debugging

###############################################################################
# Start taking time
###############################################################################

time_start = time.time()

###############################################################################
# Check whether input values are valid
###############################################################################

if type(target_name) is not str:
    raise TypeError('\'target_name\' should be of type string.')

if type(date_obs) is not str:
    raise TypeError('\'date_obs\' should be of type string.')

if type(path_main_dir) is not str:
    raise TypeError('\'path_main_dir\' should be of type string.')

if type(path_static_flat_badpixelmap) is not str:
    raise TypeError('\'path_static_flat_badpixelmap\' should be of type string.')

if double_difference_type not in ['standard', 'normalized']:
    raise ValueError('\'double_difference_type\' should be either \'standard\' or \'normalized\'.')
    
if remove_vertical_band_detector_artefact not in [True, False]:
    raise ValueError('\'remove_vertical_band_detector_artefact\' should be either True or False.')   

#TODO: Add checks of specific values of param_annulus_star/background
#      Check how an annulus believes when it is cut off by the edge of the image
if type(param_annulus_star) not in [str, tuple, list]:
    raise TypeError('\'param_annulus_star\' should be \'ao residuals\', \'star aperture\', a length-6 tuple of floats or a list of length-6 tuples of floats.')

elif type(param_annulus_star) is str and param_annulus_star not in ['ao residuals', 'star aperture']:
    raise ValueError('\'param_annulus_star\' should be \'ao residuals\', \'star aperture\', a length-6 tuple of floats or a list of length-6 tuples of floats.')

elif type(param_annulus_star) is tuple and len(param_annulus_star) is not 6:
    raise ValueError('\'param_annulus_star\' should be \'ao residuals\', \'star aperture\', a length-6 tuple of floats or a list of length-6 tuples of floats.')

elif type(param_annulus_star) is list:
    if any([type(x) is not tuple for x in param_annulus_star]):
        raise TypeError('\'param_annulus_star\' should be \'ao residuals\', \'star aperture\', a length-6 tuple of floats or a list of length-6 tuples of floats.')
    elif any([len(x) is not 6 for x in param_annulus_star]):
        raise ValueError('\'param_annulus_star\' should be \'ao residuals\', \'star aperture\', a length-6 tuple of floats or a list of length-6 tuples of floats.')

if type(param_annulus_background) not in [str, tuple, list]:
    raise TypeError('\'param_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or a list of length-6 tuples of floats.')

elif type(param_annulus_background) is str and param_annulus_background not in ['large annulus']:
    raise ValueError('\'param_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or a list of length-6 tuples of floats.')

elif type(param_annulus_background) is tuple and len(param_annulus_background) is not 6:
    raise ValueError('\'param_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or a list of length-6 tuples of floats.')
    
elif type(param_annulus_background) is list:
    if any([type(x) is not tuple for x in param_annulus_background]):
        raise TypeError('\'param_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or a list of length-6 tuples of floats.')
    elif any([len(x) is not 6 for x in param_annulus_background]):
        raise ValueError('\'param_annulus_background\' should be \'large annulus\', a length-6 tuple of floats or a list of length-6 tuples of floats.')

if combination_method_polarization_images not in ['least squares', 'trimmed mean', 'median']:
    raise ValueError('\'combination_method_polarization_images\' should be \'least squares\', \'trimmed mean\' or \'median\'.')
    
if type(trimmed_mean_proportiontocut_polarization_images) not in [int, float]:
    raise TypeError('\'trimmed_mean_proportiontocut_polarization_images\' should be of type int or float.')

if not 0 <= trimmed_mean_proportiontocut_polarization_images <= 1:
    raise ValueError('\'trimmed_mean_proportiontocut_polarization_images\' should be in range 0 <= trimmed_mean_proportiontocut_polarization_images <= 1.')

if combination_method_total_intensity_images not in ['mean', 'trimmed mean', 'median']:
    raise ValueError('\'combination_method_total_intensity_images\' should be \'mean\', \'trimmed mean\' or \'median\'.')

if type(trimmed_mean_proportiontocut_total_intensity_images) not in [int, float]:
    raise TypeError('\'trimmed_mean_proportiontocut_total_intensity_images\' should be of type int or float.')

if not 0 <= trimmed_mean_proportiontocut_total_intensity_images <= 1:
    raise ValueError('\'trimmed_mean_proportiontocut_total_intensity_images\' should be in range 0 <= trimmed_mean_proportiontocut_total_intensity_images <= 1.')

if images_north_up not in [True, False]:
    raise ValueError('\'images_north_up\' should be either True or False.')   

if create_images_DoLP_AoLP_q_u_norm not in [True, False]:
    raise ValueError('\'create_images_DoLP_AoLP_q_u_norm\' should be either True or False.')   

###############################################################################
# Define global variables
###############################################################################

# Define pupil-offset (deg) in pupil-tracking mode (SPHERE User Manual P99.0, 6th public release, P99 Phase 1)
pupil_offset = 135.99
    
# Define true North correction (deg) in pupil-tracking mode (SPHERE User Manual P99.0, 6th public release, P99 Phase 1)
true_north_correction = -1.75

# Define the base of the name of each file to be generated
name_file_root = target_name.replace(' ', '_') + '_' + date_obs.replace(' ', '_') + '_'

# Define directory names
path_raw_dir = os.path.join(path_main_dir, 'raw')
path_calib_dir = os.path.join(path_main_dir, 'calibration')
path_sky_dir = os.path.join(path_calib_dir, 'sky')
path_center_dir = os.path.join(path_calib_dir, 'center')
path_flux_dir = os.path.join(path_calib_dir, 'flux')
path_sky_flux_dir = os.path.join(path_calib_dir, 'sky_flux')
path_preprocessed_dir = os.path.join(path_main_dir, 'preprocessed')
path_reduced_dir = os.path.join(path_main_dir, 'reduced')
path_reduced_star_pol_subtr_dir = os.path.join(path_main_dir, 'reduced_star_pol_subtr')

###############################################################################
###############################################################################
## Function definitions
###############################################################################
###############################################################################

###############################################################################
# printandlog
###############################################################################

def printandlog(single_object):
    '''
    Print a single string or float on screen and save it to a log file. The log
    file is located at path_main_dir and the name starts with name_file_root, 
    which are strings and global variables to the funtion.
    
    Input:
        objects: string to be printen and logged
           
    File written by Rob van Holstein
    Function status: verified
    '''
            
    # Define path to write log file, using path_reduced_dir as global variable
    path_log_file = os.path.join(path_main_dir, name_file_root + 'LOG.txt')

    if not os.path.exists(path_log_file):
        # Create log file
        open(path_log_file, 'w+')
    
    # Print object in log file and on screen
    print(single_object, file=open(path_log_file, 'a'))
    print(single_object)

###############################################################################
# check_sort_data_create_directories
###############################################################################

def check_sort_data_create_directories():
    '''
    Check the FITS-headers of the data in the raw directory, sort the data and
    create directories to write processed data to.
    
    Note that path_raw_dir, path_sky_dir, path_center_dir, path_flux_dir, 
    path_sky_flux_dir, path_preprocessed_dir, path_reduced_dir, 
    path_reduced_star_pol_subtr_dir and 
    save_preprocessed_data are global variables to the function.
    
    Output:
        path_object_files: list of paths to raw OBJECT-files
        path_sky_files: list of paths to raw SKY-files for OBJECT
        path_center_files: list of paths to raw CENTER-files
        path_flux_files: list of paths to raw FLUX-files
        path_sky_flux_files: list of paths to raw SKY-files for FLUX
    
    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: verified
    '''

    # Check if raw directory exists, if not create it
    if not os.path.exists(path_raw_dir):
        os.makedirs(path_raw_dir)
        raise IOError('The raw directory {0:s} did not exist. It was created but you need to put your raw FITS-files there.'.format(path_raw_dir))
    
    # Extract paths to FITS-files in raw directory
    path_raw_files = glob.glob(os.path.join(path_raw_dir,'*.fits'))
    
    # Check if raw folder contains FITS-files
    if len(path_raw_files) == 0:
        raise IOError('The raw directory {0:s} does not contain FITS-files. You need to put your raw FITS-files in this folder.'.format(path_raw_dir))
        
    # Extract headers
    header = [pyfits.getheader(x) for x in path_raw_files]
    
    # Sort raw files and headers based on observation date in headers
    date_obs = [x['DATE-OBS'] for x in header]
    sort_index = list(np.argsort(date_obs))
    path_raw_files = [path_raw_files[i] for i in sort_index]
    header = [header[i] for i in sort_index]
    
    # Perform checks on header values that apply to all data
    printandlog('\nChecking the headers of the files in the raw directory\n{0:s}.'.format(path_raw_dir))
    
    if not all([x['ESO DPR TYPE'] in ['OBJECT', 'SKY', 'OBJECT,CENTER', 'OBJECT,FLUX'] for x in header]):
        raise IOError('One or more files are not of type OBJECT, SKY, OBJECT,CENTER or OBJECT,FLUX.')
        
    if not all([x['ESO INS4 OPTI8 NAME'] == 'H_NIR' for x in header]):
        raise IOError('One or more files do not have the NIR half-wave plate inserted.')
    
    if len(set([x['ESO INS1 FILT ID'] for x in header])) != 1:
        raise IOError('The data provided use different filters.')
    
    if not header[0]['ESO INS1 FILT ID'] in ['FILT_BBF_Y', 'FILT_BBF_J', 'FILT_BBF_H', 'FILT_BBF_Ks']:
        raise IOError('The filter used is not broadband Y, J, H or Ks. Narrowband filters will be supported in the near future.')
    
    if not all([x['ESO INS1 OPTI2 NAME'] == 'P0-90' for x in header]):
        raise IOError('One or more files do not have the P0-90 polarizer set inserted.')
    
    # Perform checks on header values that apply only to OBJECT files
    header_object = [x for x in header if x['ESO DPR TYPE'] == 'OBJECT']
        
    if not any(header_object):
        raise IOError('There are no OBJECT-files.')
        
    if len(set([x['EXPTIME'] for x in header_object])) != 1:
        raise IOError('The OBJECT-files have different exposure times.')
    
    if len(set([x['ESO INS4 FILT2 NAME'] for x in header_object])) != 1:
        raise IOError('The OBJECT-files use different NIR neutral density filters.')
    
    if header_object[0]['ESO INS4 FILT2 NAME'] != 'OPEN':
        printandlog('\nWARNING, the OBJECT-files use a NIR neutral density filter. The final data \nproduct will be less accurate because the neutral density filters are known to \nhave a depolarizing effect that is not calibrated.')
    
    if len(set([x['ESO INS COMB ICOR'] for x in header_object])) != 1:
        raise IOError('The OBJECT-files use different coronagraph settings.')
    
    # Determine exposure time and NIR neutral density filter for OBJECT files  
    object_exposure_time = header_object[0]['EXPTIME']
    object_nd_filter = header_object[0]['ESO INS4 FILT2 NAME']
    
    # Perform checks on header values that apply only to FLUX files
    header_flux = [x for x in header if x['ESO DPR TYPE'] == 'OBJECT,FLUX']
    
    if any(header_flux):
        if len(set([x['EXPTIME'] for x in header_flux])) != 1:
            raise IOError('The FLUX-files have different exposure times.')
        
        if len(set([x['ESO INS4 FILT2 NAME'] for x in header_flux])) != 1:
            raise IOError('The FLUX-files use different NIR neutral density filters.')
        
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
                raise IOError('One or more SKY-files have an exposure time different from that of the OBJECT- or FLUX-files.')
            else:
                raise IOError('One or more SKY-files have an exposure time different from that of the OBJECT-files.')
    
        if not all([x['ESO INS4 FILT2 NAME'] in [object_nd_filter, flux_nd_filter] for x in header_sky]):
            if any(header_flux):
                raise IOError('One or more SKY-files use a NIR neutral density filter different from that of the OBJECT- or FLUX-files.')
            else:
                raise IOError('One or more SKY-files use a NIR neutral density filter different from that of the OBJECT-files.')
    
    if not any([x['ESO DPR TYPE'] == 'SKY' and x['EXPTIME'] == object_exposure_time and x['ESO INS4 FILT2 NAME'] == object_nd_filter for x in header]):
        # TODO: Change statement here when pipeline is finished
        printandlog('\nWARNING, there are no SKY-files to subtract from the OBJECT-files. Although the \nbackground will be subtracted from the final images after determining it using \nthe annulus as defined by the input variable \'param_annulus_background\', the \nresult will be less accurate than when subtracting a SKY-image.')    
    
    if any(header_flux):    
        if not any([x['ESO DPR TYPE'] == 'SKY' and x['EXPTIME'] == flux_exposure_time and x['ESO INS4 FILT2 NAME'] == flux_nd_filter for x in header]):
            # TODO: Change statement here when flux part of pipeline is finished + change readme on what param_annulus_background is used for
            printandlog('\nWARNING, there are no SKY-files to subtract from the FLUX-file(s). Although the \nbackground will be subtracted after determining it using the annulus as defined \nby the input variable \'param_annulus_background\', the result will be less \naccurate than when subtracting a SKY-image.')
    
    # Perform checks on header values that apply only to CENTER files
    header_center = [x for x in header if x['ESO DPR TYPE'] == 'OBJECT,CENTER']
    
    if any(header_center):
        if not all([x['EXPTIME'] == object_exposure_time for x in header_center]):
            raise IOError('One or more CENTER-files have an exposure time different from that of the OBJECT-files.')
        
        if not all([x['ESO INS4 FILT2 NAME'] == object_nd_filter for x in header_center]):
            raise IOError('One or more CENTER-files use a NIR neutral density filter different from that of the OBJECT-files.')
    
    else:
        #TODO: Change statement here when centering part of pipeline is finished
        printandlog('\nWARNING, there are no CENTER-files. Centering of the OBJECT-files will be done \nby cross-correlation or from manual input.')
        
    # Print that headers have been checked and passed all checks
    printandlog('\nThe FITS-headers of the raw data have passed all checks.')
        
    # Create list of file paths for each file type
    path_object_files = []
    path_sky_files = []
    path_center_files = []
    path_flux_files = []
    path_sky_flux_files = [] 
    path_imcompatible_files = []
    
    # Sort file paths according to file type 
    for file_sel, header_sel in zip(path_raw_files, header):
    
        if header_sel['ESO DPR TYPE'] == 'OBJECT':
            path_object_files.append(file_sel)
    
        elif header_sel['ESO DPR TYPE'] == 'SKY' and \
           header_sel['EXPTIME'] == object_exposure_time and \
           header_sel['ESO INS4 FILT2 NAME'] == object_nd_filter:
            path_sky_files.append(file_sel)
    
        elif header_sel['ESO DPR TYPE'] == 'OBJECT,CENTER' and \
           header_sel['EXPTIME'] == object_exposure_time and \
           header_sel['ESO INS4 FILT2 NAME'] == object_nd_filter:
            path_center_files.append(file_sel)
    
        elif header_sel['ESO DPR TYPE'] == 'OBJECT,FLUX':
            path_flux_files.append(file_sel)
    
        elif header_sel['ESO DPR TYPE'] == 'SKY' and \
           header_sel['EXPTIME'] == flux_exposure_time and \
           header_sel['ESO INS4 FILT2 NAME'] == flux_nd_filter:
            path_sky_flux_files.append(file_sel)
    
        else:
            path_imcompatible_files.append(file_sel)
            
    # Print number of files for each file type
    printandlog('\nNumber of files found for each file type:')
    printandlog('OBJECT (O):         ' + str(len(path_object_files)))
    printandlog('SKY (S) for OBJECT: ' + str(len(path_sky_files)))
    printandlog('CENTER (C):         ' + str(len(path_center_files)))
    printandlog('FLUX (F):           ' + str(len(path_flux_files)))
    printandlog('SKY (S) for FLUX:   ' + str(len(path_sky_flux_files)))
    
    # Print warning when one or more raw files do not fall under any of the categories above
    if len(path_imcompatible_files) != 0:
        printandlog('\nWARNING, one or more files do not fall under any of the file type categories listed above and will be ignored:')
        for file_sel in path_imcompatible_files:
            printandlog('{0:s}'.format(file_sel))
        
    # Create directories to write processed data to
    directories_created = []
    directories_already_existing = []
    
    if any(path_sky_files):
        if not os.path.exists(path_sky_dir):
            os.makedirs(path_sky_dir)
            directories_created.append(path_sky_dir)
        else:
            directories_already_existing.append(path_sky_dir)
    
    if any(path_center_files):
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
       
    if save_preprocessed_data == True:
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
        #TODO: Change statement when pipeline is finished as then not all data will be overwritten probably, e.g. the calibration data
        printandlog('\nThe following directories already exist. Data in these directories will be \noverwritten:')
        for directory_sel in directories_already_existing:
            printandlog('{0:s}'.format(directory_sel))
    
    return path_object_files, path_sky_files, path_center_files, path_flux_files, path_sky_flux_files

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
            raise IOError('{0:s} is neither a cube nor a frame'.format(path_sel))
        data.append(data_read)
        header.append(header_read) 
        if silent is not True:
            printandlog('Read file ' + path_sel + '.')

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
                printandlog('Wrote file ' + path_sel + '.')
    else:
        # Write FITS-files with headers
        for data_sel, header_sel, path_sel in zip(data, header, path):
            hdu = pyfits.PrimaryHDU()
            hdu.data = data_sel.astype(np.float32)
            hdu.header = header_sel
            hdu.writeto(path_sel, overwrite=True, output_verify='silentfix+ignore')
            if silent is not True:
                printandlog('Wrote file ' + path_sel + '.')

###############################################################################
# remove_bad_pixels
###############################################################################

def remove_bad_pixels(cube, frame_master_bpm, sigmafiltering=True):
    ''' 
    Remove bad pixels from an image cube using the bad pixel map followed 
    by optional repeated sigma-filtering
    
    Input:
        cube: image data cube to filtered for bad pixels
        frame_master_bpm: frame indicating location of bad pixels with 0's and good
            pixels with 1's
        sigma_filtering: if True remove bad pixels remaining after applying
            master bad pixel map using repeated sigma-filtering (default = True)
    
    Output:
        cube_filtered: image data cube with bad pixels removed
    
    File written by Rob van Holstein
    Function status: verified
    '''

    # Define size of side of kernel for median filter
    filter_size_median = 5

    # Round filter size up to nearest odd number for a symmetric filter kernel
    filter_size_median = 2*(filter_size_median // 2) + 1
    
    # Remove bad pixels using the bad pixel map
    cube_median = ndimage.filters.median_filter(cube, size=(1, filter_size_median, \
                                                            filter_size_median))
    cube_filtered = cube_median + frame_master_bpm * (cube - cube_median)
    
    if sigmafiltering == True:
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
           
    return cube_filtered

###############################################################################
# process_sky_frames
###############################################################################

def process_sky_frames(path_sky_files, frame_master_bpm, sigmafiltering=True):
    '''
    Create a master sky-frame from the SKY-files
    
    Input:
        path_sky_files: string or list of strings specifying paths of SKY-files
        frame_master_bpm: frame indicating location of bad pixels with 0's and good
            pixels with 1's
        sigmafiltering: if True remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True)

    Output:
        frame_master_sky: master sky frame
    
    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: verified       
    '''
       
    # Read sky frames
    cube_sky_raw = read_fits_files(path=path_sky_files, silent=True)[0]
   
    if type(cube_sky_raw) == list:
        # Make a single image cube from list of image cubes or frames
        cube_sky_raw = np.vstack(cube_sky_raw)
  
    # Remove bad pixels of each frame
    cube_sky_filtered = remove_bad_pixels(cube=cube_sky_raw, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering)

    # Compute median of sky frames
    frame_master_sky = np.median(cube_sky_filtered, axis=0)
            
    printandlog('\nThe master sky was created out of {0:d} raw SKY-frame(s).\n'\
                .format(cube_sky_raw.shape[0]))
    
    return frame_master_sky

###############################################################################
# process_object_frames
###############################################################################

def process_object_frames(path_object_files, frame_master_flat, frame_master_bpm, frame_master_sky, offset_y, offset_x, star_x, star_y, sigmafiltering=True):
    '''
    Process the OBJECT frames by subtracting the background, flat-fielding, 
    removing bad pixels, computing the mean over the NDIT's, centering and
    computing the single sum and difference images
    
    Input:
        path_object_files: list of paths to raw OBJECT-files
        frame_master_flat: master flat frame
        frame_master_bpm: frame indicating location of bad pixels with 0's and good
            pixels with 1's
        frame_master_sky: master sky frame for OBJECT-files
# TODO: CENTERING VARIABLES HERE
        sigmafiltering: if True remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True)
       
    Output:
        cube_single_sum: cube of single-sum I_Q^+, I_Q^-, I_U^+ and I_U^- intensity images
        cube_single_difference: cube of single-difference Q^+, Q^-, U^+ and U^- images
        header: list of FITS-headers of raw science frames   
                
    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: verified
    '''
    
    # Create empty lists to store processed images and headers in
    list_single_sum = []
    list_single_difference = []
    header = []
    printandlog('')

    for i, path_sel in enumerate(path_object_files):
        # Read data and header from file
        cube_sel, header_sel = read_fits_files(path=path_sel, silent=True)

        # Filter master flat for zeros and NaN's
        frame_master_flat = np.nan_to_num(frame_master_flat)
        frame_master_flat[frame_master_flat == 0] = 1
            
        # Subtract background and divide by master flat
        cube_bgsubtr_flatfielded = (cube_sel - frame_master_sky) / frame_master_flat

        # Remove bad pixels of each frame
        cube_badpixel_filtered = remove_bad_pixels(cube=cube_bgsubtr_flatfielded, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering)

        # Compute mean over NDIT frames
        frame_mean = np.mean(cube_badpixel_filtered, axis=0)

        # Separate left and right part of image
        frame_left = frame_mean[:, :1024]
        frame_right = frame_mean[:, 1024:]

        # Retrieve dithering shifts in x- and y-direction from header
        x_dith = header_sel["HIERARCH ESO INS1 DITH POSX"]
        y_dith = header_sel["HIERARCH ESO INS1 DITH POSY"]
        
        #TODO: The centering routines should immediately output the shifts (but without the dithering offsets)
        # Compute shift for left and right images
        shift_x_left = 511.5 - star_x - x_dith
        shift_y_left = 511.5 - star_y - y_dith
        shift_x_right = shift_x_left + offset_x
        shift_y_right = shift_y_left + offset_y
        
        # Shift left and right images to center
        frame_left = ndimage.shift(frame_left, [shift_y_left, shift_x_left], order=3, mode='constant', cval=0.0, prefilter=True)   
        frame_right = ndimage.shift(frame_right, [shift_y_right, shift_x_right], order=3, mode='constant', cval=0.0, prefilter=True)   

        # Create single difference and sum image
        frame_single_sum = frame_left + frame_right
        frame_single_difference = frame_left - frame_right
        
        # Append single sum and difference images and header
        list_single_sum.append(frame_single_sum)
        list_single_difference.append(frame_single_difference)
        header.append(header_sel)
        
        # Print which file has been processed
        printandlog('Processed file ' + str(i + 1) + '/' + str(len(path_object_files)) + ': {0:s}'.format(path_sel))
    
    # Convert lists of single sum and difference images to image cubes
    cube_single_sum = np.stack(list_single_sum)
    cube_single_difference = np.stack(list_single_difference)    

    return cube_single_sum, cube_single_difference, header

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@
#@    Functions being worked on by Rob
#@   
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




###############################################################################
# process_flux_frames
###############################################################################

def process_flux_frames(path_flux_files, frame_master_flat, frame_master_bpm, frame_master_sky_flux, param_annulus_background, sigmafiltering=True):
    '''
    Create a master flux-frame from the FLUX-files
    
    Input:
        path_flux_files: list of paths to raw FLUX-files
        frame_master_flat: master flat frame
        frame_master_bpm: frame indicating location of bad pixels with 0's and good
            pixels with 1's
        frame_master_sky_flux: master sky frame for FLUX-files
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
        sigmafiltering: if True remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True)
       
    Output:
        frame_master_flux: master flux frame
                
    File written by Rob van Holstein
    Function status: 
    '''
      
    # Read flux frames
    cube_flux_raw = read_fits_files(path=path_flux_files, silent=True)[0]
   
    if type(cube_flux_raw) == list:
        # Make a single image cube from list of image cubes or frames
        cube_flux_raw = np.vstack(cube_flux_raw)
        
    # Filter master flat for zeros and NaN's
    frame_master_flat = np.nan_to_num(frame_master_flat)
    frame_master_flat[frame_master_flat == 0] = 1
        
    # Subtract background and divide by master flat
    cube_bgsubtr_flatfielded = (cube_flux_raw - frame_master_sky_flux) / frame_master_flat

    # Remove bad pixels of each frame
    cube_badpixel_filtered = remove_bad_pixels(cube=cube_bgsubtr_flatfielded, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering)

    # Separate left and right part of images
    cube_left = cube_badpixel_filtered[:, :, :1024]
    cube_right = cube_badpixel_filtered[:, :, 1024:]

    #TODO: The flux images should be centered using cross-correlation, similar to non-coro data
    # Compute shift for left and right images
    shift_x_left = 67.0
    shift_y_left = 21.0
    shift_x_right = 66.0
    shift_y_right = 32.0       
    
    # Shift left and right images to center
    cube_left_centered = np.zeros(cube_left.shape)
    cube_right_centered = np.zeros(cube_right.shape)
    
    for i, (frame_left, frame_right) in enumerate(zip(cube_left, cube_right)):
        cube_left_centered[i, :, :] = ndimage.shift(frame_left, [shift_y_left, shift_x_left], order=3, mode='constant', cval=0.0, prefilter=True)   
        cube_right_centered[i, :, :] = ndimage.shift(frame_right, [shift_y_right, shift_x_right], order=3, mode='constant', cval=0.0, prefilter=True)   
    
    # Compute mean over NDIT frames
    frame_left = np.mean(cube_left_centered, axis=0)
    frame_right = np.mean(cube_right_centered, axis=0)
    
    # Sum two images
    frame_flux = frame_left + frame_right

    # Determine background and subtract it
    frame_master_flux, background = subtract_background(cube=frame_flux, param_annulus_background=param_annulus_background)
    printandlog('\nSubtracted background in master flux image = ' + str(background))
       
    printandlog('\nThe master flux was created out of {0:d} raw FLUX-frame(s).\n'\
                .format(cube_flux_raw.shape[0]))

    return frame_master_flux



###############################################################################
# perform_preprocessing
###############################################################################

def perform_preprocessing(param_annulus_background_flux, sigmafiltering_sky=True, sigmafiltering_object=True, sigmafiltering_flux=True, save_preprocessed_data=False):
    '''
    
    Input:
        
        
    Output:
        
    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: 
    '''
    
    ###############################################################################
    # 
    ###############################################################################
    
    #TODO: Make master flats for Y and Ks.
    
    # Process paths
    static_flat = os.path.join(path_static_flat_badpixelmap, 'masterflat_H.fits')
    static_bpm = os.path.join(path_static_flat_badpixelmap, 'master_badpix.fits')

    # Loading Static Calibration Frames
    frame_master_flat = pyfits.getdata(static_flat)    
    frame_master_bpm = pyfits.getdata(static_bpm)

    ###############################################################################
    # Checking and sorting data and creating directories 
    ###############################################################################
    
    # Check and sort data and create directories
    printandlog('\n###############################################################################')
    printandlog('# Checking and sorting data and creating directories')
    printandlog('###############################################################################') 

    path_object_files, path_sky_files, path_center_files, path_flux_files, \
    path_sky_flux_files = check_sort_data_create_directories()

    ###############################################################################
    # Computing master sky for object images
    ###############################################################################
    
    if any(path_sky_files):
        # Process the sky files for the object files
        printandlog('\n###############################################################################')
        printandlog('# Processing SKY files for OBJECT-files')
        printandlog('###############################################################################') 

        frame_master_sky = process_sky_frames(path_sky_files=path_sky_files, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering_sky)
        
        # Write master sky-frame
        write_fits_files(data=frame_master_sky, path=os.path.join(path_sky_dir, name_file_root + 'master_sky.fits'), header=False, silent=False)
    
    else:
        # Create a master sky frame with only zeros
        frame_master_sky = np.zeros((1024, 2048))

    ###############################################################################
    # Reducing center files and extracting center coordinates
    ###############################################################################
    
    # Creating reduced center frames and obtaining center values
    printandlog('\n###############################################################################')
    printandlog('# Processing CENTER-files')
    printandlog('###############################################################################') 
    
    create_clean_center_data(path_center_files,frame_master_flat,frame_master_bpm,frame_master_sky)
    
    #TODO: I think it is best to remove the variable recompute, as we can now save
    # the single sum and single difference images to be used in post processing.
    # Centering is the most difficult part, so the one thing you want to re-run
    # is the centering so no need to save the center coordinates then.
    recompute=True
    offset_x, offset_y, star_x, star_y = do_centering(recompute=recompute,\
                                        centering_method = "Python native")
        
    ###############################################################################
    # Creating reduced and centered single-sum and -difference images
    ###############################################################################
    
    # Create reduced and centerd single-sum and -difference images
    printandlog('\n###############################################################################')
    printandlog('# Processing OBJECT-files')
    printandlog('###############################################################################') 

    cube_single_sum, cube_single_difference, header = process_object_frames(path_object_files=path_object_files, frame_master_flat=frame_master_flat, frame_master_bpm=frame_master_bpm, frame_master_sky=frame_master_sky, offset_y=offset_y, offset_x=offset_x, star_x=star_x, star_y=star_y, sigmafiltering=sigmafiltering_object)
       
    if save_preprocessed_data == True:
        # Write preprocessed cubes of single-sum and single-difference images 
        printandlog('\nSaving pre-processed data so that pre-processing can be skipped the next time.\n')
        write_fits_files(data=cube_single_sum, path=os.path.join(path_preprocessed_dir, 'cube_single_sum.fits'), header=False, silent=False)
        write_fits_files(data=cube_single_difference, path=os.path.join(path_preprocessed_dir, 'cube_single_difference.fits'), header=False, silent=False)

        # Write path of object files to a .txt-file to be able to read headers
        with open(os.path.join(path_preprocessed_dir, 'path_object_files.txt'), 'w') as fh:
            for path_sel in path_object_files:
                fh.write('%s\n' % path_sel)
        printandlog('Wrote file ' + os.path.join(path_preprocessed_dir, 'path_object_files.txt') + '.')


    ###############################################################################
    # Computing master sky for flux images
    ###############################################################################
    
    if any(path_sky_flux_files):
        # Process the sky files for the flux files
        printandlog('\n###############################################################################')
        printandlog('# Processing SKY-files for FLUX-files')
        printandlog('###############################################################################') 

        frame_master_sky_flux = process_sky_frames(path_sky_files=path_sky_flux_files, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering_sky)
        
        # Write master sky-frame
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
        if type(param_annulus_background_flux) == tuple or type(param_annulus_background_flux) == list:
            printandlog('\nThe background will be determined with a user-defined annulus or several \nuser-defined annuli:')
            if type(param_annulus_background_flux) == tuple:
                printandlog(param_annulus_background_flux)
            elif type(param_annulus_background_flux) == list:
                for x in param_annulus_background_flux:
                    printandlog(x)
        elif param_annulus_background_flux == 'large annulus':
            #TODO: When process_flux_frames is finished, check if background annulus is not too large
            param_annulus_background_flux = (511.5, 511.5, 320, 60, 0, 360)
            printandlog('\nThe background will be determined with a star-centered annulus located \nfar away from the star:')
            printandlog(param_annulus_background_flux)

        # Processthe flux files
        frame_master_flux = process_flux_frames(path_flux_files=path_flux_files, frame_master_flat=frame_master_flat, frame_master_bpm=frame_master_bpm, frame_master_sky_flux=frame_master_sky_flux, param_annulus_background=param_annulus_background_flux, sigmafiltering=sigmafiltering_flux)

        # Write master flux-frame
        write_fits_files(data=frame_master_flux, path=os.path.join(path_flux_dir, name_file_root + 'master_flux.fits'), header=False, silent=False)

#TODO: update function to process flux files to include centering by cross-correlation
#TODO: conversion of final images to mJansky/arcsec^2 should be part of post-processing part and optional. Also add possibility to express as contrast wrt central star?

    return cube_single_sum, cube_single_difference, header





















#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@
#@    Functions to reduce center frames and find center positions
#@   
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def do_centering(recompute=False,centering_method = "Python native",path_pipeline=None):
    
    switch = True
    
    if centering_method == "Python native":
        switch = False    
        offset_x, offset_y, star_x, star_y = IRDIS_centering_folder(recompute=recompute)
    
    if centering_method == "Linux SPHERE Pipeline":
        switch = False
        offset_x, offset_y, star_x, star_y = get_stellar_position_sphere_pipeline(path_pipeline)
        star_x = star_x - 0.5
        star_x = star_y - 0.5
    
    if centering_method == "Manual/external":
        switch = False
        offset_x = -1.41281338429
        offset_y = 10.6449086743
        star_x = 456.5
        star_y = 550.6 
    
    if switch:
        print("Warning, no valid centering method specified")
    
    return offset_x, offset_y, star_x, star_y

def create_clean_center_data(raw_center_images,frame_master_flat,frame_master_bpm,master_dark):
    """
    Input:
        - List of OBJECT,CENTER filenames.
    Output:
        -
    """
    print('Reading and reducing the OBJECT,CENTER frames')
    for f in raw_center_images:
        cube = pyfits.getdata(f)
        header = pyfits.getheader(f)
        if len(cube.shape) == 3:# this is a cube
            average_cube = np.median(cube,axis=0)
        elif len(cube.shape) == 2:# this is a frame
            average_cube = cube            
        if len(np.shape(master_dark)) == 3:    
            master_dark = np.squeeze(master_dark,axis=0)    

        average_cube = (average_cube - master_dark)/frame_master_flat
        average_cube = np.nan_to_num(average_cube)
        median = ndimage.filters.median_filter(average_cube, size=5)    
        masked = average_cube * frame_master_bpm
        median = median - median * frame_master_bpm
        average_cube = median + masked    
        new_filename = os.path.join(path_center_dir,\
                    os.path.basename(f).replace('.fits','_reduced.fits'))
        pyfits.writeto(new_filename,average_cube, header, \
                       overwrite=True, output_verify = "ignore")    
        print('Wrote {0:s}'.format(new_filename))
    return


def get_stellar_position_sphere_pipeline(path,dark,flat):
    
    path_center_dir = os.path.join(path,"center_frames") # make new folder for images
    subprocess.call('mkdir ' + path_center_dir, shell=True)    
    
    center_images = []

    for f in sorted(os.listdir(path)):

            # das muss sortiert werden da listdir unsortierten inhalt gibt lst = os.listdir(whatever_directory)
            #lst.sort()

            f = os.path.join(path,f)

            if os.path.isdir(f):
                continue                # skips all directories

            froot, ext = os.path.splitext(f)        #get the extension
            # (froot is d:\pictures\2008\2008_11_01_family\img0001)
            # (ext is fits)

            if ext.lower() != '.fits':
                continue                # skips all non fits files

            dir_root, fname = os.path.split(froot)  #get the dir and base names
            # (dir_root is d:\pictures\2008\2008_11_01_family, this is identical to path)
            # (fname is img0001)

            # we determine the size of the datacubes:

            raw_file = pyfits.open(f)
            header = raw_file[0].header
           
            if "HIERARCH ESO DPR TYPE" in header:
                if header["HIERARCH ESO DPR TYPE"] == 'OBJECT,CENTER':
                    center_images.append(f)
            
            raw_file.close()
    
    for k in range(0,len(center_images)):
        #print "mv " + center_images[k] + " " + path_center_dir
        subprocess.call("mv " + center_images[k] + " " + path_center_dir, shell=True)    
    
    center_images = []    
    
    for f in sorted(os.listdir(path_center_dir)):

            # das muss sortiert werden da listdir unsortierten inhalt gibt lst = os.listdir(whatever_directory)
            #lst.sort()

            f = os.path.join(path_center_dir,f)

            if os.path.isdir(f):
                continue                # skips all directories

            froot, ext = os.path.splitext(f)        #get the extension
            # (froot is d:\pictures\2008\2008_11_01_family\img0001)
            # (ext is fits)

            if ext.lower() != '.fits':
                continue                # skips all non fits files

            dir_root, fname = os.path.split(froot)  #get the dir and base names
            # (dir_root is d:\pictures\2008\2008_11_01_family, this is identical to path)
            # (fname is img0001)

            # we determine the size of the datacubes:
           
            center_images.append(f)

    drh_file = open(os.path.join(path_center_dir,"input_center"),"w")   
    
    for i in range(0,len(center_images)):
        drh_file.writelines(center_images[i] + " IRD_STAR_CENTER_WAFFLE_RAW\n")
    
    drh_file.close()    
    
    subprocess.call("esorex sph_ird_star_center input_center", cwd=path_center_dir ,shell=True)
    
    table_data = pyfits.open(os.path.join(path_center_dir,"star_center.fits"))
    table = table_data[1].data
    
    left_x = []
    left_y = []
    right_x = []
    right_y = []    
    
    for r in range(0,len(table)):
        left_x.append(table[r][1])
        left_y.append(table[r][2])
        right_x.append(table[r][3])
        right_y.append(table[r][4])
        
    left_x_average = np.average(np.array(left_x))
    left_y_average = np.average(np.array(left_y))
    right_x_average = np.average(np.array(right_x))    
    right_y_average = np.average(np.array(right_y))
    
    left_x_deviation = np.ptp(np.array(left_x))
    left_y_deviation = np.ptp(np.array(left_y))    
    right_x_deviation = np.ptp(np.array(right_x))    
    right_y_deviation = np.ptp(np.array(right_y))
    
    print(left_x_average - right_x_average) 
    print(left_y_average - right_y_average)        
    
    print(left_x_average, left_y_average)
    print(np.max([left_x_deviation,right_x_deviation]))    
    print(np.max([left_y_deviation,right_y_deviation]))

    print(left_x_deviation)  
    print(left_y_deviation)
    print(right_x_deviation)
    print(right_y_deviation)         
    
    return left_x_average - right_x_average, left_y_average - right_y_average, left_x_average, left_y_average


def fit_four_spots_array_input(scidata, spot_locations):

    x_values = []    
    y_values = []    
    
    sigma_x = 5
    sigma_y = 2     
    
    for k in range(0,len(spot_locations)):

        x_0 = spot_locations[k][0]
        y_0 = spot_locations[k][1]

        #assuming sorted input, upper left spot is 1, lower left spot is 4

        if k == 0:        
            theta = math.radians(135)
        if k == 1:        
            theta = math.radians(45)       
        if k == 2:        
            theta = math.radians(135)
        if k == 3:        
            theta = math.radians(45)       
        
        amplitude = np.max(scidata[int(y_0)-2:int(y_0)+2,int(x_0)-2:int(x_0)+2])
        fit_x, fit_y, fit_array = fit_elliptical_gaussian(scidata,amplitude,x_0,y_0,sigma_x,sigma_y,theta)  
        x_values.append(fit_x)
        y_values.append(fit_y)
        
    
    #print give_stellar_center(x_values, y_values)    
    
    return give_stellar_center(x_values, y_values)


def IRDIS_get_center_get_offset(scidata, header):

    #left:
    x_corse, y_corse = blob_detection_center_of_mass_array(scidata[0:1024,0:1024],refine_box=150)
    filtered_image = high_pass_header_dependent(scidata[0:1024,0:1024],header)
    spot_positions = give_spot_positions(x_corse,y_corse,header)
    refined_spot_positions = refine_spot_positions(filtered_image,spot_positions,box=10)
    left_center = fit_four_spots_array_input(filtered_image, refined_spot_positions)
    print("Centering left done")
    
    #right:    
    x_corse, y_corse = blob_detection_center_of_mass_array(scidata[0:1024,1024:2048],refine_box=150)
    filtered_image = high_pass_header_dependent(scidata[0:1024,1024:2048],header)
    spot_positions = give_spot_positions(x_corse,y_corse,header)
    refined_spot_positions = refine_spot_positions(filtered_image,spot_positions,box=10)
    right_center = fit_four_spots_array_input(filtered_image, refined_spot_positions)
    print("Centering right done")    
    
    offset_array = np.array(left_center) - np.array(right_center)    
    
    return offset_array[0], offset_array[1], left_center[0], left_center[1]


def IRDIS_centering_folder(recompute=False):
    """
    Compute the center of each frame. It saves the result in a csv file.
    Input:
        - recompute: boolean. If False, then if the centering has already 
        been computed, it only reads the result. If True, it recomputes the center.
    Output:
        - a dictionnary with the centers and offsets
    """
    
    if os.path.isfile(os.path.join(path_center_dir,'centers.csv')) and recompute ==False:
        pd_center = pd.read_csv(os.path.join(path_center_dir,'centers.csv'))
        print('The center file {0:s} already exists: the centering will not be recomputed'.format(os.path.join(path_center_dir,'centers.csv')))
    else:
        offset_x_array = []
        offset_y_array = []
        star_x_array = []
        star_y_array = []
        center_reduced_files = glob.glob(os.path.join(path_center_dir,'*_reduced.fits'))
        for f in center_reduced_files:
            print('Centering the frame {0:s}'.format(f))
            center_data = pyfits.getdata(f)
            header = pyfits.getheader(f)
            offset_x, offset_y, star_x, star_y = IRDIS_get_center_get_offset(center_data, header)
            offset_x_array.append(offset_x)
            offset_y_array.append(offset_y)
            star_x_array.append(star_x)
            star_y_array.append(star_y)
        dico_center = {'xoffset':offset_x_array,'yoffset':offset_y_array,\
                       'xstar':star_x_array,'ystar':star_y_array}
        pd_center = pd.DataFrame(dico_center)
        pd_center.to_csv(os.path.join(path_center_dir,'centers.csv'))

    x_err = np.std(np.array(pd_center['xstar']))
    y_err = np.std(np.array(pd_center['ystar']))
    star_x_final = np.median(np.array(pd_center['xstar']))
    star_y_final = np.median(np.array(pd_center['ystar']))
    offset_x_final = np.median(np.array(pd_center['xoffset']))
    offset_y_final = np.median(np.array(pd_center['yoffset']))
    
    if (x_err > 0.5) or (y_err > 0.5):
        print("Warning the deviation of stellar position is larger than 0.5 pixel between images!")
        print("x deviation: ", x_err)
        print("y deviation: ", y_err)
    print('The position of the central star was computed to be:')
    print('X={0:.2f} px +/- {1:.2f}'.format(star_x_final,x_err))
    print('Y={0:.2f} px +/- {1:.2f}'.format(star_y_final,y_err))    
    return offset_x_final, offset_y_final, star_x_final, star_y_final


def give_spot_positions(x_center,y_center,header):
    
    filter_band = header["HIERARCH ESO INS COMB IFLT"]

    if filter_band == "BB_Ks": sep = 64.5
    if filter_band == "DP_0_BB_K": sep = 64.5
    if filter_band == "DP_0_BB_H": sep = 48.5    
    if filter_band == "DP_0_BB_J": sep = 37.4
    if filter_band == "DP_0_BB_Y": sep = 31.0    
    
    angle = math.radians(45)    
        
    spot1_x = x_center - math.cos(angle) * sep 
    spot1_y = y_center + math.sin(angle) * sep
    
    spot2_x = x_center + math.cos(angle) * sep 
    spot2_y = y_center + math.sin(angle) * sep
    
    spot3_x = x_center + math.cos(angle) * sep 
    spot3_y = y_center - math.sin(angle) * sep
    
    spot4_x = x_center - math.cos(angle) * sep 
    spot4_y = y_center - math.sin(angle) * sep
    
    return [(spot1_x,spot1_y),(spot2_x,spot2_y),(spot3_x,spot3_y),(spot4_x,spot4_y)]


def refine_spot_positions(scidata,spot_positions,box=10):

    refined_array = []

    for r in range(0,len(spot_positions)):
        y = spot_positions[r][1]
        x = spot_positions[r][0]
        y_spot_refined,x_spot_refined = ndimage.measurements.center_of_mass(scidata[int(y-box):int(y+box),int(x-box):int(x+box)]-np.amin(scidata[int(y-box):int(y+box),int(x-box):int(x+box)]))   
        refined_array.append((x_spot_refined+(x-box),y_spot_refined+(y-box)))
        #print y_spot_refined,x_spot_refined
        #pyfits.writeto(os.path.join("/home/ginski/Desktop/test/python_sphere_cneter_test/","center_mass" + str(r) + ".fits"), scidata[y-box:y+box,x-box:x+box], output_verify="ignore")    
        #imshow(scidata[y-box:y+box,x-box:x+box])
        #waitforbuttonpress()

    return refined_array   


def blob_detection_center_of_mass_array(scidata,refine_box=150):
    
    y,x =  ndimage.measurements.center_of_mass(scidata)
    
    #print y,x    
    
    y_refined,x_refined = ndimage.measurements.center_of_mass(scidata[int(y-refine_box):int(y+refine_box),int(x-refine_box):int(x+refine_box)])   
    
    x_final = x_refined + (x-refine_box)     
    y_final = y_refined + (y-refine_box)    
    
    return x_final, y_final


def high_pass_header_dependent(image,header):
    
    filter_band = header["HIERARCH ESO INS COMB IFLT"]

    if filter_band == "BB_Ks": n = 10
    if filter_band == "DP_0_BB_K": n = 10
    if filter_band == "DP_0_BB_H": n = 6    
    if filter_band == "DP_0_BB_J": n = 5
    if filter_band == "DP_0_BB_Y": n = 4    
    
    
    blurred = gaussian_filter(image, sigma=n, truncate=6.1)    
    
    filtered_image = image - blurred
    
    return filtered_image

def give_stellar_center(x_values, y_values):
    
    
    m1 = (y_values[2] - y_values[0]) / (x_values[2] - x_values[0])
    m2 = (y_values[1] - y_values[3]) / (x_values[1] - x_values[3])    
    
    n1 = y_values[2] - m1*x_values[2]
    n2 = y_values[3] - m2*x_values[3]
    
    x_center = (n2 - n1) / (m1 - m2)
    y_center = m1*x_center + n1    
    
    return x_center, y_center

@custom_model
def elliptical_gaussian(peak,x_0,y_0,theta,sigma_x,sigma_y):
    
    a = np.cos(theta) * np.cos(theta) / (2 * sigma_x * sigma_x) + np.sin(theta) * np.sin(theta) / (2 * sigma_y * sigma_y)
    b = -1* np.sin(2 * theta) / (4 * sigma_x * sigma_x) + np.sin(2 * theta) / (4 * sigma_y * sigma_y)
    c = np.sin(theta) * np.sin(theta) / (2 * sigma_x * sigma_x) + np.cos(theta) * np.cos(theta) / (2 * sigma_y * sigma_y)    
    
    exponent = -1*(a*(x-x_0)*(x-x_0) - 2*b*(x-x_0)*(y-y_0) + c*(y-y_0)*(y-y_0))    
    gaussian = peak * np.exp(exponent)
    
    return gaussian
    

def fit_elliptical_gaussian(z,peak,x_0,y_0,sigma_x,sigma_y,theta):
    
    x_length = z.shape[1]
    center_x = x_length/2
    y_length = z.shape[0]    
    center_y = y_length/2    
    
    y, x = np.mgrid[:y_length, :x_length]
    
    #amplitude is the max count value of the PSF, x_0 and y_0 is the first guess, gamma is the fwhm of the moffat, alpha is something like the slope
    p_init = models.Gaussian2D(peak,x_0,y_0,sigma_x,sigma_y,theta)
    #p_init.amplitude.fixed = True
    p_init.theta.fixed = True
    fit_p = fitting.LevMarLSQFitter()
        
    with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        p = fit_p(p_init, x, y, z)
    
    #print fit_p.fit_info["cov_x"]
    #print p
    
    x_fit = p.x_mean[0]
    y_fit = p.y_mean[0]    
    
    return x_fit, y_fit, p(x,y)












#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@
#@    Post-processing functions
#@   
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###############################################################################
# compute_double_sum_double_difference
###############################################################################

def compute_double_sum_double_difference(cube_single_sum, cube_single_difference, header, double_difference_type='standard'):
    '''  
    Compute double-sum I_Q- and I_U-images and double-difference Q- and U-images
    
    Input:
        cube_single_sum: cube of single-sum I_Q^+, I_Q^-, I_U^+ and I_U^- intensity images
        cube_single_difference: cube of single-difference Q^+, Q^-, U^+ and U^- images
        header: list of FITS-headers of raw science frames    
        double_difference_type: type of double difference to be computed, either
        'standard' or 'normalized' (see van Holstein et al. 2019; default = 'standard')
               
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
       
    if double_difference_type == 'standard':
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
# check_header_matching
###############################################################################
  
def check_header_matching(header, header_name_to_check, print_if_matching=True):   
    '''
    Check matching of FITS-header values
    
    Input:
        header: headers or list of headers
        header_name_to_check: string of header name to check
        print_if_matching: if True print message when header matches; if False
            only print message if not matching (default = True)
    
    File written by Rob van Holstein
    Function status: verified
    '''    

    if type(header) == pyfits.header.Header:
        header = [header]

    header_to_check = [x[header_name_to_check] for x in header]

    if len(set(header_to_check)) == 1:
        if print_if_matching == True:
            printandlog('\nValues of header ' + header_name_to_check + ' match.')
    else:
        raise ValueError('Values of header ' + header_name_to_check + ' do not match.' 
                 '\nHeader values: ' + ', '.join([str(x) for x in header_to_check]))
 
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
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 deg is up and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 deg is up and rotating counterclockwise)
            
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
        coord_center_x, coord_center_y, inner_radius, width, start_angle, end_angle = param_sel
        outer_radius = inner_radius + width

        start_angle = np.mod(start_angle, 360)
        end_angle = np.mod(end_angle, 360)
        
        # Of each pixel calculate radius and angle in range [0, 360)
        radius = np.sqrt((xm - coord_center_x)**2 + (ym - coord_center_y)**2)
        angle = np.mod(np.rad2deg(np.arctan2(coord_center_x - xm, ym - coord_center_y)), 360)
           
        # Select pixels that satisfy provided requirements
        if start_angle < end_angle:
            coord_y, coord_x = np.nonzero(np.logical_and(np.logical_and(radius >= inner_radius, radius < outer_radius), np.logical_and(angle >= start_angle, angle < end_angle)))
        else:
            coord_y1, coord_x1 = np.nonzero(np.logical_and(np.logical_and(radius >= inner_radius, radius < outer_radius), np.logical_and(angle >= start_angle, angle < 360)))
            coord_y2, coord_x2 = np.nonzero(np.logical_and(np.logical_and(radius >= inner_radius, radius < outer_radius), np.logical_and(angle >= 0, angle < end_angle)))
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
# determine_star_polarization
###############################################################################

def determine_star_polarization(cube_I_Q, cube_I_U, cube_Q, cube_U, param_annulus_star, param_annulus_background):   
    '''
    Determine polarization of star in annulus
     
    Input:
        cube_I_Q: cube of I_Q-images
        cube_I_U: cube of I_U-images
        cube_Q: cube of Q-images
        cube_U: cube of U-images
        param_annulus_star: (list of) length-6-tuple(s) with parameters to generate annulus to measure polarization of star:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
              
    Output:
        q: normalized Stokes q measured in annulus
        u: normalized Stokes u measured in annulus
        DoLP: degree of linear polarization measured in annulus
        AoLP: angle of linear polarization measured in annulus
   
    File written by Rob van Holstein
    Function status: verified   
    '''

    # Set axis parameter of mean and median depending on whether the cube is a 3D or 2D object
    if cube_I_Q.ndim == 2:
        mean_median_axis = None
    elif cube_I_Q.ndim == 3:
        mean_median_axis = 1
           
    # Compute flux in I_Q, I_U, Q and U in an annulus minus the background in an annulus
    I_Q = np.mean(compute_annulus_values(cube=cube_I_Q, param=param_annulus_star)[0], axis=mean_median_axis) - \
                  np.median(compute_annulus_values(cube=cube_I_Q, param=param_annulus_background)[0], axis=mean_median_axis)       
    I_U = np.mean(compute_annulus_values(cube=cube_I_U, param=param_annulus_star)[0], axis=mean_median_axis) - \
                  np.median(compute_annulus_values(cube=cube_I_U, param=param_annulus_background)[0], axis=mean_median_axis)      
    Q = np.mean(compute_annulus_values(cube=cube_Q, param=param_annulus_star)[0], axis=mean_median_axis) - \
                np.median(compute_annulus_values(cube=cube_Q, param=param_annulus_background)[0], axis=mean_median_axis)       
    U = np.mean(compute_annulus_values(cube=cube_U, param=param_annulus_star)[0], axis=mean_median_axis) - \
                np.median(compute_annulus_values(cube=cube_U, param=param_annulus_background)[0], axis=mean_median_axis)    
    
    # Compute normalized Stokes q and u, DoLP and AoLP
    q = Q / I_Q
    u = U / I_U
    DoLP = np.sqrt(q**2 + u**2)
    AoLP = np.mod(np.rad2deg(0.5 * np.arctan2(u, q)), 180)

    return q, u, DoLP, AoLP

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
    printandlog('\nUsing model parameters corresponding to filter ' + filter_used[5:] + '.')
    delta_hwp = -0.613158589269
    delta_der = 0.500072483779
            
    if filter_used == 'FILT_BBF_Y':
        Delta_UT = 171.891576898
        Delta_M4 = 171.891576898
        epsilon_hwp = -0.00021492286258857182
        Delta_hwp = 184.24489040687035
        epsilon_der = -0.0009423930026144313                                     
        Delta_der = 126.11957538036766
        d_CI = 0.9801695109615369
                                     
    elif filter_used == 'FILT_BBF_J':
        Delta_UT = 173.414169049                
        Delta_M4 = 173.414169049
        epsilon_hwp = -0.00043278581895049085
        Delta_hwp = 177.52027378388257
        epsilon_der = -0.008303978181252019                                      
        Delta_der = 156.0584333408133
        d_CI = 0.9894796343284551
                                                                
    elif filter_used == 'FILT_BBF_H':
        Delta_UT = 174.998748608
        Delta_M4 = 174.998748608
        epsilon_hwp = -0.00029657803108325395
        Delta_hwp = 170.67214967707864
        epsilon_der = -0.002260131403393225                                       
        Delta_der = 99.32313652084311
        d_CI = 0.9955313968849352
                                                                           
    elif filter_used == 'FILT_BBF_Ks':
        Delta_UT = 176.302288996           
        Delta_M4 = 176.302288996
        epsilon_hwp = -0.00041456866069250524
        Delta_hwp = 177.61874393785442
        epsilon_der = 0.0035517563420643166                                        
        Delta_der = 84.13439892002613
        d_CI = 0.9841908773870153
    
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
        printandlog('\nUsing epsilon_UT and epsilon_M4 from before the recoating of M1 and M3 that took \nplace on April 16 2017.')
    elif not any(dates_before_recoating):
        printandlog('\nUsing epsilon_UT and epsilon_M4 from after the recoating of M1 and M3 that took \nplace on April 16 2017.')
    else:
        printandlog('\nUsing epsilon_UT and epsilon_M4 from both before and after the recoating of M1 \nand M3 that took place on April 16 2017.')

    # Convert observation dates to allow comparison to date of recoating
    dates = [time.strptime(x, '%Y-%m-%d') for x in dates]

    # Compute coefficient matrix describing the measurements
    n = len(p1) 
    X = np.zeros((n, 4))

    for i in range(n):
        
        # Define model parameters that depend on the date of the observations
        if filter_used == 'FILT_BBF_Y':
            if dates[i] < date_recoating:
                epsilon_UT = 0.023607413903534567 
                epsilon_M4 = 0.018211735858186456
            else:
                epsilon_UT = 0.01745394681183012
                epsilon_M4 = 0.018194769704367342
            
        elif filter_used == 'FILT_BBF_J':
            if dates[i] < date_recoating:
                epsilon_UT = 0.016685701811847004
                epsilon_M4 = 0.012844478639635984
            else:
                epsilon_UT = 0.01213513552053676
                epsilon_M4 = 0.013046513475544473
   
        elif filter_used == 'FILT_BBF_H':
            if dates[i] < date_recoating:
                epsilon_UT = 0.012930082215499676
                epsilon_M4 = 0.009845229155837451 
            else:
                epsilon_UT = 0.009032205030412622
                epsilon_M4 = 0.009220985704954044   
                
        elif filter_used == 'FILT_BBF_Ks':
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

def correct_instrumental_polarization_effects(cube_I_Q_double_sum, cube_I_U_double_sum, cube_Q_double_difference, cube_U_double_difference, header, path_reduced_dir, param_annulus_star, param_annulus_background, combination_method_polarization_images='trimmed mean', trimmed_mean_proportiontocut_polarization_images=0.1, combination_method_total_intensity_images='trimmed mean', trimmed_mean_proportiontocut_total_intensity_images=0.1, images_north_up=True):
    '''
    Calculate incident I_Q-, I_U-, Q- and U-images by correcting for the instrumental polarization effects of IRDIS using the polarimetric instrument model

    Input:
        cube_I_Q_double_sum: cube of double-sum intensity I_Q-images in order of HWP cycles
        cube_I_U_double_sum: cube of double-sum intensity I_U-images in order of HWP cycles  
        cube_Q_double_difference: cube of double-difference Stokes Q-images in order of HWP cycles
        cube_U_double_difference: cube of double-difference Stokes U-images in order of HWP cycles 
        header: list of FITS-headers of raw science frames    
        param_annulus_star: (list of) length-6-tuple(s) with parameters to generate annulus to measure polarization of star:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
        combination_method_polarization_images: method to be used to produce the incident Q- and U-images, 
            'least squares', 'trimmed mean' or 'median' (default = 'trimmed mean')
        trimmed_mean_proportiontocut_polarization_images: fraction to cut off of both tails of the distribution if 
            combination_method_polarization_images = 'trimmed mean' (default = 0.1) 
        combination_method_total_intensity_images: method to be used to produce the incident I_Q- and I_U-images, 
            'mean', 'trimmed mean' or 'median' (default = 'trimmed mean')
        trimmed_mean_proportiontocut_total_intensity_images: fraction to cut off of both tails of the distribution if 
            trimmed_mean_proportiontocut_total_intensity_images = 'trimmed mean' (default = 0.1) 
        images_north_up: if True the images produced are oriented with North up; if False the images have the image orientation of the
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

#TODO:    printandlog('\nExtracting (component) angles from headers.')

    # Determine filter used    
    check_header_matching(header=header, header_name_to_check='ESO INS1 FILT ID', print_if_matching=False)
    filter_used = header[0]['ESO INS1 FILT ID']   
    
    # Determine tracking mode used
    check_header_matching(header=header, header_name_to_check='ESO INS4 COMB ROT', print_if_matching=False)
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

    # Define mean solar day (s)
    msd = 86400
    
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
        printandlog('\nInterpolating the altitude angle using a single spline for this single night of \nobservations.' )
    else:
        printandlog('\nInterpolating the altitude angles using a separate spline for each of the ' + str(number_of_nights) + ' \nnights.')
    
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
    mjd_half1 = mjd_half[indices_QUplus]
    mjd_half2 = mjd_half[indices_QUminus]
    
    # Plot header values as a function of time
    plot_name = name_file_root + 'header_angles.png'
    printandlog('\nCreating plot ' + plot_name + ' showing the parallactic, \naltitude, HWP and derotator angles of the observations.')
    angles = np.hstack([p1, p2, a1, a2, theta_hwp1, theta_hwp2, theta_der1, theta_der2])
    font_size = 10
    plt.figure(figsize = (12, 8))
    plt.plot(mjd_half1, p1, 'ob', label = 'Parallactic angle 1')    
    plt.plot(mjd_half2, p2, 'sb', label = 'Parallactic angle 2')    
    plt.plot(mjd_half1, a1, 'or', label = 'Altitude angle 1')       
    plt.plot(mjd_half2, a2, 'sr', label = 'Altitude angle 2')       
    plt.plot(mjd_half1, theta_hwp1, 'og', label = 'HWP angle 1')       
    plt.plot(mjd_half2, theta_hwp2, 'sg', label = 'HWP angle 2')    
    plt.plot(mjd_half1, theta_der1, 'ok', label = 'Derotator angle 1') 
    plt.plot(mjd_half2, theta_der2, 'sk', label = 'Derotator angle 2')
    ax = plt.gca()
    ax.set_xlabel(r'Modified Julien Date', fontsize = font_size)
    ax.tick_params(axis = 'x', labelsize = font_size)
    ax.set_ylim([np.min(angles) - 10, np.max(angles) + 10])
    ax.set_ylabel(r'Angle ($^\circ$)', fontsize = font_size)
    ax.tick_params(axis = 'y', labelsize = font_size)  
    ax.grid()
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.savefig(os.path.join(path_reduced_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.show()     
    
    ###############################################################################
    # Compute model coefficient matrix, IP and cross-talk elements
    ###############################################################################

    # Construct coefficient matrix describing the measurements
#TODO:    printandlog('\nComputing the coefficient matrix describing the measurements.')
    X = compute_irdis_model_coefficient_matrix(p1=p1, p2=p2, a1=a1, a2=a2, theta_hwp1=theta_hwp1, theta_hwp2=theta_hwp2, theta_der1=theta_der1, theta_der2=theta_der2, dates=dates, filter_used=filter_used) 
 
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

    if tracking_mode_used == 'FIELD' and number_derotator_position_angles == 1 and images_north_up == False:
        # Set rotation angles of images equal to zero
        printandlog('\nSetting the image rotation angles for field-tracking mode equal to zero so \nthat North is not up in the final images.')
        rotation_angles_Q = np.zeros(len(indices_Qplus))
        rotation_angles_U = np.zeros(len(indices_Uplus))

    elif tracking_mode_used == 'FIELD':
        # Compute rotation angle of images
        printandlog('\nComputing the image rotation angles for field-tracking mode to orient the \nimages with North up. The true North correction used equals %.2f deg.' % true_north_correction)
        rotation_angles_Q = -derotator_position_angle[indices_Qplus] - true_north_correction
        rotation_angles_U = -derotator_position_angle[indices_Uplus] - true_north_correction
        
    elif tracking_mode_used == 'PUPIL':    
        # Compute mean parallactic angles of first and second measurements used to compute I_Q and Q, and I_U and U
        printandlog('\nComputing the image rotation angles for pupil-tracking mode to orient the \nimages with North up. The true North correction used equals %.2f deg.' % true_north_correction)
        p_Q = compute_mean_angle(np.vstack([p[indices_Qplus], p[indices_Qminus]]), axis = 0)
        p_U = compute_mean_angle(np.vstack([p[indices_Uplus], p[indices_Uminus]]), axis = 0)
        
        # Compute rotation angle of images
        rotation_angles_Q = -p_Q - pupil_offset - true_north_correction
        rotation_angles_U = -p_U - pupil_offset - true_north_correction

    ###############################################################################
    # Measure polarization signal in annulus on star
    ###############################################################################


    # If the annuli for the star and background are centered on the star and rotationally symmetric
    if np.all(np.array(param_annulus_star, ndmin=2)[:, np.array([0, 1, 4, 5])] == np.array([511.5, 511.5, 0, 360])) and \
       np.all(np.array(param_annulus_background, ndmin=2)[:, np.array([0, 1, 4, 5])] == np.array([511.5, 511.5, 0, 360])):
        # Do not rotate the cubes of double-sum and double-difference images and compute normalized Stokes q and u in an annulus   
#TODO:        printandlog('\nDetermining the polarization signal in an annulus before correcting for the \ninstrumental polarization effects.')
        printandlog('\nNot rotating the images used to determine the polarization signal in an \nannulus because the annuli used for the star and the background are \ncentered on the star and rotationally symmetric.')

        q_annulus, u_annulus = determine_star_polarization(cube_I_Q=cube_I_Q_double_sum, cube_I_U=cube_I_U_double_sum, cube_Q=cube_Q_double_difference, cube_U=cube_U_double_difference, param_annulus_star=param_annulus_star, param_annulus_background=param_annulus_background)[:2]

    else:
        # Rotate the cubes of double-sum and double-difference images so that the annuli are at the right position
        printandlog('\nRotating the images used to determine the polarization signal in annulus \nbecause the annuli used for the star and the background are not centered \non the star and/or not rotationally symmetric.')
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
#TODO:        printandlog('\nDetermining the polarization signal in an annulus before correcting for the \ninstrumental polarization effects.')
        q_annulus, u_annulus = determine_star_polarization(cube_I_Q=cube_I_Q_annulus, cube_I_U=cube_I_U_annulus, cube_Q=cube_Q_annulus, cube_U=cube_U_annulus, param_annulus_star=param_annulus_star, param_annulus_background=param_annulus_background)[:2]
    
    ###############################################################################
    # Fit polarization of star from measurements using model coefficient matrix
    ###############################################################################
      
#TODO:    printandlog('\nFitting the star polarization from the measured polarization signal using the \nmodel coefficient matrix:')
    
    # Create y-data array                               
    ydata = np.zeros(len(np.hstack([q_annulus, u_annulus])))
    ydata[indices_Q] = q_annulus
    ydata[indices_U] = u_annulus
  
    # Define function that computes the sum of squared residuals of fit
    def residuals(parameters, ydata, X):
        return np.sum((ydata - X.dot(np.array([1, parameters[0], parameters[1], 0])))**2) 

    # Fit data to model
    fit_result = optimize.minimize(residuals, np.array([0, 0]), args=(ydata, X), method='SLSQP', bounds=[(-1,1), (-1,1)], tol=1e-14, options={'maxiter':100, 'disp':False})
 
    # Retrieve normalized Stokes q and u and compute DoLP and AoLP of star
    q_star, u_star = fit_result.x
    DoLP_star = np.sqrt(q_star**2 + u_star**2)
    AoLP_star = np.mod(np.rad2deg(0.5 * np.arctan2(u_star, q_star)), 180)
    
    # Print resulting star polarization
    printandlog('\nFitted star polarization:')
    printandlog('q_star = %.4f %%' % (100*q_star))
    printandlog('u_star = %.4f %%' % (100*u_star))
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
    printandlog('\nCreating plot ' + plot_name + ' showing the model-predicted IP, \nmeasured polarization signal and model-predicted IP + fitted star polarization \nvs. HWP cycle number.')
    font_size = 10  
    marker_size = 4
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[5, 2]}, sharex=True, figsize = (9.4, 10.34))  
    ax1.plot([0, len(IP_Q) + 1],[0, 0], '-k', zorder=1) 
    ph1, = ax1.plot(np.arange(1, len(IP_Q) + 1), 100*IP_Q, '--b', markersize = marker_size, label = '$q$ model', zorder=5)  
    ph2, = ax1.plot(np.arange(1, len(IP_U) + 1), 100*IP_U, '--r', markersize = marker_size, label = '$u$ model', zorder=4)
    ph3, = ax1.plot(np.arange(1, len(IP_Q) + 1), 100*q_annulus, 'ob', markersize = marker_size, label = '$q$ meas.', zorder=9)  
    ph4, = ax1.plot(np.arange(1, len(IP_U) + 1), 100*u_annulus, 'or', markersize = marker_size, label = '$u$ meas.', zorder=8) 
    ph5, = ax1.plot(np.arange(1, len(IP_Q) + 1), 100*q_fit, '-b', markersize = marker_size, label = '$q$ model + star fit.', zorder=7)  
    ph6, = ax1.plot(np.arange(1, len(IP_U) + 1), 100*u_fit, '-r', markersize = marker_size, label = '$u$ model + star fit.', zorder=6)  
    ax1.set_ylabel(r'Normalized Stokes parameter (%)', fontsize = font_size)
    ax1.tick_params(axis = 'y', labelsize = font_size)
    ax1.grid()
    l1 = ax1.legend([ph1, ph3, ph5], ['$q$ model', '$q$ meas.', '$q$ model + $q/u_\mathrm{in}$ fit.'], loc = 'center left', prop={'size': font_size}, handlelength = 2.5, ncol = 1)
    ax1.legend([ph2, ph4, ph6], ['$u$ model', '$u$ meas.', '$u$ model + $q/u_\mathrm{in}$ fit.'], loc = 'lower right', prop={'size': font_size}, handlelength = 2.5, ncol = 1)
    ax1.add_artist(l1)
    ax2.plot([0, len(IP_Q) + 1],[0, 0], '-k', zorder=1) 
    ax2.plot(np.arange(1, len(IP_Q) + 1), 100*res_Q, 'ob', markersize = marker_size, label = '$q$ model', zorder=5)  
    ax2.plot(np.arange(1, len(IP_U) + 1), 100*res_U, 'or', markersize = marker_size, label = '$u$ model', zorder=4)
    ax2.set_xlabel(r'HWP cycle', fontsize = font_size)
    ax2.tick_params(axis = 'x', labelsize = font_size)
    ax2.set_xlim([0, len(IP_Q) + 1])
    ax2.set_ylabel(r'Residuals (%)', fontsize = font_size)
    ax2.tick_params(axis = 'y', labelsize = font_size)
    ax2.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(path_reduced_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.show()

    # Plot elements Q->Q and U->Q from model as a function of HWP cycle number
    plot_name = name_file_root + 'model_crosstalk_transmission.png'
    printandlog('\nCreating plot ' + plot_name + ' showing the \nmodel-predicted polarized transmission and crosstalk/rotation elements vs. HWP \ncycle number.')
    font_size = 10
    plt.figure(figsize = (5.9, 3.8))
    plt.plot([0, len(QQ_Q) + 1],[0, 0], '-k')  
    plt.plot(np.arange(1, len(QQ_Q) + 1), QQ_Q, 'o-b', label = r'Qin $\rightarrow$ Q')  
    plt.plot(np.arange(1, len(UQ_Q) + 1), UQ_Q, 'o-', color = (0.5, 0, 0.5), label = r'Uin $\rightarrow$ Q')
    plt.plot(np.arange(1, len(QQ_U) + 1), QQ_U, 'o-', color = (1, 0.5, 0), label = r'Qin $\rightarrow$ U')
    plt.plot(np.arange(1, len(UQ_U) + 1), UQ_U, 'o-r', label = r'Uin $\rightarrow$ U')  
    ax = plt.gca()
    ax.set_xlabel(r'HWP cycle', fontsize = font_size)
    ax.tick_params(axis = 'x', labelsize = font_size)
    ax.set_xlim([0, len(QQ_Q) + 1])
    ax.set_ylabel(r'Retardance/rotation elements from model', fontsize = font_size)
    ax.tick_params(axis = 'y', labelsize = font_size)  
    ax.grid()
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.savefig(os.path.join(path_reduced_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name), dpi = 300, bbox_inches = 'tight')
    plt.show()    
      
    ###############################################################################
    # Compute incident Q- and U-images by correcting for instrumental polarization effects
    ###############################################################################

    # Subtract intrumental polarization from Q- and U-images
#TODO:    printandlog('\nCorrecting the Q- and U-images for the instrumental polarization effects and \nrotating them to the desired final orientation.')
    cube_Q_IP_subtracted = cube_Q_double_difference - IP_Q[:, np.newaxis, np.newaxis]*cube_I_Q_double_sum
    cube_U_IP_subtracted = cube_U_double_difference - IP_U[:, np.newaxis, np.newaxis]*cube_I_U_double_sum
    
    # Derotate IP-subtracted Q- and U-images
    cube_Q_derotated = np.zeros(cube_Q_IP_subtracted.shape)
    cube_U_derotated = np.zeros(cube_U_IP_subtracted.shape)

    for i, (frame_Q, rotation_angle_Q) in enumerate(zip(cube_Q_IP_subtracted, rotation_angles_Q)):
        cube_Q_derotated[i, :, :] = rotate(frame_Q, rotation_angle_Q, reshape=False)
    
    for i, (frame_U, rotation_angle_U) in enumerate(zip(cube_U_IP_subtracted, rotation_angles_U)):
        cube_U_derotated[i, :, :] = rotate(frame_U, rotation_angle_U, reshape=False)
         
    # Calculate Q- and U-images incident on telescope by solving system of equations per HWP cycle
    cube_Q_incident = np.zeros(cube_Q_derotated.shape)
    cube_U_incident = np.zeros(cube_U_derotated.shape)
    
    for i, (frame_Q, frame_U) in enumerate(zip(cube_Q_derotated, cube_U_derotated)):
        # Apply model correction of Q->Q and U->Q elements by solving equations per HWP cycle
        Y = np.stack([frame_Q, frame_U])    
        Y_stretched = Y.reshape(Y.shape[0], Y.shape[1]*Y.shape[2])
        X = X_QU[2*i:2*i + 2, :]
        cube_QU_incident_stretched = np.linalg.solve(X, Y_stretched)
        cube_QU_incident = cube_QU_incident_stretched.reshape(cube_QU_incident_stretched.shape[0], Y.shape[1], Y.shape[2]) 
        cube_Q_incident[i, :, :] = cube_QU_incident[0, :, :]
        cube_U_incident[i, :, :] = cube_QU_incident[1, :, :]
        
    # Create incident Q- and U-images
    if combination_method_polarization_images == 'least squares':
        # Obtain incident Q- and U-images from the least squares solution
        printandlog('\nComputing the incident Q- and U-images using least squares.')
          
        # Create cube Y with Q/U measurements in order of HWP cycles
        Y = np.zeros(((len(p1),) + cube_Q_derotated.shape[-2:]))
        Y[indices_Q, :, :] = cube_Q_derotated
        Y[indices_U, :, :] = cube_U_derotated
    
        # Make images 1D, compute linear least-squares solution and reshape final images to 2D again   
        Y_stretched = Y.reshape(Y.shape[0], Y.shape[1]*Y.shape[2])     
        QU_in_stretched = np.linalg.lstsq(X_QU, Y_stretched)[0]
        QU_in = QU_in_stretched.reshape(QU_in_stretched.shape[0], Y.shape[1], Y.shape[2])    
        frame_Q_incident = QU_in[0, :, :] 
        frame_U_incident = QU_in[1, :, :]
    
    elif combination_method_polarization_images == 'trimmed mean': 
        # Compute incident Q- and U-images from the trimmed mean of incident cubes
        printandlog('\nComputing the incident Q- and U-images using the trimmed mean with a proportion \nto cut equal to ' + str(trimmed_mean_proportiontocut_polarization_images) + '.')
        frame_Q_incident = trim_mean(cube_Q_incident, proportiontocut=trimmed_mean_proportiontocut_polarization_images, axis=0)
        frame_U_incident = trim_mean(cube_U_incident, proportiontocut=trimmed_mean_proportiontocut_polarization_images, axis=0)
    
    elif combination_method_polarization_images == 'median': 
        # Compute incident Q- and U-images from the median of incident cubes 
        printandlog('\nComputing the incident Q- and U-images using the median.')
        frame_Q_incident = np.median(cube_Q_incident, axis=0)
        frame_U_incident = np.median(cube_U_incident, axis=0)

    ###############################################################################
    # Compute incident I_Q- and I_U-images
    ###############################################################################

    # Derotate I_Q- and I_U-images
#TODO:    printandlog('\nRotating the I_Q- and I_U-images to the desired final orientation.')
    cube_I_Q_incident = np.zeros(cube_I_Q_double_sum.shape)
    cube_I_U_incident = np.zeros(cube_I_U_double_sum.shape)
    
    for i, (frame_I_Q, rotation_angle_Q) in enumerate(zip(cube_I_Q_double_sum, rotation_angles_Q)):
        cube_I_Q_incident[i, :, :] = rotate(frame_I_Q, rotation_angle_Q, reshape=False)
    
    for i, (frame_I_U, rotation_angle_U) in enumerate(zip(cube_I_U_double_sum, rotation_angles_U)):
        cube_I_U_incident[i, :, :] = rotate(frame_I_U, rotation_angle_U, reshape=False)

    # Create incident I_Q- and I_U-images
    if combination_method_total_intensity_images == 'mean':
        # Compute incident I_Q- and I_U-images from the mean 
        printandlog('\nComputing the incident I_Q- and I_U-images using the mean.')
        frame_I_Q_incident = np.mean(cube_I_Q_incident, axis = 0)
        frame_I_U_incident = np.mean(cube_I_U_incident, axis = 0)  
    
    elif combination_method_total_intensity_images == 'trimmed mean': 
        # Compute incident I_Q- and I_U-images from the trimmed mean
        printandlog('\nComputing the incident I_Q- and I_U-images using the trimmed mean with a proportion \nto cut equal to ' + str(trimmed_mean_proportiontocut_total_intensity_images) + '.')
        frame_I_Q_incident = trim_mean(cube_I_Q_incident, proportiontocut=trimmed_mean_proportiontocut_total_intensity_images, axis=0)
        frame_I_U_incident = trim_mean(cube_I_U_incident, proportiontocut=trimmed_mean_proportiontocut_total_intensity_images, axis=0)
    
    elif combination_method_total_intensity_images == 'median': 
        # Compute incident I_Q- and I_U-images from the median
        printandlog('\nComputing the incident I_Q- and I_U-images using the median.')
        frame_I_Q_incident = np.median(cube_I_Q_incident, axis=0)
        frame_I_U_incident = np.median(cube_I_U_incident, axis=0)
  
    return frame_I_Q_incident, frame_I_U_incident, frame_Q_incident, frame_U_incident, cube_I_Q_incident, cube_I_U_incident, cube_Q_incident, cube_U_incident

###############################################################################
# subtract_background
###############################################################################

def subtract_background(cube, param_annulus_background):   
    '''
    Subtract background from cube or frame
     
    Input:
        cube: image cube or frame to subtract background from
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)

    Output:
        cube_background_subtracted: image cube or frame with background subtracted
        background: float or array of background values for each image
   
    File written by Rob van Holstein
    Function status: verified   
    '''

    # Determine background in frame or cube and subtract it
    if cube.ndim == 2:
        background = np.median(compute_annulus_values(cube=cube, param=param_annulus_background)[0])
        cube_background_subtracted = cube - background
    elif cube.ndim == 3:
        background = np.median(compute_annulus_values(cube=cube, param=param_annulus_background)[0], axis=1)
        cube_background_subtracted = cube - background[:, np.newaxis, np.newaxis] 
   
    return cube_background_subtracted, background

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
    
def compute_final_images(frame_I_Q, frame_I_U, frame_Q, frame_U, header, images_north_up):
    ''' 
    Compute final images of polarimetric data reduction
    
    Input:
        frame_I_Q: I_Q-image
        frame_I_U: I_U-image
        frame_Q: Q-image
        frame_U: U-image
        header: list of headers of raw science frames    
        images_north_up: if True the images produced are oriented with North up; if False the images have the image orientation of the
            raw frames (default = True); only valid for observations taken in field-tracking mode with a single derotator 
            position angle; parameter is ignored for pupil-tracking observations or field-tracking observations with more 
            than one derotator position angle, because in these cases the final images produced always have North up
        
    Output:
        frame_I_tot: total-intensity image
        frame_Q_phi: image of Q_phi
        frame_U_phi: image of U_phi
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
    
    # Determine tracking mode used
    tracking_mode_used = header[0]['ESO INS4 COMB ROT']
    
    # Determine number of unique derotator position angles
    derotator_position_angle = [x['ESO INS4 DROT2 POSANG'] for x in header]
    number_derotator_position_angles = len(np.unique(derotator_position_angle))

    if tracking_mode_used == 'FIELD' and number_derotator_position_angles == 1 and images_north_up == False:
        # Compute image position angle
        rotation_angle = -derotator_position_angle[0] - true_north_correction
    else:
        rotation_angle = 0.0
    
    frame_Q_phi, frame_U_phi, frame_azimuthal_angle = compute_azimuthal_stokes_parameters(frame_Q, frame_U, rotation_angle=-rotation_angle)

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

def perform_postprocessing(cube_single_sum, cube_single_difference, header, param_annulus_star, param_annulus_background, double_difference_type='standard', remove_vertical_band_detector_artefact=True, combination_method_polarization_images='trimmed mean', trimmed_mean_proportiontocut_polarization_images=0.1, combination_method_total_intensity_images='trimmed mean', trimmed_mean_proportiontocut_total_intensity_images=0.1, images_north_up=True, create_images_DoLP_AoLP_q_u_norm=False):
    '''
    Perform post-processing of data and save final images to FITS-files
    
    Input:
        cube_single_sum: cube of pre-processed single-sum images
        cube_single_difference: cube of pre-processed single-difference images
        header: list of FITS-headers of raw science frames    
        param_annulus_star: (list of) length-6-tuple(s) with parameters to generate annulus to measure polarization of star:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
        double_difference_type: type of double difference to be computed, either
        'standard' or 'normalized' (see van Holstein et al. 2019; default = 'standard')
        remove_vertical_band_detector_artefact: If True remove the vertical band detector artefact seen in 
            the double-difference Q- and U-images. If False don't remove it.
        combination_method_polarization_images: method to be used to produce the incident Q- and U-images, 
            'least squares', 'trimmed mean' or 'median' (default = 'trimmed mean')
        trimmed_mean_proportiontocut_polarization_images: fraction to cut off of both tails of the distribution if 
            combination_method_polarization_images = 'trimmed mean' (default = 0.1) 
        combination_method_total_intensity_images: method to be used to produce the incident I_Q- and I_U-images, 
            'mean', 'trimmed mean' or 'median' (default = 'trimmed mean')
        trimmed_mean_proportiontocut_total_intensity_images: fraction to cut off of both tails of the distribution if 
            trimmed_mean_proportiontocut_total_intensity_images = 'trimmed mean' (default = 0.1) 
        images_north_up: if True the images produced are oriented with North up; if False the images have the image orientation of the
            raw frames (default = True); only valid for observations taken in field-tracking mode with a single derotator 
            position angle; parameter is ignored for pupil-tracking observations or field-tracking observations with more 
            than one derotator position angle, because in these cases the final images produced always have North up
        create_images_DoLP_AoLP_q_u_norm: if True create final images of degree of linear polarization, normalized Stokes q and u
            and degree and angle of linear polarization computed from q and u; such images only have meaning if all flux in the images
            originates from the astrophysical source of interest (default = False)
        
    File written by Rob van Holstein
    Function status: verified      
    '''

    ###############################################################################
    # Define annulus for star polarization and background
    ###############################################################################
    
    printandlog('\n###############################################################################')
    printandlog('# Defining the annulus for star polarization and background')
    printandlog('###############################################################################') 
                
    # Determine filter used    
    check_header_matching(header=header, header_name_to_check='ESO INS1 FILT ID', print_if_matching=False)
    filter_used = header[0]['ESO INS1 FILT ID']   

    # Define and print annulus to determine the star polarization from
    if type(param_annulus_star) == tuple or type(param_annulus_star) == list:
        printandlog('\nThe star polarization will be determined with a user-defined annulus or \nseveral user-defined annuli:')
        if type(param_annulus_star) == tuple:
            printandlog(param_annulus_star)
        elif type(param_annulus_star) == list:
            for x in param_annulus_star:
                printandlog(x)
    elif param_annulus_star == 'ao residuals':
        if filter_used == 'FILT_BBF_Y':
            param_annulus_star = (511.5, 511.5, 40, 25, 0, 360)
        elif filter_used == 'FILT_BBF_J':
            param_annulus_star = (511.5, 511.5, 45, 30, 0, 360)
        elif filter_used == 'FILT_BBF_H':
            param_annulus_star = (511.5, 511.5, 60, 35, 0, 360)
        elif filter_used == 'FILT_BBF_Ks':
            param_annulus_star = (511.5, 511.5, 80, 40, 0, 360)  
        printandlog('\nThe star polarization will be determined with a star-centered annulus \nlocated over the AO residuals:')
        printandlog(param_annulus_star)
    elif param_annulus_star == 'star aperture':
        param_annulus_star = (511.5, 511.5, 0, 11, 0, 360)
        printandlog('\nThe star polarization will be determined with an aparture located at \nthe position of the central star:')
        printandlog(param_annulus_star)
    
    # Define and print annulus to determine the background from
    if type(param_annulus_background) == tuple or type(param_annulus_background) == list:
        printandlog('\nThe background will be determined with a user-defined annulus or several \nuser-defined annuli:')
        if type(param_annulus_background) == tuple:
            printandlog(param_annulus_background)
        elif type(param_annulus_background) == list:
            for x in param_annulus_background:
                printandlog(x)
    elif param_annulus_background == 'large annulus':
        param_annulus_background = (511.5, 511.5, 360, 60, 0, 360)
        printandlog('\nThe background will be determined with a star-centered annulus located \nfar away from the central star:')
        printandlog(param_annulus_background)
    
    ###############################################################################
    # Compute double sum and difference
    ###############################################################################
    
    printandlog('\n###############################################################################')
    printandlog('# Computing the double sum and double difference')
    printandlog('###############################################################################') 
                
    cube_I_Q_double_sum, cube_I_U_double_sum, cube_Q_double_difference, cube_U_double_difference \
    = compute_double_sum_double_difference(cube_single_sum, cube_single_difference, header, double_difference_type=double_difference_type)
    
    ###############################################################################
    # Remove vertical band detector artefact 
    ###############################################################################
    
    if remove_vertical_band_detector_artefact == True:
        # Remove IRDIS vertical band detector artefact seen in double-difference Q- and U-images
        number_pixels_artefact = 60
        printandlog('\nRemoving the vertical band detector artefact in the double-difference Q- and \nU-images by subtracting the median of the top and bottom ' + str(number_pixels_artefact) + ' pixels of each \ncolumn.')
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
    = correct_instrumental_polarization_effects(cube_I_Q_double_sum=cube_I_Q_double_sum, cube_I_U_double_sum=cube_I_U_double_sum, cube_Q_double_difference=cube_Q_artefact_removed, cube_U_double_difference=cube_U_artefact_removed, header=header, path_reduced_dir=path_reduced_dir, param_annulus_star=param_annulus_star, param_annulus_background=param_annulus_background, combination_method_polarization_images=combination_method_polarization_images, trimmed_mean_proportiontocut_polarization_images=trimmed_mean_proportiontocut_polarization_images, combination_method_total_intensity_images=combination_method_total_intensity_images, trimmed_mean_proportiontocut_total_intensity_images=trimmed_mean_proportiontocut_total_intensity_images, images_north_up=images_north_up)
    
    ###############################################################################
    # Subtract background in I_Q-, I_U-, Q- and U-frames
    ###############################################################################
 
    printandlog('\n###############################################################################')
    printandlog('# Subtracting the backgrounds and measuring and removing the star polarization')
    printandlog('###############################################################################') 

    # Determine background in corrected I_Q-, I_U-, Q- and U-frames and subtract it
#TODO:    printandlog('\nSubtracting the backgrounds in the incident I_Q-, I_U-, Q- and U-images:')
    frame_I_Q_background_subtracted, background_frame_I_Q = subtract_background(cube=frame_I_Q_incident, param_annulus_background=param_annulus_background)
    frame_I_U_background_subtracted, background_frame_I_U = subtract_background(cube=frame_I_U_incident, param_annulus_background=param_annulus_background)
    frame_Q_background_subtracted, background_frame_Q = subtract_background(cube=frame_Q_incident, param_annulus_background=param_annulus_background)
    frame_U_background_subtracted, background_frame_U = subtract_background(cube=frame_U_incident, param_annulus_background=param_annulus_background)
   
    # Print resulting background values
    printandlog('\nSubtracted backgrounds in the incident I_Q-, I_U-, Q- and U-images:')
    printandlog('Background I_Q = ' + str(background_frame_I_Q))
    printandlog('Background I_U = ' + str(background_frame_I_U))
    printandlog('Background Q = ' + str(background_frame_Q))
    printandlog('Background U = ' + str(background_frame_U))
                 
    ###############################################################################
    # Determine star polarization in frames and subtract it 
    ###############################################################################    
    
    # Compute normalized Stokes q and u, DoLP and AoLP in an annulus on the star
#TODO:    printandlog('\nMeasuring the star polarization in the background-subtracted I_Q-, I_U-, Q- and \nU-images:')
    q_star, u_star, DoLP_star, AoLP_star = determine_star_polarization(cube_I_Q=frame_I_Q_background_subtracted, cube_I_U=frame_I_U_background_subtracted, cube_Q=frame_Q_background_subtracted, cube_U=frame_U_background_subtracted, param_annulus_star=param_annulus_star, param_annulus_background=param_annulus_background)
       
    # Print resulting star polarization
    printandlog('\nMeasured star polarization in the background-subtracted I_Q-, I_U-, Q- and \nU-images:')
    printandlog('q_star = %.4f %%' % (100*q_star))
    printandlog('u_star = %.4f %%' % (100*u_star))
    printandlog('DoLP_star = %.4f %%' % (100*DoLP_star))
    printandlog('AoLP_star = %.2f deg' % AoLP_star)
    printandlog('\nThe measured star polarization should be very similar to that fitted above when \ncorrecting the instrumental polarization effects. The difference in q, u and \nDoLP should be <0.1% or <<0.1% depending on the noise in images.')

    # Subtract star polarization
#TODO:    printandlog('\nSubtracting the star polarization in the incident Q- and U-images.')
    frame_Q_star_polarization_subtracted = frame_Q_background_subtracted - q_star*frame_I_Q_background_subtracted
    frame_U_star_polarization_subtracted = frame_U_background_subtracted - u_star*frame_I_U_background_subtracted
    
    # Subtract very small residual background
    frame_Q_star_polarization_subtracted, background_frame_Q_star_polarization_subtracted = subtract_background(cube=frame_Q_star_polarization_subtracted, param_annulus_background=param_annulus_background)
    frame_U_star_polarization_subtracted, background_frame_U_star_polarization_subtracted = subtract_background(cube=frame_U_star_polarization_subtracted, param_annulus_background=param_annulus_background)
    
    # Print resulting background values
#TODO:    printandlog('\nSubtracting the residual background in the Q- and U-images after subtracting \nthe star polarization:')
    printandlog('\nSubtracted residual backgrounds in the star-polarization-subtracted Q- and U-\nimages:')
    printandlog('Background Q = ' + str(background_frame_Q_star_polarization_subtracted))
    printandlog('Background U = ' + str(background_frame_U_star_polarization_subtracted))

    ###############################################################################
    # Subtract background in I_Q-, I_U-, Q- and U-cubes
    ###############################################################################
    
    # Determine background in corrected I_Q-, I_U-, Q- and U-frames and subtract it
#TODO:    printandlog('\nSubtracting the backgrounds in the incident I_Q-, I_U-, Q- and U-image cubes:')
    cube_I_Q_background_subtracted, background_cube_I_Q = subtract_background(cube=cube_I_Q_incident, param_annulus_background=param_annulus_background)
    cube_I_U_background_subtracted, background_cube_I_U = subtract_background(cube=cube_I_U_incident, param_annulus_background=param_annulus_background)
    cube_Q_background_subtracted, background_cube_Q = subtract_background(cube=cube_Q_incident, param_annulus_background=param_annulus_background)
    cube_U_background_subtracted, background_cube_U = subtract_background(cube=cube_U_incident, param_annulus_background=param_annulus_background)
   
    # Print resulting background values
    printandlog('\nSubtracted mean backgrounds in the incident I_Q-, I_U-, Q- and U-image cubes:')
    printandlog('Mean background cube I_Q = ' + str(np.mean(background_cube_I_Q)))
    printandlog('Mean background cube I_U = ' + str(np.mean(background_cube_I_U)))
    printandlog('Mean background cube Q = ' + str(np.mean(background_cube_Q)))
    printandlog('Mean background cube U = ' + str(np.mean(background_cube_U)))
    
    ###############################################################################
    # Determine star polarization in cubes and plot it as function of HWP cycle number
    ###############################################################################    

    # Compute normalized Stokes q and u, DoLP and AoLP in an annulus on the star as a function of HWP cycle number
#TODO:    printandlog('\nDetermining the spread (standard deviation) of the measured star polarization \nwith HWP cycle number:')
    q_star_HWP_cycle, u_star_HWP_cycle, DoLP_star_HWP_cycle, AoLP_star_HWP_cycle = determine_star_polarization(cube_I_Q=cube_I_Q_background_subtracted, cube_I_U=cube_I_U_background_subtracted, cube_Q=cube_Q_background_subtracted, cube_U=cube_U_background_subtracted, param_annulus_star=param_annulus_star, param_annulus_background=param_annulus_background)
    
    # Compute spread of q_star, u_star, DoLP_star and AoLP_star
    sigma_q_star = np.std(q_star_HWP_cycle, ddof=1)
    sigma_u_star = np.std(u_star_HWP_cycle, ddof=1)
    sigma_DoLP_star = np.std(DoLP_star_HWP_cycle, ddof=1)
    sigma_AoLP_star = np.std(AoLP_star_HWP_cycle, ddof=1)
   
    # Print resulting spread of star polarization
    printandlog('\nMeasured spread (standard deviation) of the star polarization with HWP cycle \nnumber:')
    printandlog('sigma_q_star = %.4f %%' % (100*sigma_q_star))
    printandlog('sigma_u_star = %.4f %%' % (100*sigma_u_star))
    printandlog('sigma_DoLP_star = %.4f %%' % (100*sigma_DoLP_star))
    printandlog('sigma_AoLP_star = %.2f deg' % sigma_AoLP_star)
    
    # Plot q, u and DoLP from annulus as function of HWP cycle number
    plot_name_star_quDoLP = name_file_root + 'star_pol_quDoLP.png'
    printandlog('\nCreating plot ' + plot_name_star_quDoLP + ' showing the measured star \npolarization as a function of HWP cycle number.')          
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
    ax.set_ylabel(r'Meas. star polarization (%)', fontsize = font_size)
    ax.tick_params(axis = 'y', labelsize = font_size)  
    ax.grid()
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.savefig(os.path.join(path_reduced_dir, plot_name_star_quDoLP), dpi = 300, bbox_inches = 'tight')
    plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name_star_quDoLP), dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    # Plot AoLP from annulus as function of HWP cycle number
    plot_name_star_AoLP = name_file_root + 'star_pol_AoLP.png'
    printandlog('\nCreating plot ' + plot_name_star_AoLP + ' showing the measured \nangle of linear polarization of the star as a function of HWP cycle number.')
    plt.figure(figsize = (4.7, 3.0))
    plt.plot([0, x_max],[0, 0], '-k')     
    plt.plot(np.arange(1, x_max), AoLP_star_HWP_cycle, 'ok', label = 'AoLP') 
    plt.plot([0, x_max], [AoLP_star, AoLP_star], '-k')
    ax = plt.gca()
    ax.set_xlabel(r'HWP cycle', fontsize = font_size)
    ax.tick_params(axis = 'x', labelsize = font_size)
    ax.set_xlim([0, x_max])
    ax.set_ylabel(r'Meas. star AoLP ($^\circ$)', fontsize = font_size)
    ax.tick_params(axis = 'y', labelsize = font_size)  
    ax.grid()
    plt.legend(loc = 'best')
    plt.tight_layout()
    plt.savefig(os.path.join(path_reduced_dir, plot_name_star_AoLP), dpi = 300, bbox_inches = 'tight')
    plt.savefig(os.path.join(path_reduced_star_pol_subtr_dir, plot_name_star_AoLP), dpi = 300, bbox_inches = 'tight')
    plt.show()
    
    printandlog('\nHorizontal trends in the data points of the plots ' + plot_name_star_quDoLP + ' \nand ' + plot_name_star_AoLP + ' indicate that the instrumental \npolarization effects have been \nremoved successfully. However, this is only true \nprovided that:\n\n\
1) the observations are taken with a sufficiently large range of parallactic \nand altitude angles,\n\
2) the observations are taken with a sufficiently high signal-to-noise ratio, \nand \n\
3) the annulus for the star is placed in a region where there is only \nstarlight.\n\n\
A non-zero measured star polarization then indicates the star is truly \npolarized, which is often caused by the presence micron-sized particles in the \nline of sight. This star polarization can therefore indicate the presence of an \nunresolved (inner) circumstellar disk, starlight passing through a resolved \n(outer) part of a circumstellar disk or the presence of interstellar dust \nbetween the star and the Earth.')    

    ###############################################################################
    # Compute final images
    ###############################################################################
    
    printandlog('\n###############################################################################')
    printandlog('# Computing the final images and writing them to FITS-files')
    printandlog('###############################################################################') 

#TODO:    printandlog('\nComputing the final images both with the star polarization still present and \nwith the star polarization subtracted.')

    # Compute final images with the star polarization still present
    frame_I_tot, frame_Q_phi, frame_U_phi, frame_I_pol, frame_AoLP, frame_DoLP, frame_q, frame_u, frame_AoLP_norm, frame_DoLP_norm \
    = compute_final_images(frame_I_Q=frame_I_Q_background_subtracted, frame_I_U=frame_I_U_background_subtracted, frame_Q=frame_Q_background_subtracted, frame_U=frame_U_background_subtracted, header=header, images_north_up=images_north_up)

    # Compute final images with the star polarization subtracted   
    frame_Q_phi_star_polarization_subtracted, frame_U_phi_star_polarization_subtracted, frame_I_pol_star_polarization_subtracted, frame_AoLP_star_polarization_subtracted, frame_DoLP_star_polarization_subtracted, frame_q_star_polarization_subtracted, frame_u_star_polarization_subtracted, frame_AoLP_norm_star_polarization_subtracted, frame_DoLP_norm_star_polarization_subtracted \
    = compute_final_images(frame_I_Q=frame_I_Q_background_subtracted, frame_I_U=frame_I_U_background_subtracted, frame_Q=frame_Q_star_polarization_subtracted, frame_U=frame_U_star_polarization_subtracted, header=header, images_north_up=images_north_up)[1:]

    # Create frames that show annuli used to retrieve star and background signals   
    frame_annulus_star = compute_annulus_values(cube=frame_I_Q_background_subtracted, param=param_annulus_star)[1]
    frame_annulus_background = compute_annulus_values(cube=frame_I_Q_background_subtracted, param=param_annulus_background)[1]

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
    if tracking_mode_used == 'FIELD' and number_derotator_position_angles == 1 and images_north_up == False and np.sign(rotation_angle) != 0:
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
    
    if create_images_DoLP_AoLP_q_u_norm == True:
        # Add images of DoLP, normalized Stokes q and u and AoLP and DoLP created using q- and u-images
#TODO:        printandlog('\nAlso writing images of the degree of linear polarization, normalized Stokes q \nand u, and the degree and angle of linear polarization computed from q and u. \n\nWARNING\n\nNote that these images (DoLP.fits, q_norm.fits, u_norm.fits, AoLP_norm.fits and \nDoLP_norm.fits) only have meaning if all flux in the images originates from the \nastrophysical source of interest. This is generally the case for observations \nof for example solar system objects or galaxies. These images are generally not \nvalid for observations of circumstellar disks or companions because in that \ncase a large part of the flux in the total intensity images originates from the \ncentral star.')
        printandlog('\nWARNING, the images DoLP.fits, q_norm.fits, u_norm.fits, AoLP_norm.fits and \nDoLP_norm.fits are only valid if all flux in the images originates from the \nastrophysical source of interest. This is generally the case for observations \nof for example solar system objects or galaxies. The images are generally not \nvalid for observations of circumstellar disks or companions because in that \ncase a large part of the flux in the total intensity images originates from the \ncentral star.\n')
        frames_to_write += [frame_DoLP, frame_q, frame_u, frame_AoLP_norm, frame_DoLP_norm]       
        file_names += ['DoLP', 'q_norm', 'u_norm', 'AoLP_norm', 'DoLP_norm']
    
#TODO:    printandlog('\nWriting images to FITS-files.\n')

    # Write files of the images with the star polarization present
    for frame, file_name in zip(frames_to_write, file_names):
        write_fits_files(data=frame, path=os.path.join(path_reduced_dir, name_file_root + file_name + '.fits'), header=False)

    # Write frames that show annuli used to retrieve star and background signals in reduced directory
    write_fits_files(data=frame_annulus_star, path=os.path.join(path_reduced_dir, name_file_root + 'annulus_star.fits'), header=False)
    write_fits_files(data=frame_annulus_background, path=os.path.join(path_reduced_dir, name_file_root + 'annulus_background.fits'), header=False)    

    # List files of the images with the star polarization subtracted and define their file names
    frames_to_write = [frame_I_Q_background_subtracted, frame_I_U_background_subtracted, frame_I_tot, frame_Q_star_polarization_subtracted, frame_U_star_polarization_subtracted, frame_Q_phi_star_polarization_subtracted, frame_U_phi_star_polarization_subtracted, frame_I_pol_star_polarization_subtracted, frame_AoLP_star_polarization_subtracted]
    file_names = ['I_Q', 'I_U', 'I_tot', 'Q_star_pol_subtr', 'U_star_pol_subtr', 'Q_phi_star_pol_subtr', 'U_phi_star_pol_subtr', 'I_pol_star_pol_subtr', 'AoLP_star_pol_subtr']

    if create_images_DoLP_AoLP_q_u_norm == True:
        # Add images of DoLP, normalized Stokes q and u and AoLP and DoLP created using q- and u-images
        frames_to_write += [frame_DoLP_star_polarization_subtracted, frame_q_star_polarization_subtracted, frame_u_star_polarization_subtracted, frame_AoLP_norm_star_polarization_subtracted, frame_DoLP_norm_star_polarization_subtracted]       
        file_names += ['DoLP_star_pol_subtr', 'q_norm_star_pol_subtr', 'u_norm_star_pol_subtr', 'AoLP_norm_star_pol_subtr', 'DoLP_norm_star_pol_subtr']

    # Write files of the images with the star polarization subtracted
    for frame, file_name in zip(frames_to_write, file_names):
        write_fits_files(data=frame, path=os.path.join(path_reduced_star_pol_subtr_dir, name_file_root + file_name + '.fits'), header=False)

    # Write frames that show annuli used to retrieve star and background signals in reduced_star_pol_subtr directory
    write_fits_files(data=frame_annulus_star, path=os.path.join(path_reduced_star_pol_subtr_dir, name_file_root + 'annulus_star.fits'), header=False)
    write_fits_files(data=frame_annulus_background, path=os.path.join(path_reduced_star_pol_subtr_dir, name_file_root + 'annulus_background.fits'), header=False)    

###############################################################################
###############################################################################
## Run main functions
###############################################################################
###############################################################################

# Turn off pyfits FITS-warnings
warnings.filterwarnings('ignore', category=UserWarning)

# Define path of log file and delete it if it already exists
path_log_file = os.path.join(path_main_dir, name_file_root + 'LOG.txt')

if os.path.exists(path_log_file) == False:
    # Start log file
    printandlog('\n###############################################################################')
    printandlog('# Starting pre-processing')
    printandlog('###############################################################################') 
elif os.path.exists(path_log_file) == True and skip_preprocessing == False:
    # Empty existing log file and start writing
    open(path_log_file, 'w').close()
    printandlog('\n###############################################################################')
    printandlog('# Starting pre-processing')
    printandlog('###############################################################################') 
    printandlog('\nWARNING, the previous log file has been deleted.')
elif os.path.exists(path_log_file) == True and skip_preprocessing == True:
    # Remove all lines concerned with the post-processing from the log file
    log_file_lines = [x.rstrip('\n') for x in open(path_log_file, 'r')]
    log_file_lines = log_file_lines[:log_file_lines.index('# Starting post-processing') - 1]
    open(path_log_file, 'w').close()
    for line in log_file_lines:                            
        print(line, file=open(path_log_file, 'a'))
   
if skip_preprocessing == False:
    # Pre-process raw data
    cube_single_sum, cube_single_difference, header = perform_preprocessing(param_annulus_background_flux=param_annulus_background_flux, sigmafiltering_sky=sigmafiltering_sky, sigmafiltering_object=sigmafiltering_object, sigmafiltering_flux=sigmafiltering_flux, save_preprocessed_data=save_preprocessed_data)
    
    # Print that post-processing starts
    printandlog('\n###############################################################################')
    printandlog('# Starting post-processing')
    printandlog('###############################################################################') 
                
elif skip_preprocessing == True:
    # Define paths to read pre-processed data and headers from
    path_cube_single_sum = os.path.join(path_preprocessed_dir, 'cube_single_sum.fits')
    path_cube_single_difference = os.path.join(path_preprocessed_dir, 'cube_single_difference.fits')
    path_object_files_text = os.path.join(path_preprocessed_dir, 'path_object_files.txt')
   
    if os.path.exists(path_cube_single_sum) and os.path.exists(path_cube_single_difference) and os.path.exists(path_object_files_text):
        # Print that post-processing starts
        printandlog('\n###############################################################################')
        printandlog('# Starting post-processing')
        printandlog('###############################################################################') 
        printandlog('\nSkipping pre-processing and reading pre-processed data and headers.\n')
        
        # Read pre-processed single-sum and difference- images        
        cube_single_sum = read_fits_files(path=path_cube_single_sum, silent=False)[0]
        cube_single_difference = read_fits_files(path=path_cube_single_difference, silent=False)[0] 

        # Read headers
        header = [pyfits.getheader(x.rstrip('\n')) for x in open(path_object_files_text, 'r')]
        printandlog('Read headers from OBJECT-files specified in ' + path_object_files_text + '.')
    
    else:
        raise IOError('The files ' + path_cube_single_sum + ', ' + path_cube_single_difference + ' and/or ' + path_object_files_text + ' do not exist. Set skip_preprocessing to False and save_preprocessed_data to True to perform the pre-processing of the raw data and save the results.')

# Perform post-processing of data           
perform_postprocessing(cube_single_sum=cube_single_sum, cube_single_difference=cube_single_difference, header=header, param_annulus_star=param_annulus_star, param_annulus_background=param_annulus_background, double_difference_type=double_difference_type, remove_vertical_band_detector_artefact=remove_vertical_band_detector_artefact, combination_method_polarization_images=combination_method_polarization_images, trimmed_mean_proportiontocut_polarization_images=trimmed_mean_proportiontocut_polarization_images, combination_method_total_intensity_images=combination_method_total_intensity_images, trimmed_mean_proportiontocut_total_intensity_images=trimmed_mean_proportiontocut_total_intensity_images, images_north_up=images_north_up, create_images_DoLP_AoLP_q_u_norm=create_images_DoLP_AoLP_q_u_norm)

# Print time elapsed
time_end = time.time()
d = datetime.datetime(1, 1, 1) + datetime.timedelta(seconds = time_end - time_start)
printandlog('\nTime elapsed: %d h %d min %d s' % (d.hour, d.minute, d.second)) 

               















###############################################################################
###############################################################################
## Post-process CT Cha data
###############################################################################
###############################################################################


##TODO:TODO:TODO:TODO: CT Cha
## Paths to read from
#file_science_reduced_single_sum = [x for x in next(os.walk(path_raw_dir + '\\Single Sum\\'))[2] if x.endswith('.fits')]
#path_science_reduced_single_sum = [path_raw_dir + '\\Single Sum\\' + x for x in file_science_reduced_single_sum]
#file_science_reduced_single_difference = [x for x in next(os.walk(path_raw_dir + '\\Single Difference\\'))[2] if x.endswith('.fits')]
#path_science_reduced_single_difference = [path_raw_dir + '\\Single Difference\\' + x for x in file_science_reduced_single_difference]
#
##TODO: Remove first 2 HWP cycles because they are only Qplus/minus. Maybe later add them again.
#path_science_reduced_single_sum = path_science_reduced_single_sum[4:]
#path_science_reduced_single_difference = path_science_reduced_single_difference[4:]
#
################################################################################
## ReadFits
################################################################################
#
#def ReadFits(path):
#    ''' 
#    Read .fits-files
#    
#    Input:
#        path: string or list of strings specifying paths to read .fits-files from 
#    
#    Output:
#        data: image data cube (when 1 path specified) or list of image data 
#              cubes (when more than 1 path specified) read from .fits-files
#        header: headers (when 1 path specified) or list of headers 
#                (when more than 1 path specified) read from .fits-files
#    
#    Note:
#        If the data read is a single frame, an additional dimension is added to
#        make it a 3D image data cube (third dimension of length 1).
#                
#    File written by Rob van Holstein
#    Function status: verified
#    '''
#    
#    if type(path) == list and len(path) == 1:
#        path = path[0]
#    
#    if type(path) == str:
#        hdulist = pyfits.open(path)
#        data = hdulist[0].data
#        if data.ndim == 2:
#            data = np.expand_dims(data, axis = 0)
#        header = hdulist[0].header
#        hdulist.close()
#        printandlog('Read file', path)
#            
#    elif type(path) == list:
#        data = []
#        header = []
#        for path_sel in path:
#            hdulist = pyfits.open(path_sel)
#            data0 = hdulist[0].data
#            if data0.ndim == 2:
#                data0 = np.expand_dims(data0, axis = 0)
#            data.append(data0)
#            header.append(hdulist[0].header) 
#            hdulist.close()
#            printandlog('Read file ' + path_sel)
#            
#    return data, header
#
## Read single sum and difference frames and their headers and stack frames  
#printandlog('')
#data_single_sum, header_array = ReadFits(path=path_science_reduced_single_sum)
#cube_single_sum = np.vstack(data_single_sum)
#data_single_difference, _ = ReadFits(path=path_science_reduced_single_difference)
#cube_single_difference = np.vstack(data_single_difference)
#
## Filter for zeros and NaN's
#single_sum = np.nan_to_num(cube_single_sum)
#single_sum[cube_single_sum == 0] = 1
#single_diff = np.nan_to_num(cube_single_difference)
#single_diff[cube_single_difference == 0] = 1
#
####TODO: Make CT Cha data set have just one derotator offset
##single_sum = single_sum[:24, :, :]
##single_diff = single_diff[:24, :, :]
##header_array = header_array[:24]
#
#perform_postprocessing(cube_single_sum=single_sum, cube_single_difference=single_diff, header=header_array, param_annulus_star=param_annulus_star, param_annulus_background=param_annulus_background, double_difference_type=double_difference_type, remove_vertical_band_detector_artefact=remove_vertical_band_detector_artefact, combination_method_polarization_images=combination_method_polarization_images, trimmed_mean_proportiontocut_polarization_images=trimmed_mean_proportiontocut_polarization_images, combination_method_total_intensity_images=combination_method_total_intensity_images, trimmed_mean_proportiontocut_total_intensity_images=trimmed_mean_proportiontocut_total_intensity_images, images_north_up=images_north_up, create_images_DoLP_AoLP_q_u_norm=create_images_DoLP_AoLP_q_u_norm)
#
#time_end = time.time()
#d = datetime.datetime(1, 1, 1) + datetime.timedelta(seconds = time_end - time_start)
#printandlog('\nTime elapsed: %d h %d min %d s' % (d.hour, d.minute, d.second)) 






################################################################################
## Create folder and copy script
################################################################################
#
#def create_folder_copy_script(path_write_folder):
#    '''
#    Create a folder to write files to and copy used script to this folder
#    
#    Input:
#        path_write_folder: path of folder to write .png- and .pdf images to; if None do not write files
#        
#    Output:
#        path_write_folder_new: new path of write folder if folder already existed
#            
#    Written by Rob van Holstein
#    Function status: verified   
#    '''
#
##TODO: Check if path_current_script and name_current_script still work if this function is put in separate file
#    
#    # Create folder to write files to
#    if not os.path.exists(path_write_folder):
#        # If the folder does not exist yet, create it
#        os.makedirs(path_write_folder)
#        path_write_folder_new = path_write_folder
#    else:
#        # If the folder exists, check if folders with added numbers exist and create a folder with a higher number
#        name_folder_addition = 2
#        path_write_folder_new = path_write_folder + ' ' + str(name_folder_addition)
#        
#        while os.path.exists(path_write_folder_new):
#            name_folder_addition += 1
#            path_write_folder_new = path_write_folder + ' ' + str(name_folder_addition)
#            
#        os.makedirs(path_write_folder_new)
#    
#    print('\nWriting files in', path_write_folder_new)
#    
#    # Copy script that is run to output folder    
#    path_current_script = os.path.abspath(__file__)
#    name_current_script = os.path.basename(__file__)
#    path_copy_script = path_write_folder_new + '\\' + name_current_script[:name_current_script.rfind('.')] + \
#                       '_as_run' + name_current_script[name_current_script.rfind('.'):]
#    
#    shutil.copy2(path_current_script, path_copy_script)
#    
#    return path_write_folder_new
#
