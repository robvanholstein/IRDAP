'''

@author: Rob the Modellist and Christian the Artist
'''

###############################################################################
###############################################################################
## Input
###############################################################################
###############################################################################

# Definition of paths
path_main_dir = r'C:\Users\Rob\Desktop\IP Measurement Test\GQ Lup'
#path_main_dir = r'C:\Users\Rob\Desktop\IP Measurement Test\GG Tau'
path_static_flat_badpixelmap = r'C:\Users\Rob\Documents\PhD\CentralFiles\irdis_polarimetry_pipeline'

# Options for pre-processing
skip_preprocessing = False
sigmafiltering = False
centering_method_object = 'center frames' # 'center frames', 'gaussian', 'cross-correlation', 'manual'
centering_subtract_object = True
center_coordinates_object = (478, 522, 1504, 512) # 'default'
#param_centering_object = (12, 7, None) # 'default' Coords need to be accurate within a few pixels; for manual without dithering
param_centering_object = (12, None, None) # 'default' Coords need to be accurate within a few pixels; for manual without dithering
collapse_ndit_object = True
plot_centering_sub_images = True
center_coordinates_flux = () #TODO: or list of center coordinates if different for multiple frames

# Flux
centering_method_flux = 'gaussian' # 'gaussian', 'manual'
center_coordinates_flux = (445, 491, 1470, 480)
param_centering_flux = (12, None, 30000)
param_annulus_background_flux = 'large annulus'

save_preprocessed_data = True

frames_to_remove = [(1, 3),
                    (2, 2),
                    (15, 1),     
                    4, 
                    7,
                    16,  
                    (17, 2),
                    (18, 1),
                    (19, 2),
                    (19, 6),
                    (20, 7)] # GQ Lup

# Options for post-processing
double_difference_type = 'standard'
remove_vertical_band_detector_artefact = True
param_annulus_star = 'ao residuals'
param_annulus_star = [(512.5, 512.5, 50, 20, -65, 130),  # 'ao residuals', 'star aperture'
                      (512.5, 512.5, 50, 20, -175, -95)] # GQ Lup
#param_annulus_star = (516, 478, 0, 11, 0, 360)
#param_annulus_background = 'large annulus'
combination_method_polarization_images = 'trimmed mean'
#combination_method_polarization_images = 'least squares'
trimmed_mean_proportiontocut_polarization_images = 0.10
combination_method_total_intensity_images = 'mean'
trimmed_mean_proportiontocut_total_intensity_images = 0.10
images_north_up = True
create_images_DoLP_AoLP_q_u_norm = True

#TODO: turn combination_method_polarization_images and trimmed_mean_proportiontocut_polarization_images
# into a single variable (also for the total intensity ones)? Instead of trimmed mean, give a number, 
# also for mean can probably give 0 (should test it though).

#frames_to_remove = [(1, 1), 
#                    (1, 2), 
#                    (1, 4), 
#                    (1, 10),
#                    (2, 6), 
#                    (2, 7),
#                    3, 
#                    (4, 6),
#                    5,
#                    (7, 3),
#                    (7, 6),
#                    (7, 15), 
#                    (8, 3)]

#frames_to_remove = [(1, 1),
#                    (1, 4),
#                    (8, 3), 
#                    (9, 1),
#                    (10, 2),
#                    (11, 3), 
#                    (12, 4),
#                    (12, 5)]

#frames_to_remove = []

###############################################################################
###############################################################################
## README
###############################################################################
###############################################################################

'''
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
- x-coordinate of center (pixels, 1-based, center of image is at 512.5)
- y-coordinate of center (pixels, 1-based, center of image is at 512.5)
- inner radius (pixels)
- width (pixels)
- start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
- end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
For example, when using 'star aperture' the annulus used will be: 
(512.5, 512.5, 0, 11, 0, 360). The coordinates and angles of the annulus are
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
- x-coordinate of center (pixels, 1-based, center of image is at 512.5)
- y-coordinate of center (pixels, 1-based, center of image is at 512.5)
- inner radius (pixels)
- width (pixels)
- start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
- end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
For example, when using 'large annulus' the annulus used will be: 
(512.5, 512.5, 360, 60, 0, 360). The coordinates and angles of the annulus are
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
###############################################################################
## Function definitions
###############################################################################
###############################################################################

# Import packages
import os
import glob
import time
import datetime
import warnings
import textwrap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from scipy.ndimage.interpolation import rotate
from scipy import interpolate
from scipy import optimize
from scipy import ndimage
from scipy.stats import trim_mean
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
from skimage.transform import rotate as rotateskimage
from skimage.feature import register_translation

###############################################################################
# printandlog
###############################################################################

def printandlog(single_object, wrap=True):
    '''
    Print a single object (string) on screen and save it to a log file. The log
    file is located at path_main_dir and the name starts with name_file_root, 
    which are strings and global variables to the funtion.
    
    Input:
        single_object: single object (string) to be printed and logged
        wrap: if True, wrap print statements to a maximum of 80 character. If 
            False, do not wrap print statements (default = True).
           
    File written by Rob van Holstein
    Function status: verified
    '''
            
    # Define path to write log file, using path_reduced_dir as global variable
    path_log_file = os.path.join(path_main_dir, name_file_root + 'LOG.txt')

    if not os.path.exists(path_log_file):
        # Create log file
        open(path_log_file, 'w+')
    
    # Wrap string to not exceed 79 characters
    if type(single_object) == str and wrap == True:
        single_object = textwrap.fill(single_object, width=80, 
                                      replace_whitespace=False, 
                                      break_long_words=False)
    
    # Print object in log file and on screen
    print(single_object, file=open(path_log_file, 'a'))
    print(single_object) 

###############################################################################
# create_overview_headers
###############################################################################

def create_overview_headers():
    '''
    Create an overview of relevant FITS-headers and write it to a text-file
    
    Note that path_raw_diris a global variable to the function.
        
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
    path_overview = os.path.join(path_main_dir, name_file_root + 'headers.txt')
    np.savetxt(path_overview, print_array, fmt = '%s', newline= '\r\n')
    printandlog('\nWrote file ' + path_overview + ' showing an overview of relevant headers for each file in the raw directory.')
    
###############################################################################
# check_sort_data_create_directories
###############################################################################

def check_sort_data_create_directories(frames_to_remove, 
                                       combination_method_polarization_images='trimmed mean', 
                                       centering_method_object='center frames', 
                                       save_preprocessed_data=True, 
                                       plot_centering_sub_images=True):
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
        combination_method_polarization_images: method to be used to produce the 
            incident Q- and U-images, 'least squares', 'trimmed mean' or 'median' 
            (default = 'trimmed mean')
        centering_method_object: method to center the OBJECT-frames. If 
            'center frames' or 'manual', use fixed coordinates as provided by 
            center_coordinates. If 'gaussian', fit a 2D Gaussian to each frame. 
            If 'cross-correlation', fit a 2D Gaussian to the first frame and then
            use cross-correlation to align (register) the other frames onto the 
            centered first  frame. For 'gaussian' and 'cross-correlation' 
            center_coordinates is used as initial guess of the center 
            coordinates and the determined center coordinates are plotted for
            each image (default = 'center frames').
        save_preprocessed_data: If True, save preprocessed cubes of single-sum 
            and single-difference images in the 'preprocessed' folder so that 
            the preprocessing can be skipped when re-running the pipeline
            (default = True).
        plot_centering_sub_images: If True, plot the sub-images showing the 
            center coordinates for each frame. The plots allow for checking
            whether the centering is correct and to scan the data for frames 
            with bad quality (default = True). 
    
    Note that combination_method_polarization_images, centering_method_object, 
    save_preprocessed_data and plot_centering_sub_images are input to this 
    function as they are required for the sorting of the data, not to actually
    perform centering for example.
    
    Note that path_raw_dir, path_sky_dir, path_center_dir, path_flux_dir, 
    path_sky_flux_dir, path_preprocessed_dir, path_reduced_dir and 
    path_reduced_star_pol_subtr_dir are global variables to the function.
    
    Output:
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
        combination_method_polarization_images: method to be used to produce the 
            incident Q- and U-images, 'least squares', 'trimmed mean' or 'median'
     
    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: verified
    '''

    ###############################################################################
    # Obtain raw files and headers
    ###############################################################################
   
    # Extract paths to FITS-files in raw directory
    path_raw_files = glob.glob(os.path.join(path_raw_dir,'*.fits'))
    
    # Check if raw folder contains FITS-files
    if len(path_raw_files) == 0:
        raise IOError('The raw directory {0:s} does not contain FITS-files. You need to put your raw FITS-files in this folder.'.format(path_raw_dir))
    
    # Create overview of relevant FITS-headers in a text file
    printandlog('\nReading raw directory ' + path_raw_dir + '.')
    create_overview_headers()

    # Extract headers
    header = [pyfits.getheader(x) for x in path_raw_files]
    
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
        raise ValueError('One or more files to be removed are listed more than once.')
        
    # Raise error if files to remove do not exist
    if not all([x in file_index for x in files_to_remove]):
        raise ValueError('One or more files to be removed do not exist.')
    
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

    if len(set([x['ESO OBS TARG NAME'] for x in header])) != 1:
        raise IOError('The data provided have different targets.')
    
    if not all([x['ESO DPR TYPE'] in ['OBJECT', 'SKY', 'OBJECT,CENTER', 'OBJECT,FLUX'] for x in header]):
        raise IOError('One or more files are not of type OBJECT, SKY, OBJECT,CENTER or OBJECT,FLUX.')
        
    if len(set([x['ESO INS4 COMB ROT'] for x in header])) != 1:
        raise IOError('The data provided use different tracking modes.')

    if not header[0]['ESO INS4 COMB ROT'] in ['FIELD', 'PUPIL']:
        raise IOError('The tracking mode used is not field-tracking or pupil-tracking.')
        
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
        printandlog('\nWARNING, the OBJECT-files use a NIR neutral density filter. The final data product will be less accurate because the neutral density filters are known to have a depolarizing effect that is not calibrated.')
    
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
        printandlog('\nWARNING, there are no SKY-files to subtract from the OBJECT-files. Although the background will be subtracted from the final images after determining it using the annulus as defined by the input variable \'param_annulus_background\', the result will be less accurate than when subtracting a SKY-image.')    
    
    if any(header_flux):    
        if not any([x['ESO DPR TYPE'] == 'SKY' and x['EXPTIME'] == flux_exposure_time and x['ESO INS4 FILT2 NAME'] == flux_nd_filter for x in header]):
            printandlog('\nWARNING, there are no SKY-files to subtract from the FLUX-file(s). Although the background will be subtracted after determining it using the annulus as defined by the input variable \'param_annulus_background_flux\', the result will be less accurate than when subtracting a SKY-image.')
    
    # Perform checks on header values that apply only to CENTER files
    header_center = [x for x in header if x['ESO DPR TYPE'] == 'OBJECT,CENTER']
    
    if any(header_center):
        if not all([x['EXPTIME'] == object_exposure_time for x in header_center]):
            raise IOError('One or more CENTER-files have an exposure time different from that of the OBJECT-files.')
        
        if not all([x['ESO INS4 FILT2 NAME'] == object_nd_filter for x in header_center]):
            raise IOError('One or more CENTER-files use a NIR neutral density filter different from that of the OBJECT-files.')
        
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
        raise ValueError('One or more frames to be removed are listed more than once.')

    # Raise error if frames need to be removed from a file that is already completely removed
    if any([x in files_to_remove for x in files_to_remove_frames_from]):
        raise ValueError('One or more files to remove frames from are already completely removed.')  

    # Raise error if files to remove do not exist
    if not all([x in file_index for x in files_to_remove_frames_from]):
        raise ValueError('One or more files to remove frames from do not exist.')

    # Raise error if for one or more files the indices of frames to remove are negative or zero
    if any([x < 0 for x in frames_to_remove_per_file]):
        raise ValueError('One or more of the frame numbers to be removed are negative or zero.')

    # Create list of arrays with indices of frames to remove for each file
    indices_to_remove = [np.sort(frames_to_remove_per_file[files_to_remove_frames_from == x]).astype(np.int) for x in file_index]

    # Extract NDIT of each file
    NDIT = [x['ESO DET NDIT'] for x in header]
    
    # Raise error if one or more files have all the frames removed
    if any([np.array_equal(x, np.arange(0, y)) for x,y in zip(indices_to_remove, NDIT)]):
        raise ValueError('One or more files have all frames removed. If you want to remove the complete file, please specify just the file number in frames_to_remove and not a tuple of frames.')
    
    # Raise error if for one or more files the indices of frames to remove are outside of the NDIT
    if any([any(x >= y) for x,y in zip(indices_to_remove, NDIT)]):
        raise ValueError('One or more of the frame numbers to be removed are higher than the NDIT of the corresponding file.')

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
    path_object_files = []
    path_sky_files = []
    path_center_files = []
    path_flux_files = []
    path_sky_flux_files = [] 
    path_imcompatible_files = []
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
    
        if header_sel['ESO DPR TYPE'] == 'OBJECT':
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
            separation = 9 - len(str(file_index_object[file_to_remove]))
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
            raise IOError('The data has no U-measurements and therefore cannot be reduced.')
        if any([x in stokes_parameter for x in ['Uplus', 'Uminus']]):
            raise IOError('The data has no Q-measurements and therefore cannot be reduced.')
        
    # Force 'least squares' for combining the polarization images if the number of Q- and U-measurements is unequal
    if combination_method_polarization_images != 'least squares':
        if stokes_parameter.count('Qplus') != stokes_parameter.count('Uplus'):
            combination_method_polarization_images = 'least squares'
            printandlog('\nWARNING, combination_method_polarization_images is forced to \'least squares\' because the number of Q- and U-measurements are unequal. The plots showing the measured star polarization as a function of HWP cycle number will only show the complete HWP cycles.')
 
    ###############################################################################
    # Print number of files for each file type
    ###############################################################################

    # Print number of files for each file type
    printandlog('\nNumber of files for each file type:')
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
    
    # Raise error when there are no center files, but they are required by the selected centering method
    if centering_method_object == 'center frames' and not any(path_center_files):
        raise IOError('centering_method_object = \'{0:s}\', but there are no CENTER-files provided.'.format(centering_method_object))

    ###############################################################################
    # Create directories to write processed data to
    ###############################################################################

    # Create directories to write processed data to
    directories_created = []
    directories_already_existing = []
    
    if any(path_sky_files):
        if not os.path.exists(path_sky_dir):
            os.makedirs(path_sky_dir)
            directories_created.append(path_sky_dir)
        else:
            directories_already_existing.append(path_sky_dir)
    
    if any(path_center_files) and centering_method_object == 'center frames':
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
       
    if save_preprocessed_data == True or plot_centering_sub_images == True or \
       centering_method_object in ['gaussian', 'cross-correlation']:
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
    
    return path_object_files, path_sky_files, path_center_files, path_object_center_files, \
           path_flux_files, path_sky_flux_files, indices_to_remove_object, indices_to_remove_sky, \
           indices_to_remove_center, indices_to_remove_object_center, indices_to_remove_flux, \
           indices_to_remove_sky_flux, file_index_object, file_index_flux, \
           combination_method_polarization_images

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

def remove_bad_pixels(cube, frame_master_bpm, sigmafiltering=True):
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
    
    # If the input is a frame turn the resulting cube back into a frame
    if cube_ndim == 2:
        cube_filtered = np.squeeze(cube_filtered)
       
    return cube_filtered

###############################################################################
# process_sky_frames
###############################################################################

def process_sky_frames(path_sky_files, indices_to_remove_sky, frame_master_bpm, sigmafiltering=True):
    '''
    Create a master sky-frame from the SKY-files
    
    Input:
        path_sky_files: string or list of strings specifying paths of SKY-files
        indices_to_remove_sky: list of arrays with indices of frames to remove for each SKY-file
        frame_master_bpm: frame indicating location of bad pixels with 0's and good
            pixels with 1's
        sigmafiltering: if True remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True)

    Output:
        frame_master_sky: master sky frame
    
    File written by Rob van Holstein; based on a function by Christian Ginski
    Function status: verified       
    '''
      
    # Read sky frames
    list_cube_sky_raw = []
        
    for path_sel, indices_sel in zip(path_sky_files, indices_to_remove_sky):
        # Read data and header from file
        cube_sel = read_fits_files(path=path_sel, silent=True)[0]

        # Remove frames based on list of indices provided by the user
        cube_sel = np.delete(cube_sel, indices_sel, axis=0)

        # Append cube to list
        list_cube_sky_raw.append(cube_sel)
    
    # Make a single image cube from list of image cubes or frames
    cube_sky_raw = np.vstack(list_cube_sky_raw)

    # Remove bad pixels of each frame
    cube_sky_filtered = remove_bad_pixels(cube=cube_sky_raw, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering)

    # Compute median of sky frames
    frame_master_sky = np.median(cube_sky_filtered, axis=0)
            
    printandlog('\nThe master sky frame was created out of ' + str(len(path_sky_files)) + 
                ' raw SKY-file(s) comprising a total of ' + str(cube_sky_raw.shape[0]) + ' frame(s).')
    
    return frame_master_sky

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
    
    # Define pixel scale of IRDIS detector ("/px)
    pixel_scale = 12.25e-3

    # Define approximate separation of satellite spots (pixels) and filter central wavelength and bandwidth (m)
    if filter_used == 'FILT_BBF_Y':
        separation = 31.0
        lambda_c = 1043e-9
        delta_lambda = 140e-9
    elif filter_used == 'FILT_BBF_J':
        separation = 37.4
        lambda_c = 1245e-9
        delta_lambda = 240e-9
    elif filter_used == 'FILT_BBF_H':
        separation = 48.5
        lambda_c = 1625e-9
        delta_lambda = 290e-9
    elif filter_used == 'FILT_BBF_Ks':
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
                          centering_subtract_object=True, 
                          center_coordinates=(477, 521, 1503, 511), 
                          sigmafiltering=True):
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
        centering_subtract_object: if True subtract the OBJECT-file taken
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
        sigmafiltering: if True remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True)
       
    Output:
        list_frame_center_processed: list of processed center frames
        header: list of headers of center files
                
    File written by Rob van Holstein
    Function status: verified
    '''

    # Print if and how background is subtracted
    if centering_subtract_object == True:
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
   
        if centering_subtract_object == True:
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
                printandlog('\nRotating OBJECT-frame around the initial guess of center before subtracting it from the CENTER-file.')
                printandlog('')
                frame_background = remove_bad_pixels(cube=frame_background, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering)
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
        cube_badpixel_filtered = remove_bad_pixels(cube=cube_bgsubtr_flatfielded, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering)
                
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
        raise ValueError('frame should be a 2-dimensional array.')
    if type(crop_radius) != int and crop_radius is not None:
        raise TypeError('crop_radius must be integer or None.') 
    if type(x0) not in [int, np.int32, np.int64] and y0 is not None:
        raise TypeError('x0 must be integer or None.')
    if type(y0) not in [int, np.int32, np.int64] and x0 is not None:
        raise TypeError('y0 must be integer or None.')

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
                            param_centering=(12, 7, 30000)):
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
        param_centering: length-3-tuple with parameters for centering:
            crop_radius: half the length of side of square cropped sub-images used 
                to fit Gaussian to (pixels). Must be integer. If None, the complete 
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
            The default value of param_centering is (12, 7, 30000).
                
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
    crop_radius = param_centering[0]
    sigfactor = param_centering[1]
    saturation_level = param_centering[2]
    
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
        title_main = 'center_coordinates = %s' % (center_coordinates,)
        if center_coordinates == (477, 521, 1503, 511):
            title_main += ' (default)'
        title_main += '\nparam_centering = %s' % (param_centering,)         
        if param_centering == (12, 7, 30000):
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
        plt.show()
               
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
        printandlog('\nWrote file ' + path_reg_file + ' which can be loaded as a region in DS9 and shows the fitted coordinates of the satellite spots and the center.')

    # Print center coordinates found
    max_path_length = max([len(os.path.basename(x)) for x in path_processed_center_files])
    separator_top = ' '*(max_path_length - 9)
    printandlog('\nMeasured center coordinates per file (pixels):')
    printandlog('File name' + separator_top + '    x_left    y_left    x_right    y_right', wrap=False)
    for i, path_sel in enumerate(path_processed_center_files):
        separator = ' '*(max_path_length - len(os.path.basename(path_sel)))
        printandlog(os.path.basename(path_sel) + separator + '    %.2f    %.2f    %.2f    %.2f' 
                    % (x_center_fit[i, 0], y_center_fit[i, 0], x_center_fit[i, 1], y_center_fit[i, 1]), wrap=False)    
    
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
                    % (x_center[0], x_center_std[0], y_center[0], y_center_std[0], x_center[1], x_center_std[1], y_center[1], y_center_std[1]), wrap=False)    
        
        # Print warning if there is significant deviation among the center coordinates found
        if any(np.append(x_center_std, y_center_std) > 0.5):
            printandlog('\nWARNING, for at least one of the fitted center coordinates the standard deviation is larger than 0.5 pixels.')
    else:
        # Print mean center coordinates without error
        printandlog('\nFinal center coordinates (pixels):')
        printandlog('x_left    y_left    x_right    y_right')
        printandlog('%.2f    %.2f    %.2f    %.2f' % (x_center[0], y_center[0], x_center[1], y_center[1]))    

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
                          sigmafiltering=True, 
                          centering_method='center frames', 
                          center_coordinates=(477, 521, 1503, 511), 
                          param_centering=(12, 7, 30000), 
                          collapse_ndit=True, 
                          plot_centering_sub_images=True):
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
        sigmafiltering: if True remove bad pixels remaining after applying
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
        param_centering: length-3-tuple with parameters for centering:
            crop_radius: half the length of the sides of the square cropped 
                sub-images used to fit the 2D Gaussian to and used for 
                cross-correlating the images (pixels). Must be integer. The 
                sub-image is centered on the coordinates as provided by
                center_coordinates. If None, the complete frame is used for the 
                fitting and center_coordinates is ignored. The value of
                crop_radius is also used to create the sub-images when
                plot_centering_sub_images = True.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of param_centering is (12, 7, 30000).
        collapse_ndit: If True, compute the mean over the (NDIT) frames of a
            file before subtracting the background, flat-fielding, bad pixel
            removal and centering. If False, perform the above steps for each
            frame and after that compute the mean over the frames 
            (default = True).
        plot_centering_sub_images: If True, plot the sub-images showing the 
            center coordinates for each frame. The plots allow for checking
            whether the centering is correct and to scan the data for frames 
            with bad quality (default = True). 
        
    Output:
        cube_single_sum: cube of single-sum I_Q^+, I_Q^-, I_U^+ and I_U^- intensity images
        cube_single_difference: cube of single-difference Q^+, Q^-, U^+ and U^- images
        header: list of FITS-headers of raw science frames   
                
    File written by Rob van Holstein; based on function by Christian Ginski
    Function status: verified
    '''
    
    # Print centering method selected
    if centering_method == 'center frames':
        printandlog('\nCentering frames with center coordinates found from CENTER-file(s):')
        printandlog('(%.2f, %.2f, %.2f, %.2f)' % center_coordinates)
    elif centering_method == 'gaussian':
        printandlog('\nCentering frames by fitting a 2D Gaussian.')
    elif centering_method == 'cross-correlation':
        printandlog('\nCentering frames using cross-correlation.')
    elif centering_method == 'manual':
        printandlog('\nCentering frames with user-provided center coordinates:')
        printandlog('(%.2f, %.2f, %.2f, %.2f)' % center_coordinates)

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
    list_single_sum = []
    list_single_difference = []
    header = []

    if centering_method in ['gaussian', 'cross-correlation']:
        # Determine filter used and compute theoretical full width half maximum
        filter_used = pyfits.getheader(path_object_files[0])['ESO INS1 FILT ID']
        fwhm = compute_fwhm_separation(filter_used)[0]
    
    if plot_centering_sub_images == True:
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
        cube_badpixel_filtered = remove_bad_pixels(cube=cube_bgsubtr_flatfielded, frame_master_bpm=frame_master_bpm, sigmafiltering=sigmafiltering)

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

                    if plot_centering_sub_images == True:
                        # Create sub-image to show center coordinates in
                        x_center_0_rounded = int(np.round(x_center_0_sel))
                        y_center_0_rounded = int(np.round(y_center_0_sel))
                        sub_image = create_sub_image(frame=frame_half, x0=x_center_0_rounded + x_dith, y0=y_center_0_rounded + y_dith, crop_radius=crop_radius)
    
                        # Compute center position in coordinates of sub-image 
                        x_fit_sub_image = x_center_0_sel - x_center_0_rounded + crop_radius
                        y_fit_sub_image = y_center_0_sel - y_center_0_rounded + crop_radius

                elif centering_method == 'gaussian':
                    # Determine accurate coordinates of star position fitting a Gaussian           
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

                    if plot_centering_sub_images == True:
                        # Compute fit position in coordinates of sub-image                   
                        x_fit_sub_image = x_fit_sub_image_template[k] - x_shift_fit
                        y_fit_sub_image = y_fit_sub_image_template[k] - y_shift_fit
                        
                if plot_centering_sub_images == True:
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

        # Create single difference and sum image
        frame_single_sum = frame_left_centered + frame_right_centered
        frame_single_difference = frame_left_centered - frame_right_centered
        
        # Append single sum and difference images and header
        list_single_sum.append(frame_single_sum)
        list_single_difference.append(frame_single_difference)
        header.append(header_sel)
        
        # Print which file has been processed
        printandlog('Processed file ' + str(i + 1) + '/' + str(len(path_object_files)) + ': {0:s}'.format(os.path.basename(path_sel)))

    # Convert lists of single sum and difference images to image cubes
    cube_single_sum = np.stack(list_single_sum)
    cube_single_difference = np.stack(list_single_difference)    
        
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
                printandlog('\nCreating plot ' + path_plot + ' showing the center coordinates in the ' 
                            + x_y + '-direction of the ' + left_right + ' frame halves.')
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
                plt.show()

            # Plot center coordinates in x- and y-direction of left and right frame halves
            plot_center_coordinates(data=x_center[:, 0], x_y='x', left_right='left')
            plot_center_coordinates(data=y_center[:, 0], x_y='y', left_right='left')
            plot_center_coordinates(data=x_center[:, 1], x_y='x', left_right='right')
            plot_center_coordinates(data=y_center[:, 1], x_y='y', left_right='right') 
        
            # Compute standard deviation of determined center coordinates
            x_center_std = np.std(x_center, ddof=1, axis=0)
            y_center_std = np.std(y_center, ddof=1, axis=0)
    
            # Print mean center coordinates with error
            separator = ' '*max([len('%.2f' % x) for x in np.append(x_center_std, y_center_std)])
            printandlog('\nMean center coordinates (pixels):')
            printandlog('x_left' + separator + '         y_left' + separator + '         x_right' + separator + '         y_right')
            printandlog('%.2f +/- %.2f    %.2f +/- %.2f    %.2f +/- %.2f    %.2f +/- %.2f' 
                        % (x_center_mean[0], x_center_std[0], y_center_mean[0], y_center_std[0], x_center_mean[1], x_center_std[1], y_center_mean[1], y_center_std[1]))    
            
            # Print warning if there is significant deviation among the center coordinates found
            if any(np.append(x_center_std, y_center_std) > 0.5):
                printandlog('\nWARNING, for at least one of the fitted center coordinates the standard deviation is larger than 0.5 pixels.')
        else:
            # Print mean center coordinates without error
            printandlog('\nCenter coordinates (pixels):')
            printandlog('x_left    y_left    x_right    y_right')
            printandlog('%.2f    %.2f    %.2f    %.2f' % (x_center_mean[0], y_center_mean[0], x_center_mean[1], y_center_mean[1]))   

    if plot_centering_sub_images == True:
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
        for i in range(number_figures):
            if number_figures == 1:
                plot_name = name_file_root + 'centering_sub_images.png'            
            else:
                plot_name = name_file_root + 'centering_sub_images_' + str(i + 1) + '.png'            
            path_plot = os.path.join(path_plots_dir, plot_name)
            printandlog('\nCreating plot ' + path_plot + ' showing the sub-images of each frame and the center coordinates fitted.')
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
            plt.show()
    
    return cube_single_sum, cube_single_difference, header

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

def subtract_background(cube, param_annulus_background):   
    '''
    Subtract background from cube or frame
     
    Input:
        cube: image cube or frame to subtract background from
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
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
# process_flux_frames
###############################################################################

def process_flux_frames(path_flux_files, 
                        file_index_flux, 
                        indices_to_remove_flux, 
                        frame_master_flat, frame_master_bpm, 
                        frame_master_sky_flux, 
                        param_annulus_background, 
                        sigmafiltering=True, 
                        centering_method='gaussian', 
                        center_coordinates=(444, 490, 1469, 479), 
                        param_centering=(12, None, 30000), 
                        collapse_ndit=False, 
                        plot_centering_sub_images=True):
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
        param_annulus_background: (list of) length-6-tuple(s) with parameters 
            to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and
                rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and 
                rotating counterclockwise)
        sigmafiltering: if True remove bad pixels remaining after applying
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
            is (444, 490, 1469, 479).         
        param_centering: length-3-tuple with parameters for centering:
            crop_radius: half the length of the sides of the square cropped 
                sub-images used to fit the 2D Gaussian to (pixels). Must be 
                integer. The sub-image is centered on the coordinates as 
                provided by center_coordinates. If None, the complete frame is 
                used for the fitting and center_coordinates is ignored. The
                value of crop_radius is also used to create the sub-images when
                plot_centering_sub_images = True.
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
        collapse_ndit: If True, compute the mean over the (NDIT) frames of a
            file before subtracting the background, flat-fielding, bad pixel
            removal and centering. If False, perform the above steps for each
            frame and after that compute the mean over the frames 
            (default = False).
        plot_centering_sub_images: If True, plot the sub-images showing the 
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
    cube_single_sum = process_object_frames(path_object_files=path_flux_files, 
                                            file_index_object=file_index_flux, 
                                            indices_to_remove_object=indices_to_remove_flux, 
                                            frame_master_flat=frame_master_flat, 
                                            frame_master_bpm=frame_master_bpm, 
                                            frame_master_sky=frame_master_sky_flux, 
                                            sigmafiltering=sigmafiltering, 
                                            centering_method=centering_method, 
                                            center_coordinates=center_coordinates, 
                                            param_centering=param_centering, 
                                            collapse_ndit=collapse_ndit, 
                                            plot_centering_sub_images=plot_centering_sub_images)[0]

    # Compute flux frame as mean of single sum cube
    frame_flux = np.mean(cube_single_sum, axis=0)

    # Determine background and subtract it
    frame_master_flux, background = subtract_background(cube=frame_flux, 
                                                        param_annulus_background=param_annulus_background)

    printandlog('\nSubtracted background in master flux image = %.3f' % background)
      
    # Create frame showing annulus used to determine background
    frame_annulus_background = compute_annulus_values(cube=frame_flux, param=param_annulus_background)[1]
    
    # Print number of frames used to create master flux frame
    number_frames_total = sum([pyfits.getheader(x)['ESO DET NDIT'] for x in path_flux_files])
    number_frames_removed = sum([len(x) for x in indices_to_remove_flux])
    number_frames_used = number_frames_total - number_frames_removed
    
    printandlog('\nThe master flux frame was created out of ' + str(len(path_flux_files)) + 
            ' raw FLUX-file(s) comprising a total of ' + str(number_frames_used) + ' frame(s).')

    return frame_master_flux, frame_annulus_background

###############################################################################
# perform_preprocessing
###############################################################################

def perform_preprocessing(frames_to_remove=[], 
                          sigmafiltering=True, 
                          collapse_ndit_object=True, 
                          plot_centering_sub_images=True, 
                          centering_method_object='center frames', 
                          centering_subtract_object=True, 
                          center_coordinates_object=(477, 521, 1503, 511), 
                          param_centering_object=(12, 7, 30000), 
                          centering_method_flux='gaussian', 
                          center_coordinates_flux=(444, 490, 1469, 479), 
                          param_centering_flux=(12, None, 30000), 
                          param_annulus_background_flux='large annulus',
                          save_preprocessed_data=True, 
                          combination_method_polarization_images='trimmed mean'):
    '''
    Perform pre-processing of OBJECT, CENTER, SKY and FLUX-files, i.e. sorting data, 
    background subtraction, flat-fielding, bad pixel removal, centering and compution
    of the single-sum and -difference images used in the post-processing

    Input:
        frames_to_remove: list of integers and length-2-tuples of integers 
            indicating which files and frames to remove (0-based). A complete
            file can be removed by specifying its integer index, while a frame
            of specific file can be removed by specifying a tuple
            (file_index, frame_index). If no files or frames should be removed, 
            use an empty list [] (default = []). The files are sorted in
            chronological order from oldest to newest.
        sigmafiltering: if True, remove bad pixels remaining after applying
            master bad pixel map using sigma-filtering (default = True). Applies
            to all file-types (OBJECT, CENTER, SKY and FLUX).
        collapse_ndit_object: If True, compute the mean over the (NDIT) frames of
            the OBJECT-files before subtracting the background, flat-fielding, bad 
            pixel removal and centering to speed up the preprocessing. If False, 
            perform the above steps for each frame and after that compute the 
            mean over the frames (default = True).
        plot_centering_sub_images: If True, plot the sub-images showing the 
            center coordinates for each frame of the OBJECT- and FLUX-files. 
            The plots allow for checking whether the centering is correct and 
            to scan the data for frames with bad quality (default = True).    
        centering_method_object: method to center the OBJECT-frames. If 
            'center frames' or 'manual', use fixed coordinates as provided by 
            center_coordinates_object. If 'gaussian', fit a 2D Gaussian to each frame. 
            If 'cross-correlation', fit a 2D Gaussian to the first frame and then
            use cross-correlation to align (register) the other frames onto the 
            centered first  frame. For 'gaussian' and 'cross-correlation' 
            center_coordinates_object is used as initial guess of the center 
            coordinates and the determined center coordinates are plotted for
            each image (default = 'center frames').
        centering_subtract_object: if True subtract the OBJECT-file(s) taken
            closest in time from the CENTER-file(s) (default = True). This generally
            results in a more accurate determination of the center coordinates
            as the background and any other celestial objects in the field of view
            are suppressed. When the difference between the image orientations of 
            the CENTER- and OBJECT-frames is large (i.e. difference in derotator 
            position angle for field-tracking and difference in parallactic angle 
            for pupil-tracking), the OBJECT-frame will be rotated around the initial 
            guess of the centers as defined by center_coordinates_object before 
            subtracting it from the CENTER-file. 
        center_coordinates_object: length-4-tuple with center coordinates of OBJECT-frames:
            x_left: x-coordinate of center of left frame half
            y_left: y-coordinate of center of left frame half
            x_right: x-coordinate of center of right frame half
            y_right: y-coordinate of center of right frame half
            Note that the center coordinates are defined in the complete frame, 
            i.e. with both detector halves (pixels; 0-based). The default value 
            is (477, 521, 1503, 511). 
        param_centering_object: length-3-tuple with parameters for centering of OBJECT-frames:
            crop_radius: half the length of the sides of the square cropped 
                sub-images used to fit the 2D Gaussian to and used for 
                cross-correlating the images (pixels). Must be integer. The 
                sub-image is centered on the coordinates as provided by
                center_coordinates. If None, the complete frame is used for the 
                fitting and center_coordinates is ignored. The value of
                crop_radius is also used to create the sub-images when
                plot_centering_sub_images = True.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of param_centering_object is (12, 7, 30000).            
        centering_method_flux: method to center the FLUX-frames. If 'manual', use fixed 
            coordinates as provided by center_coordinates_flux. If 'gaussian', fit 
            a 2D Gaussian to each frame. For 'gaussian' center_coordinates_flux is 
            used as initial guess of the center coordinates and the determined 
            center coordinates are plotted for each image (default = 'gaussian').
        center_coordinates_flux: length-4-tuple with center coordinates of FLUX-frames:
            x_left: x-coordinate of center of left frame half
            y_left: y-coordinate of center of left frame half
            x_right: x-coordinate of center of right frame half
            y_right: y-coordinate of center of right frame half
            Note that the center coordinates are defined in the complete frame, 
            i.e. with both detector halves (pixels; 0-based). The default value 
            is (444, 490, 1469, 479).         
        param_centering_flux: length-3-tuple with parameters for centering if FLUX-frames:
            crop_radius: half the length of the sides of the square cropped 
                sub-images used to fit the 2D Gaussian to (pixels). Must be 
                integer. The sub-image is centered on the coordinates as 
                provided by center_coordinates. If None, the complete frame is 
                used for the fitting and center_coordinates is ignored. The
                value of crop_radius is also used to create the sub-images when
                plot_centering_sub_images = True.
            sigfactor: all sub-image pixels with values smaller than 
                sigfactor*standard deviation are replaced by random Gaussian noise 
                to mask them for fitting the 2D Gaussian. If None, no pixels are
                replaced by Gaussian noise.
            saturation_level: all pixels within the smallest circle encompassing 
                the pixels with a value equal to or higher than saturation_level 
                are ignored when fitting the 2D Gaussian. We use a circle because
                strongly saturated pixels in the peak of the PSF often have values 
                lower than saturation_level. If None, no pixels are ignored.
            The default value of param_centering_flux is (12, None, 30000).            
        param_annulus_background_flux: (list of) length-6-tuple(s) with parameters 
            to generate annulus to measure and subtract background in master flux frame:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and
                rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and 
                rotating counterclockwise)
            If string 'large annulus' the annulus will be star-centered and 
            located far away from the star with an inner radius of 320 pixels
            and a width of 60 pixels (default = 'large annulus').
        save_preprocessed_data: If True, save preprocessed cubes of single-sum 
            and single-difference images in the 'preprocessed' folder so that 
            the preprocessing can be skipped when re-running the pipeline
            (default = True). 
        combination_method_polarization_images: method to be used to produce the 
            incident Q- and U-images, 'least squares', 'trimmed mean' or 'median' 
            (default = 'trimmed mean')
        
    Output:
        cube_single_sum: cube of single-sum I_Q^+, I_Q^-, I_U^+ and I_U^-intensity images
        cube_single_difference: cube of single-difference Q^+, Q^-, U^+ and U^-images
        header: list of FITS-headers of raw science frames 
        file_index_object: list of file indices of OBJECT-files (0-based)
        combination_method_polarization_images: method to be used to produce the 
            incident Q- and U-images, 'least squares', 'trimmed mean' or 'median'
            
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
    path_object_files, path_sky_files, path_center_files, path_object_center_files, path_flux_files, path_sky_flux_files, \
    indices_to_remove_object, indices_to_remove_sky, indices_to_remove_center, indices_to_remove_object_center, indices_to_remove_flux, \
    indices_to_remove_sky_flux, file_index_object, file_index_flux, combination_method_polarization_images \
    = check_sort_data_create_directories(frames_to_remove=frames_to_remove, 
                                         combination_method_polarization_images=combination_method_polarization_images, 
                                         centering_method_object=centering_method_object, 
                                         save_preprocessed_data=save_preprocessed_data, 
                                         plot_centering_sub_images=plot_centering_sub_images)

    ###############################################################################
    # Read static master flat and bad pixel map
    ###############################################################################
    
    #TODO: Make master flats for Y and Ks. 
    # Determine filter used
    filter_used = pyfits.getheader(path_object_files[0])['ESO INS1 FILT ID']

    # Read static master flat    
    if filter_used == 'FILT_BBF_Y':
        path_static_flat = os.path.join(path_static_flat_badpixelmap, 'masterflat_Y.fits')
    elif filter_used == 'FILT_BBF_J':
        path_static_flat = os.path.join(path_static_flat_badpixelmap, 'masterflat_J.fits')
    elif filter_used == 'FILT_BBF_H':
        path_static_flat = os.path.join(path_static_flat_badpixelmap, 'masterflat_H.fits')
    elif filter_used == 'FILT_BBF_Ks':
        path_static_flat = os.path.join(path_static_flat_badpixelmap, 'masterflat_Ks.fits')
    
    frame_master_flat = np.squeeze(read_fits_files(path=path_static_flat, silent=True)[0])
    
    #TODO: Perform the the 2 lines below on the master flats and save them again so that these lines can be removed from the code
    # Filter master flat for zeros and NaN's
    frame_master_flat = np.nan_to_num(frame_master_flat)
    frame_master_flat[frame_master_flat == 0] = 1

    # Read static bad pixel map
    frame_master_bpm = np.squeeze(read_fits_files(path=os.path.join(path_static_flat_badpixelmap, 'master_badpix.fits'), silent=True)[0])

    #TODO: Here read static master dark,background frames (in all filters) in case there are no sky frames
    
    #TODO: Somewhere we should add the possibility to create master flats and master dark,backgrounds yourself
    
    ###############################################################################
    # Computing master sky for object images
    ###############################################################################
    
    if any(path_sky_files):
        # Process the sky files for the object files
        printandlog('\n###############################################################################')
        printandlog('# Processing SKY-files for OBJECT-files')
        printandlog('###############################################################################') 

        frame_master_sky = process_sky_frames(path_sky_files=path_sky_files, 
                                              indices_to_remove_sky=indices_to_remove_sky, 
                                              frame_master_bpm=frame_master_bpm, 
                                              sigmafiltering=sigmafiltering)
        
        # Write master sky-frame
        printandlog('')
        write_fits_files(data=frame_master_sky, path=os.path.join(path_sky_dir, name_file_root + 'master_sky.fits'), header=False, silent=False)
    
    else:
        # Create a master sky frame with only zeros
        frame_master_sky = np.zeros((1024, 2048))

    ###############################################################################
    # Processing center files and extracting center coordinates
    ###############################################################################
    
    if centering_method_object == 'center frames':
        # Print that we process the center files
        printandlog('\n###############################################################################')
        printandlog('# Processing CENTER-files')
        printandlog('###############################################################################')       
        
#        # Define or read centering parameters and print them on screen
#        if param_centering == 'default':
#            printandlog('\nUsing the default centering parameters:')
#            param_centering = (477, 521, 1503, 511, 12, 7, 30000)
#        else:
#            printandlog('\nUsing user-defined centering parameters:')
#        printandlog(param_centering)
#TODO: Remove 3 lines below
#        separator_1 = ' '*(14 - len('%s'% param_centering[4]))
#        printandlog('x_left    y_left    x_right    y_right    crop_radius   sigfactor')
#        printandlog('%.2f    %.2f    %.2f    %.2f     %s    ' % param_centering[:5] + separator_1 + '%.1f' % param_centering[5])
            
        # Process the center files            
        list_frame_center_processed, header_center = process_center_frames(path_center_files=path_center_files, 
                                                                           indices_to_remove_center=indices_to_remove_center, 
                                                                           path_object_center_files=path_object_center_files, 
                                                                           indices_to_remove_object_center=indices_to_remove_object_center, 
                                                                           frame_master_flat=frame_master_flat, 
                                                                           frame_master_bpm=frame_master_bpm, 
                                                                           frame_master_sky=frame_master_sky, 
                                                                           centering_subtract_object=centering_subtract_object, 
                                                                           center_coordinates=center_coordinates_object, 
                                                                           sigmafiltering=sigmafiltering)
    
        # Write processed center frames
        path_processed_center_files = [os.path.join(path_center_dir, os.path.splitext(os.path.basename(x))[0] + '_processed.fits') for x in path_center_files]
        printandlog('')
        write_fits_files(data=list_frame_center_processed, path=path_processed_center_files, header=header_center, silent=False)
        
        # Find center coordinates and replace values of center_coordinates
        center_coordinates_object = find_center_coordinates(list_frame_center_processed=list_frame_center_processed, 
                                                            path_processed_center_files=path_processed_center_files, 
                                                            center_coordinates=center_coordinates_object, 
                                                            param_centering=param_centering_object)

    ###############################################################################
    # Creating processed and centered single-sum and -difference images
    ###############################################################################
    
    # Create reduced and centerd single-sum and -difference images
    printandlog('\n###############################################################################')
    printandlog('# Processing OBJECT-files')
    printandlog('###############################################################################') 
      
    cube_single_sum, cube_single_difference, header = process_object_frames(path_object_files=path_object_files, 
                                                                            file_index_object=file_index_object, 
                                                                            indices_to_remove_object=indices_to_remove_object, 
                                                                            frame_master_flat=frame_master_flat, 
                                                                            frame_master_bpm=frame_master_bpm, 
                                                                            frame_master_sky=frame_master_sky, 
                                                                            sigmafiltering=sigmafiltering, 
                                                                            centering_method=centering_method_object, 
                                                                            center_coordinates=center_coordinates_object, 
                                                                            param_centering=param_centering_object, 
                                                                            collapse_ndit=collapse_ndit_object, 
                                                                            plot_centering_sub_images=plot_centering_sub_images)
    
    if save_preprocessed_data == True:
        # Write preprocessed cubes of single-sum and single-difference images 
        printandlog('\nSaving pre-processed data so that pre-processing can be skipped the next time.')
        printandlog('')
        write_fits_files(data=cube_single_sum, path=os.path.join(path_preprocessed_dir, 'cube_single_sum.fits'), header=False, silent=False)
        write_fits_files(data=cube_single_difference, path=os.path.join(path_preprocessed_dir, 'cube_single_difference.fits'), header=False, silent=False)

        # Write path of object files to a .txt-file to be able to read headers
        with open(os.path.join(path_preprocessed_dir, 'path_object_files.txt'), 'w') as fh:
            for path_sel in path_object_files:
                fh.write('%s\n' % path_sel)
        printandlog('Wrote file ' + os.path.join(path_preprocessed_dir, 'path_object_files.txt') + '.', wrap=False)

    ###############################################################################
    # Computing master sky for flux images
    ###############################################################################
    
    if any(path_sky_flux_files):
        # Process the sky files for the flux files
        printandlog('\n###############################################################################')
        printandlog('# Processing SKY-files for FLUX-files')
        printandlog('###############################################################################') 

        frame_master_sky_flux = process_sky_frames(path_sky_files=path_sky_flux_files, 
                                                   indices_to_remove_sky=indices_to_remove_sky_flux, 
                                                   frame_master_bpm=frame_master_bpm, 
                                                   sigmafiltering=sigmafiltering)
        
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
        if type(param_annulus_background_flux) == tuple or type(param_annulus_background_flux) == list:
            printandlog('\nThe background will be determined with a user-defined annulus or several user-defined annuli:')
            if type(param_annulus_background_flux) == tuple:
                printandlog(param_annulus_background_flux)
            elif type(param_annulus_background_flux) == list:
                for x in param_annulus_background_flux:
                    printandlog(x)
        elif param_annulus_background_flux == 'large annulus':
            param_annulus_background_flux = (511.5, 511.5, 320, 60, 0, 360)
            printandlog('\nThe background will be determined with a star-centered annulus located far away from the star:')
            printandlog(param_annulus_background_flux)

        # Processthe flux files
        frame_master_flux, frame_annulus_background_flux = process_flux_frames(path_flux_files=path_flux_files, 
                                                                               file_index_flux=file_index_flux, 
                                                                               indices_to_remove_flux=indices_to_remove_flux, 
                                                                               frame_master_flat=frame_master_flat, 
                                                                               frame_master_bpm=frame_master_bpm, 
                                                                               frame_master_sky_flux=frame_master_sky_flux, 
                                                                               param_annulus_background=param_annulus_background_flux, 
                                                                               sigmafiltering=sigmafiltering, 
                                                                               centering_method=centering_method_flux, 
                                                                               center_coordinates=center_coordinates_flux, 
                                                                               param_centering=param_centering_flux, 
                                                                               collapse_ndit=False, 
                                                                               plot_centering_sub_images=plot_centering_sub_images)

        # Write master flux-frame and frame showing annulus used to determine background
        printandlog('')
        write_fits_files(data=frame_master_flux, path=os.path.join(path_flux_dir, name_file_root + 'master_flux.fits'), header=False, silent=False)
        write_fits_files(data=frame_annulus_background_flux, path=os.path.join(path_flux_dir, name_file_root + 'annulus_background_flux.fits'), header=False)    

#TODO: conversion of final images to mJansky/arcsec^2 should be part of post-processing part and optional. Also add possibility to express as contrast wrt central star? 
        
    return cube_single_sum, cube_single_difference, header, file_index_object, combination_method_polarization_images

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
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
              
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
    I_Q = np.mean(compute_annulus_values(cube=cube_I_Q, param=param_annulus_star)[0], axis=mean_median_axis) - \
                  np.median(compute_annulus_values(cube=cube_I_Q, param=param_annulus_background)[0], axis=mean_median_axis)       
    I_U = np.mean(compute_annulus_values(cube=cube_I_U, param=param_annulus_star)[0], axis=mean_median_axis) - \
                  np.median(compute_annulus_values(cube=cube_I_U, param=param_annulus_background)[0], axis=mean_median_axis)      
    Q = np.mean(compute_annulus_values(cube=cube_Q, param=param_annulus_star)[0], axis=mean_median_axis) - \
                np.median(compute_annulus_values(cube=cube_Q, param=param_annulus_background)[0], axis=mean_median_axis)       
    U = np.mean(compute_annulus_values(cube=cube_U, param=param_annulus_star)[0], axis=mean_median_axis) - \
                np.median(compute_annulus_values(cube=cube_U, param=param_annulus_background)[0], axis=mean_median_axis)    
    
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

def correct_instrumental_polarization_effects(cube_I_Q_double_sum, 
                                              cube_I_U_double_sum, 
                                              cube_Q_double_difference, 
                                              cube_U_double_difference, 
                                              header, 
                                              file_index_object, 
                                              param_annulus_star, 
                                              param_annulus_background, 
                                              combination_method_polarization_images='trimmed mean', 
                                              trimmed_mean_proportiontocut_polarization_images=0.1, 
                                              combination_method_total_intensity_images='trimmed mean', 
                                              trimmed_mean_proportiontocut_total_intensity_images=0.1, 
                                              images_north_up=True):
    '''
    Calculate incident I_Q-, I_U-, Q- and U-images by correcting for the instrumental polarization effects of IRDIS using the polarimetric instrument model

    Input:
        cube_I_Q_double_sum: cube of double-sum intensity I_Q-images in order of HWP cycles
        cube_I_U_double_sum: cube of double-sum intensity I_U-images in order of HWP cycles  
        cube_Q_double_difference: cube of double-difference Stokes Q-images in order of HWP cycles
        cube_U_double_difference: cube of double-difference Stokes U-images in order of HWP cycles 
        header: list of FITS-headers of OBJECT-files
        file_index_object: list of file indices of OBJECT-files (0-based)
        param_annulus_star: (list of) length-6-tuple(s) with parameters to generate annulus to measure polarization of star:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
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
    plt.show()     
    
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

    if tracking_mode_used == 'FIELD' and number_derotator_position_angles == 1 and images_north_up == False:
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
    if np.all(np.array(param_annulus_star, ndmin=2)[:, np.array([0, 1, 4, 5])] == np.array([511.5, 511.5, 0, 360])) and \
       np.all(np.array(param_annulus_background, ndmin=2)[:, np.array([0, 1, 4, 5])] == np.array([511.5, 511.5, 0, 360])):
        # Do not rotate the cubes of double-sum and double-difference images and compute normalized Stokes q and u in an annulus   
        printandlog('\nNot rotating the images used to determine the polarization signal in an annulus because the annuli used for the star and the background are centered on the star and rotationally symmetric.')

        q_annulus, u_annulus = determine_star_polarization(cube_I_Q=cube_I_Q_double_sum, 
                                                           cube_I_U=cube_I_U_double_sum, 
                                                           cube_Q=cube_Q_double_difference, 
                                                           cube_U=cube_U_double_difference, 
                                                           param_annulus_star=param_annulus_star, 
                                                           param_annulus_background=param_annulus_background)

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
                                                           param_annulus_star=param_annulus_star, 
                                                           param_annulus_background=param_annulus_background)
    
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
    plt.show()

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
    plt.show()    
      
    ###############################################################################
    # Compute incident Q- and U-images by correcting for instrumental polarization effects
    ###############################################################################

    # Subtract intrumental polarization from Q- and U-images
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
    if combination_method_polarization_images == 'least squares':
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
    
    elif combination_method_polarization_images == 'trimmed mean': 
        # Compute incident Q- and U-images from the trimmed mean of incident cubes
        printandlog('\nComputing the incident Q- and U-images using the trimmed mean with a proportion to cut equal to ' + str(trimmed_mean_proportiontocut_polarization_images) + '.')
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
        printandlog('\nComputing the incident I_Q- and I_U-images using the trimmed mean with a proportion to cut equal to ' + str(trimmed_mean_proportiontocut_total_intensity_images) + '.')
        frame_I_Q_incident = trim_mean(cube_I_Q_incident, proportiontocut=trimmed_mean_proportiontocut_total_intensity_images, axis=0)
        frame_I_U_incident = trim_mean(cube_I_U_incident, proportiontocut=trimmed_mean_proportiontocut_total_intensity_images, axis=0)
    
    elif combination_method_total_intensity_images == 'median': 
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

def perform_postprocessing(cube_single_sum, 
                           cube_single_difference, 
                           header, 
                           file_index_object, 
                           param_annulus_star='ao residuals', 
                           param_annulus_background='large annulus', 
                           double_difference_type='standard', 
                           remove_vertical_band_detector_artefact=True, 
                           combination_method_polarization_images='trimmed mean', 
                           trimmed_mean_proportiontocut_polarization_images=0.1, 
                           combination_method_total_intensity_images='trimmed mean', 
                           trimmed_mean_proportiontocut_total_intensity_images=0.1, 
                           images_north_up=True, 
                           create_images_DoLP_AoLP_q_u_norm=False):
    '''
    Perform post-processing of data, including applying the model-based correction
    for the instrumental polarization effects, and save final images to FITS-files
    
    Input:
        cube_single_sum: cube of pre-processed single-sum images
        cube_single_difference: cube of pre-processed single-difference images
        header: list of FITS-headers of OBJECT-files 
        file_index_object: list of file indices of OBJECT-files (0-based)
        param_annulus_star: (list of) length-6-tuple(s) with parameters to generate annulus to measure polarization of star:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            If string 'ao residuals' the annulus will automatically determined and 
            be star-centered and located over the AO residuals. The inner radius 
            and width of the annulus will depend on the filter used. If 
            'star aperture' a small aparture located at the position of
            the central star will be used (default = 'ao residuals').
        param_annulus_background: (list of) length-6-tuple(s) with parameters to generate annulus to measure and subtract background:
            coord_center_x: x-coordinate of center (pixels; 0-based)
            coord_center_y: y-coordinate of center (pixels; 0-based)
            inner_radius: inner radius (pixels)
            width: width (pixels)
            start_angle: start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            end_angle: end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
            If string 'large annulus' the annulus will be star-centered and 
            located far away from the star with an inner radius of 360 pixels
            and a width of 60 pixels (default = 'large annulus').
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
                
    # Determine filter and coronagraph used    
    filter_used = header[0]['ESO INS1 FILT ID']   
    coronagraph_used = header[0]['ESO INS COMB ICOR']

    # Define and print annulus to determine the star polarization from
    if type(param_annulus_star) == tuple or type(param_annulus_star) == list:
        printandlog('\nThe star polarization will be determined with a user-defined annulus or several user-defined annuli:')
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
        printandlog('\nThe star polarization will be determined with a star-centered annulus located over the AO residuals:')
        printandlog(param_annulus_star)
        if coronagraph_used == 'N_NS_CLEAR':
            printandlog('\nWARNING, the data is non-coronagraphic so there might be little flux at the AO residuals. Determining the star polarization using an aperture at the position of the central star (\'star aperture\') will probably yield better results.')
    elif param_annulus_star == 'star aperture':
        param_annulus_star = (511.5, 511.5, 0, 11, 0, 360)
        printandlog('\nThe star polarization will be determined with an aparture located at the position of the central star:')
        printandlog(param_annulus_star)
    
    # Define and print annulus to determine the background from
    if type(param_annulus_background) == tuple or type(param_annulus_background) == list:
        printandlog('\nThe background will be determined with a user-defined annulus or several user-defined annuli:')
        if type(param_annulus_background) == tuple:
            printandlog(param_annulus_background)
        elif type(param_annulus_background) == list:
            for x in param_annulus_background:
                printandlog(x)
    elif param_annulus_background == 'large annulus':
        param_annulus_background = (511.5, 511.5, 360, 60, 0, 360)
        printandlog('\nThe background will be determined with a star-centered annulus located far away from the central star:')
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
                                                param_annulus_star=param_annulus_star, 
                                                param_annulus_background=param_annulus_background, 
                                                combination_method_polarization_images=combination_method_polarization_images, 
                                                trimmed_mean_proportiontocut_polarization_images=trimmed_mean_proportiontocut_polarization_images, 
                                                combination_method_total_intensity_images=combination_method_total_intensity_images, 
                                                trimmed_mean_proportiontocut_total_intensity_images=trimmed_mean_proportiontocut_total_intensity_images, 
                                                images_north_up=images_north_up)
    
    ###############################################################################
    # Subtract background in I_Q-, I_U-, Q- and U-frames
    ###############################################################################
 
    printandlog('\n###############################################################################')
    printandlog('# Subtracting the backgrounds and measuring and removing the star polarization')
    printandlog('###############################################################################') 

    # Determine background in corrected I_Q-, I_U-, Q- and U-frames and subtract it
    frame_I_Q_background_subtracted, background_frame_I_Q = subtract_background(cube=frame_I_Q_incident, param_annulus_background=param_annulus_background)
    frame_I_U_background_subtracted, background_frame_I_U = subtract_background(cube=frame_I_U_incident, param_annulus_background=param_annulus_background)
    frame_Q_background_subtracted, background_frame_Q = subtract_background(cube=frame_Q_incident, param_annulus_background=param_annulus_background)
    frame_U_background_subtracted, background_frame_U = subtract_background(cube=frame_U_incident, param_annulus_background=param_annulus_background)
   
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
                                                 param_annulus_star=param_annulus_star, 
                                                 param_annulus_background=param_annulus_background)
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
                                                                                                                param_annulus_background=param_annulus_background)
    frame_U_star_polarization_subtracted, background_frame_U_star_polarization_subtracted = subtract_background(cube=frame_U_star_polarization_subtracted, 
                                                                                                                param_annulus_background=param_annulus_background)
    
    # Print resulting background values
    printandlog('\nSubtracted residual backgrounds in the star-polarization-subtracted Q- and U-images:')
    printandlog('Background Q = %.3f' % background_frame_Q_star_polarization_subtracted)
    printandlog('Background U = %.3f' % background_frame_U_star_polarization_subtracted)

    if len(cube_I_Q_incident) > 1 and len(cube_I_U_incident) > 1:
        ###############################################################################
        # Subtract background in I_Q-, I_U-, Q- and U-cubes
        ###############################################################################
        
        # Determine background in corrected I_Q-, I_U-, Q- and U-frames and subtract it
        cube_I_Q_background_subtracted, background_cube_I_Q = subtract_background(cube=cube_I_Q_incident, param_annulus_background=param_annulus_background)
        cube_I_U_background_subtracted, background_cube_I_U = subtract_background(cube=cube_I_U_incident, param_annulus_background=param_annulus_background)
        cube_Q_background_subtracted, background_cube_Q = subtract_background(cube=cube_Q_incident, param_annulus_background=param_annulus_background)
        cube_U_background_subtracted, background_cube_U = subtract_background(cube=cube_U_incident, param_annulus_background=param_annulus_background)
       
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
                                                                         param_annulus_star=param_annulus_star, 
                                                                         param_annulus_background=param_annulus_background)
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
        plt.show()
        
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
        plt.show()
        
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
                           images_north_up=images_north_up)

    # Compute final images with the star polarization subtracted   
    frame_Q_phi_star_polarization_subtracted, frame_U_phi_star_polarization_subtracted, frame_I_pol_star_polarization_subtracted, \
    frame_AoLP_star_polarization_subtracted, frame_DoLP_star_polarization_subtracted, frame_q_star_polarization_subtracted, \
    frame_u_star_polarization_subtracted, frame_AoLP_norm_star_polarization_subtracted, frame_DoLP_norm_star_polarization_subtracted \
    = compute_final_images(frame_I_Q=frame_I_Q_background_subtracted, 
                           frame_I_U=frame_I_U_background_subtracted, 
                           frame_Q=frame_Q_star_polarization_subtracted, 
                           frame_U=frame_U_star_polarization_subtracted, 
                           header=header, 
                           images_north_up=images_north_up)[1:]

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
        printandlog('\nWARNING, the images DoLP.fits, q_norm.fits, u_norm.fits, AoLP_norm.fits and DoLP_norm.fits are only valid if all flux in the images originates from the astrophysical source of interest. This is generally the case for observations of for example solar system objects or galaxies. The images are generally not valid for observations of circumstellar disks or companions because in that case a large part of the flux in the total intensity images originates from the central star.')
        frames_to_write += [frame_DoLP, frame_q, frame_u, frame_AoLP_norm, frame_DoLP_norm]       
        file_names += ['DoLP', 'q_norm', 'u_norm', 'AoLP_norm', 'DoLP_norm']
    
    # Write files of the images with the star polarization present
    printandlog('')
    for frame, file_name in zip(frames_to_write, file_names):
        write_fits_files(data=frame, path=os.path.join(path_reduced_dir, name_file_root + file_name + '.fits'), header=False)

    # Write frames that show annuli used to retrieve star and background signals in reduced directory
    write_fits_files(data=frame_annulus_star, path=os.path.join(path_reduced_dir, name_file_root + 'annulus_star.fits'), header=False)
    write_fits_files(data=frame_annulus_background, path=os.path.join(path_reduced_dir, name_file_root + 'annulus_background.fits'), header=False)    

    # List files of the images with the star polarization subtracted and define their file names
    frames_to_write = [frame_I_Q_background_subtracted, frame_I_U_background_subtracted, frame_I_tot, frame_Q_star_polarization_subtracted, frame_U_star_polarization_subtracted, \
                       frame_Q_phi_star_polarization_subtracted, frame_U_phi_star_polarization_subtracted, frame_I_pol_star_polarization_subtracted, frame_AoLP_star_polarization_subtracted]
    file_names = ['I_Q', 'I_U', 'I_tot', 'Q_star_pol_subtr', 'U_star_pol_subtr', 'Q_phi_star_pol_subtr', 'U_phi_star_pol_subtr', 'I_pol_star_pol_subtr', 'AoLP_star_pol_subtr']

    if create_images_DoLP_AoLP_q_u_norm == True:
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

###############################################################################
###############################################################################
## Run pipeline
###############################################################################
###############################################################################

# Start taking time
time_start = time.time()

###############################################################################
# Check whether input values are valid
###############################################################################

if type(path_main_dir) is not str:
    raise TypeError('\'path_main_dir\' should be of type string.')

if type(path_static_flat_badpixelmap) is not str:
    raise TypeError('\'path_static_flat_badpixelmap\' should be of type string.')

if double_difference_type not in ['standard', 'normalized']:
    raise ValueError('\'double_difference_type\' should be either \'standard\' or \'normalized\'.')
    
if remove_vertical_band_detector_artefact not in [True, False]:
    raise ValueError('\'remove_vertical_band_detector_artefact\' should be either True or False.')   

#TODO: Add checks of specific values of param_annulus_star/background
#      Check how an annulus behaves when it is cut off by the edge of the image
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
# Convert input from 1-based to 0-based indexing
###############################################################################

# center_coordinates_object
center_coordinates_object = tuple(x - 1 for x in center_coordinates_object)

# center_coordinates_flux
center_coordinates_flux = tuple(x - 1 for x in center_coordinates_flux)

# param_annulus_background_flux
if type(param_annulus_background_flux) is tuple:
    param_annulus_background_flux = (param_annulus_background_flux[0] - 1,) + (param_annulus_background_flux[1] - 1,) + param_annulus_background_flux[2:]
elif type(param_annulus_background_flux) is list:
    for i,x in enumerate(param_annulus_background_flux):
        x = (x[0] - 1,) + (x[1] - 1,) + x[2:]
        param_annulus_background_flux[i] = x

# frames_to_remove
for i,x in enumerate(frames_to_remove):
    if type(x) is tuple:
        x = tuple(y - 1 for y in x)
    else:
        x -= 1
    frames_to_remove[i] = x

# param_annulus_star
if type(param_annulus_star) is tuple:
    param_annulus_star = (param_annulus_star[0] - 1,) + (param_annulus_star[1] - 1,) + param_annulus_star[2:]
elif type(param_annulus_star) is list:
    for i,x in enumerate(param_annulus_star):
        x = (x[0] - 1,) + (x[1] - 1,) + x[2:]
        param_annulus_star[i] = x

# param_annulus_background
if type(param_annulus_background) is tuple:
    param_annulus_background = (param_annulus_background[0] - 1,) + (param_annulus_background[1] - 1,) + param_annulus_background[2:]
elif type(param_annulus_background) is list:
    for i,x in enumerate(param_annulus_background):
        x = (x[0] - 1,) + (x[1] - 1,) + x[2:]
        param_annulus_background[i] = x
      
###############################################################################
# Define global variables
###############################################################################

# Define pupil-offset (deg) in pupil-tracking mode (SPHERE User Manual P99.0, 6th public release, P99 Phase 1)
pupil_offset = 135.99
    
# Define true North correction (deg) in pupil-tracking mode (SPHERE User Manual P99.0, 6th public release, P99 Phase 1)
true_north_correction = -1.75
    
# Define mean solar day (s)
msd = 86400
   
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

# Check if raw directory exists, if not create it
if not os.path.exists(path_raw_dir):
    os.makedirs(path_raw_dir)
    raise IOError('The raw directory {0:s} did not exist. It was created but you need to put your raw FITS-files there.'.format(path_raw_dir))

# Define the base of the name of each file to be generated
header_first_file = pyfits.getheader(glob.glob(os.path.join(path_raw_dir,'*.fits'))[0])
target_name = header_first_file['ESO OBS TARG NAME']
date_obs = header_first_file['DATE-OBS']

name_file_root = target_name.replace(' ', '_') + '_' + date_obs[:10].replace(' ', '_') + '_'

###############################################################################
# Run pre-processing and post-processing functions
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
    printandlog('\nCreating log file ' + path_log_file + '.')
elif os.path.exists(path_log_file) == True and skip_preprocessing == False:
    # Empty existing log file and start writing
    open(path_log_file, 'w').close()
    printandlog('\n###############################################################################')
    printandlog('# Starting pre-processing')
    printandlog('###############################################################################') 
    printandlog('\nWARNING, the log file ' + path_log_file + ' will be overwritten.')
elif os.path.exists(path_log_file) == True and skip_preprocessing == True:
    # Remove all lines concerned with the post-processing from the log file
    log_file_lines = [x.rstrip('\n') for x in open(path_log_file, 'r')]
    log_file_lines = log_file_lines[:log_file_lines.index('# Starting post-processing') - 1]
    open(path_log_file, 'w').close()
    for line in log_file_lines:                            
        print(line, file=open(path_log_file, 'a'))
   
if skip_preprocessing == False:
    # Pre-process raw data
    cube_single_sum, cube_single_difference, header, file_index_object, combination_method_polarization_images \
    = perform_preprocessing(frames_to_remove=frames_to_remove, 
                            sigmafiltering=sigmafiltering, 
                            collapse_ndit_object=collapse_ndit_object, 
                            plot_centering_sub_images=plot_centering_sub_images, 
                            centering_method_object=centering_method_object, 
                            centering_subtract_object=centering_subtract_object, 
                            center_coordinates_object=center_coordinates_object, 
                            param_centering_object=param_centering_object, 
                            centering_method_flux=centering_method_flux, 
                            center_coordinates_flux=center_coordinates_flux, 
                            param_centering_flux=param_centering_flux, 
                            param_annulus_background_flux=param_annulus_background_flux, 
                            save_preprocessed_data=save_preprocessed_data, 
                            combination_method_polarization_images=combination_method_polarization_images)

    # Print that post-processing starts
    printandlog('\n###############################################################################')
    printandlog('# Starting post-processing')
    printandlog('###############################################################################') 
    printandlog('\nContinuing with the pre-processed data.')
                
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
        printandlog('\nSkipping pre-processing and reading pre-processed data and headers.')
        printandlog('')
        
        # Read pre-processed single-sum and difference- images        
        cube_single_sum = read_fits_files(path=path_cube_single_sum, silent=False)[0]
        cube_single_difference = read_fits_files(path=path_cube_single_difference, silent=False)[0] 

        # Read headers
        header = [pyfits.getheader(x.rstrip('\n')) for x in open(path_object_files_text, 'r')]
        printandlog('Read headers from OBJECT-files specified in ' + path_object_files_text + '.')
    
    else:
        raise IOError('The files ' + path_cube_single_sum + ', ' + path_cube_single_difference + ' and/or ' + path_object_files_text + ' do not exist. Set skip_preprocessing to False and save_preprocessed_data to True to perform the pre-processing of the raw data and save the results.')

# Perform post-processing of data           
perform_postprocessing(cube_single_sum=cube_single_sum, 
                       cube_single_difference=cube_single_difference, 
                       header=header, 
                       file_index_object=file_index_object, 
                       param_annulus_star=param_annulus_star, 
                       param_annulus_background=param_annulus_background, 
                       double_difference_type=double_difference_type, 
                       remove_vertical_band_detector_artefact=remove_vertical_band_detector_artefact, 
                       combination_method_polarization_images=combination_method_polarization_images, 
                       trimmed_mean_proportiontocut_polarization_images=trimmed_mean_proportiontocut_polarization_images, 
                       combination_method_total_intensity_images=combination_method_total_intensity_images, 
                       trimmed_mean_proportiontocut_total_intensity_images=trimmed_mean_proportiontocut_total_intensity_images,
                       images_north_up=images_north_up, 
                       create_images_DoLP_AoLP_q_u_norm=create_images_DoLP_AoLP_q_u_norm)

# Print time elapsed
time_end = time.time()
d = datetime.datetime(1, 1, 1) + datetime.timedelta(seconds = time_end - time_start)
printandlog('\nTime elapsed: %d h %d min %d s' % (d.hour, d.minute, d.second)) 

               





















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
