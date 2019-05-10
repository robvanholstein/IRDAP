'''
This file contains the top-level code of IRDAP that allows the user to execute 
it.

IRDAP is a Python package to accurately reduce SPHERE-IRDIS polarimetric data.
Copyright (C) 2019 R.G. van Holstein

Full documentation: http://www.spherepol.nl.
Feedback, questions, comments: vanholstein@strw.leidenuniv.nl.

When publishing data reduced with IRDAP you must cite van Holstein et al. 
(2019): <ADS link>. 
For data in pupil-tracking mode you must additionally cite van Holstein et al. 
(2017): http://adsabs.harvard.edu/abs/2017SPIE10400E..15V.
                                                 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''     
    
# Import packages
import sys
import os
import shutil
import argparse
import glob
import astropy.io.fits as pyfits
from argparse import RawTextHelpFormatter
from .irdap import __version__
from .irdap import create_overview_headers
from .irdap import print_wrap
from .irdap import input_wrap
from .irdap import run_pipeline

###############################################################################
# main
###############################################################################

def main(args=None):
    '''
    Main function to run IRDAP
    
    Input:
        args: user input arguments
    
    File written by Rob van Holstein
    Function status: verified       
    '''
    
    if args is None:
        # Obtain arguments that user put in
        args = sys.argv[1:]
        
    # Define the arser including the description and epilog
    parser = argparse.ArgumentParser(description='IRDAP (IRDIS Data-reduction for Accurate Polarimetry) is a pipeline for\n' +
                                                 'accurate reduction of SPHERE-IRDIS polarimetric data.\n\n' + 
                                                 'To run IRDAP, create a directory (e.g. "/home/T Cha 2016-02-20") containing a\n' +
                                                 'subdirectory called "raw" in which you place the raw FITS-files. Then in the\n' +
                                                 'shell navigate to the directory (e.g. "cd /home/T Cha 2016-02-20") and type\n' +
                                                 '"irdap --makeconfig" to create a default configuration file "config.conf" in\n' +
                                                 'this directory. You can then adjust the parameters in the configuration file\n' +
                                                 'with a text editor. Finally in the shell type "irdap --run" to perform the\n' +
                                                 'data reduction.\n\n' +
                                                 'When publishing data reduced with IRDAP you must cite van Holstein et al.\n' +
                                                 '(2019): <ADS link>.\n' +
                                                 'For data in pupil-tracking mode you must additionally cite van Holstein et al.\n' +
                                                 '(2017): http://adsabs.harvard.edu/abs/2017SPIE10400E..15V.',
                                     epilog='Full documentation: http://www.spherepol.nl.\n' +
                                             'Feedback, questions, comments: vanholstein@strw.leidenuniv.nl.\n\n' + 
                                             'IRDAP Copyright (C) 2019 R.G. van Holstein.', 
                                             formatter_class=RawTextHelpFormatter)
    
    # Add parser arguments
    parser.add_argument('-v', '--version', action='version', version=('IRDAP %s' % __version__))
    parser.add_argument('-o', '--headers', action='store_true',
                        help='create overview of relevant headers of FITS-files in raw\nsubdirectory')
    parser.add_argument('-c', '--makeconfig', action='store_true',
                        help='create default configuration file in current working\ndirectory')                    
    parser.add_argument('-r', '--run', action='store_true',
                        help='run pipeline using configuration file in current working\ndirectory')                    
    
    # Use current working directory (of shell) as path of main directory of reduction    
    path_main_dir = os.getcwd()

    # Evaluate and act upon user arguments
    args = parser.parse_args()
 
    if args.headers:
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
                print_wrap('\nCreated an overview of the headers ' + path_overview + '.')
    
    elif args.makeconfig:
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

    elif args.run:
        # Run the pipeline
        run_pipeline(path_main_dir)

###############################################################################
# Run the function main
###############################################################################

# Run function when called, i.e. in the shell one can just write "irdap --run" i.o. "python -m irdap --run"
if __name__ == "__main__":
    main()