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
the Free Software Foundation, either version 3 of the License.

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
import argparse
from argparse import RawTextHelpFormatter
from .irdap import __version__
from .irdap import run_demo
from .irdap import create_overview_headers_main
from .irdap import make_config
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
                                                 'To run IRDAP, create a directory (e.g. "/home/T_Cha_2016-02-20") containing a\n' +
                                                 'subdirectory called "raw" in which you place the raw FITS-files. Then in the\n' +
                                                 'terminal navigate to the directory (e.g. "cd /home/T_Cha_2016-02-20") and type\n' +
                                                 '"irdap --makeconfig" to create a default configuration file "config.conf" in\n' +
                                                 'this directory. You can then adjust the parameters in the configuration file\n' +
                                                 'with a text editor. Finally in the terminal type "irdap --run" to perform the\n' +
                                                 'data reduction.\n\n' +
                                                 'If this is the first time you use IRDAP, it is recommended to run the demo\n' + 
                                                 'first by using the terminal to navigate to a directory of your choice and\n' +
                                                 'typing irdap --demo.\n\n' + 
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
    parser.add_argument('-d', '--demo', action='store_true',
                        help='run pipeline in current working directory with example\ndata of the circumstellar disk of T Cha')
    parser.add_argument('-o', '--headers', action='store_true',
                        help='create overview of relevant headers of FITS-files in raw\nsubdirectory')
    parser.add_argument('-c', '--makeconfig', action='store_true',
                        help='create default configuration file in current working\ndirectory')                    
    parser.add_argument('-r', '--run', action='store_true',
                        help='run pipeline using configuration file in current working\ndirectory')                    
    
    # Use current working directory (of terminal) as path of main directory of reduction    
    path_main_dir = os.getcwd()

    # Evaluate and act upon user arguments
    args = parser.parse_args()
 
    if args.demo:
        # Run example reduction
        run_demo(path_main_dir)
        
    elif args.headers:
        # Create an overview of relevant headers
        create_overview_headers_main(path_main_dir)
    
    elif args.makeconfig:
        # Create a default configuration file
        make_config(path_main_dir)
        
    elif args.run:
        # Run the pipeline 
        run_pipeline(path_main_dir)

###############################################################################
# Run the function main
###############################################################################

# Run function when called, i.e. in the terminal one can just write "irdap --run" i.o. "python -m irdap --run"
if __name__ == "__main__":
    main()