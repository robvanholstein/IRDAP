'''
This file contains the top-level code of IRDAP that allows the user to execute
it.

IRDAP is a Python package to accurately reduce SPHERE-IRDIS polarimetric data.
Copyright (C) 2019 R.G. van Holstein

Full documentation: https://irdap.readthedocs.io
Feedback, questions, comments: rob.vanholstein@eso.org

When publishing data reduced with IRDAP, please cite van Holstein et al.
(2020): https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..64V/abstract.
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

# Import packages
import sys
import os
import argparse
import webbrowser
import urllib
from argparse import RawTextHelpFormatter
from .version import __version__
from .irdap import run_demo
from .irdap import create_overview_headers_main
from .irdap import make_config
from .irdap import run_pipeline
from .irdap import mean_combine_images

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

    # Check if at least one argument is given
    if len(args) == 0:
        print('\nNo arguments were provided. Please check the help message by typing\n"irdap --help".')

    # Define the arser including the description and epilog
    parser = argparse.ArgumentParser(description='IRDAP (IRDIS Data reduction for Accurate Polarimetry) is a pipeline for\n' +
                                                 'accurate reduction of SPHERE-IRDIS polarimetric data.\n\n' +
                                                 'To run IRDAP, create a directory (e.g. "/home/T_Cha_2016-02-20") containing a\n' +
                                                 'subdirectory called "raw" in which you place the raw FITS-files. Then in the\n' +
                                                 'terminal navigate to the directory (e.g. "cd /home/T_Cha_2016-02-20") and type\n' +
                                                 '"irdap --makeconfig" to create a default configuration file "config.conf" in\n' +
                                                 'this directory. You can then adjust the parameters in the configuration file\n' +
                                                 'with a text editor. Finally in the terminal type "irdap --run" to perform the\n' +
                                                 'data reduction.\n\n' +
                                                 'The reduced images of two or more reductions can be mean-combined by typing\n' +
                                                 '"irdap --meancombine path1 path2 ... pathx", where the space-separated paths\n' +
                                                 'are absolute paths to the main directories of the reductions,\n' +
                                                 'e.g. "irdap --meancombine /home/T_Cha_2016-02-20 /home/T_Cha_2016-02-21".\n' +
                                                 'The mean-combined images will be written to the current working directory\n' +
                                                 'of the terminal.\n\n' +
                                                 'If this is the first time you use IRDAP, it is recommended to run the demo\n' +
                                                 'first by using the terminal to navigate to a directory of your choice and\n' +
                                                 'typing irdap --demo. Note that an internet connection is required as a small\n' +
                                                 'amount of raw data needs to be downloaded.\n\n' +
                                                 'When publishing data reduced with IRDAP, please cite van Holstein et al.\n' +
                                                 '(2020): https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..64V/abstract.\n' +
                                                 'For data in pupil-tracking mode please additionally cite van Holstein et al.\n' +
                                                 '(2017): https://ui.adsabs.harvard.edu/abs/2017SPIE10400E..15V.',
                                     epilog='Full documentation: https://irdap.readthedocs.io\n' +
                                            'Feedback, questions, comments: rob.vanholstein@eso.org\n\n' +
                                            'IRDAP Copyright (C) 2019 R.G. van Holstein',
                                            formatter_class=RawTextHelpFormatter)

    # Add parser arguments
    parser.add_argument('-v', '--version', action='store_true',
                        help='show program\'s version number')
    parser.add_argument('-w', '--website', action='store_true',
                        help='open IRDAP online documentation in web browser')
    parser.add_argument('-p', '--print', action='store_true',
                        help='toggle printing of log statements in the terminal')
    parser.add_argument('-d', '--demo', action='store_true',
                        help='run pipeline in current working directory with example\ndata of the circumstellar disk of T Cha (1 HWP cycle)')
    parser.add_argument('-o', '--headers', action='store_true',
                        help='create overview of relevant headers of FITS-files in raw\nsubdirectory')
    parser.add_argument('-c', '--makeconfig', action='store_true',
                        help='create default configuration file in current working\ndirectory')
    parser.add_argument('-r', '--run', action='store_true',
                        help='run pipeline using configuration file in current working\ndirectory')
    parser.add_argument('-m', '--meancombine', nargs='+', type=str, metavar='path',
                        help='mean-combine images of two or more reductions. The\n' \
                             'absolute paths to the main directories of the reductions\n' \
                             'should be supplied as arguments and be separated by\n' \
                             'spaces.')

    # Use current working directory (of terminal) as path of main directory of reduction
    path_main_dir = os.getcwd()

    # Evaluate and act upon user arguments
    args = parser.parse_args()

    if args.version:
        # Print the current version
        print('\nIRDAP version %s' % __version__)

    elif args.website:
        webbrowser.open_new_tab('https://irdap.readthedocs.io')

    if args.print:
        # Toggle printing in terminal
        path_file = os.path.join(os.path.dirname(__file__), 'print_in_terminal.txt')
        f = open(path_file, 'r')
        current_value = f.read()
        f.close()
        if current_value == 'True':
            print('\nIRDAP will not print log statements in the terminal.')
            f = open(path_file, 'w')
            f.write('False')
            f.close()
        elif current_value == 'False':
            print('\nIRDAP will print log statements in the terminal.')
            f = open(path_file, 'w')
            f.write('True')
            f.close()
        else:
            print('\nThe file ' + path_file + ' should contain either the word \'True\' or \'False\'.')

    elif args.demo:
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

    elif args.meancombine:
        # Mean-combine the images of two or more reductions
        path_read_dirs = args.meancombine

        if len(path_read_dirs) == 1:
            print('\nPlease provide at least two absolute paths to directories containing reduced data to be combined.')
        else:
            mean_combine_images(path_main_dir, path_read_dirs)

    # Check if latest version of IRDAP is used and if not suggest updating it
    url_github_version = 'https://raw.githubusercontent.com/robvanholstein/IRDAP/master/irdap/version.py'
    try:
        version_string = str(urllib.request.urlopen(url_github_version).readlines()[0], 'utf-8')
        version_github = version_string[version_string.rfind('=') + 1:].replace(' ', '').replace('\'', '')
    except:
        version_github = ''

    if version_github != '':
        if __version__ != version_github:
            print('\n\n\n\n\n\n\n\nA newer version of IRDAP is available (v' + __version__ + ' --> v' + version_github +
                  '). Please consider\n' +
                  'updating IRDAP by typing "pip install irdap --upgrade" in the terminal.')

###############################################################################
# Run the function main
###############################################################################

# Run function when called, i.e. in the terminal one can just write "irdap --run" i.o. "python -m irdap --run"
if __name__ == "__main__":
    main()