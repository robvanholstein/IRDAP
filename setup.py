'''
This file contains the parameters to setup the IRDAP package.

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

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
   name='irdap',
   version='2019.5.11',
   author='R.G. van Holstein',
   author_email='vanholstein@strw.leidenuniv.nl',
   packages=['irdap'],
   url='https://github.com/robvanholstein/IRDAP',
   download_url='https://github.com/robvanholstein/IRDAP/archive/master.zip',
   license='GNU General Public License v3.0',
   description='Pipeline for accurate reduction of SPHERE-IRDIS polarimetric data.',
   long_description=long_description,
   python_requires='>=3',
   install_requires=[
      "numpy >= 1.16.1",
      "matplotlib >= 3.0.1",
      "scipy >= 1.2.1",
      "astropy >= 3.1.1",
      "scikit_image >= 0.14.2"
   ],
   classifiers=[
   # How mature is this project? Common values are
   #   3 - Alpha
   #   4 - Beta
   #   5 - Production/Stable
   'Development Status :: 3 - Alpha',
   'Intended Audience :: Science/Research',
   'Topic :: Scientific/Engineering :: Astronomy',
   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
   'Programming Language :: Python :: 3.6',
   ],
   include_package_data=True, # So that non .py files make it onto pypi, and then back !
   entry_points={'console_scripts': ['irdap = irdap.__main__:main']},
)