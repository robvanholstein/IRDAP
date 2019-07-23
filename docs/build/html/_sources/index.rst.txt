.. IRDAP documentation master file, created by
   sphinx-quickstart on Sat May 11 17:17:28 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _keppler: https://ui.adsabs.harvard.edu/abs/2018A%26A...617A..44K

.. |keppler| replace:: *(Keppler et al. 2018),*

.. _ginski: https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..79G

.. |ginski| replace:: *(Ginski et al. 2018)*

.. _pohl: https://ui.adsabs.harvard.edu/abs/2017A%26A...605A..34P

.. |pohl| replace:: *(Pohl et al. 2017;*

.. _vanholstein: ADS link

.. |vanholstein| replace:: *van Holstein et al. 2019).*

.. |stars| image:: https://img.shields.io/github/stars/robvanholstein/IRDAP.svg?style=social&label=Stars
   :target: https://github.com/robvanholstein/IRDAP/
   
.. |watch| image:: https://img.shields.io/github/watchers/robvanholstein/IRDAP.svg?style=social&label=Watch
   :target: https://github.com/robvanholstein/IRDAP/
   
.. |pypi| image:: https://img.shields.io/pypi/v/irdap.svg?colorB=<brightgreen>
    :target: https://pypi.python.org/pypi/irdap/
	
.. |python| image:: https://img.shields.io/badge/Python-3.6%2C%203.7-yellow.svg?style=flat
    :target: https://pypi.python.org/pypi/irdap/

.. |github| image:: https://img.shields.io/github/release/robvanholstein/IRDAP.svg
   :target: https://github.com/robvanholstein/IRDAP/ 
   
.. |last-commit| image:: https://img.shields.io/github/last-commit/robvanholstein/IRDAP.svg?colorB=e6c000
   :target: https://github.com/robvanholstein/IRDAP/

.. |license| image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: https://github.com/robvanholstein/IRDAP/blob/master/LICENSE

.. |ads1| image:: https://img.shields.io/badge/ADS-%3CADS%20link%3E-blueviolet.svg
	:target: <ADS link> ALSO CHANGE IN README!!!

.. |ads2| image:: https://img.shields.io/badge/ADS-van%20Holstein%20et%20al.%20(2017)-blueviolet.svg
	:target: https://ui.adsabs.harvard.edu/abs/2017SPIE10400E..15V
	
.. Made ads-link above on https://shields.io/ with "your badge"
	
.. IRDAP |stars| |watch|

IRDAP
=================================

.. 
   |pypi| |python| |github| |last-commit| |license| |ads1| |ads2|

IRDAP (IRDIS Data reduction for Accurate Polarimetry) is a highly-automated end-to-end pipeline to reduce `SPHERE-IRDIS <https://www.eso.org/sci/facilities/paranal/instruments/sphere.html>`_ polarimetric data. Its core feature is the model-based correction method of the instrumental polarization effects as described in `van Holstein et al. (2019) <ADS link>`_. IRDAP handles data taken both in field- and pupil-tracking mode and using the broadband filters Y, J, H and K\ :sub:`s`. Data taken with the narrowband filters can be reduced as well, although with a somewhat worse accuracy. 

Reducing data with IRDAP is very straightforward and does not require the user to do any coding or have knowledge of Python (IRDAP is written for Python 3.6 and 3.7). IRDAP is simply run from a terminal with only a few commands and uses a configuration file with a limited number of input parameters. Within several minutes, IRDAP performs a complete data reduction from raw data to final data products.

.. note::
   If you use IRDAP for your publication, please :ref:`cite our papers <Citing IRDAP>`.

.. figure:: ./figs/home_3disks.png
    :width: 750px

    *Three examples of data sets reduced with IRDAP: the circumstellar disk of PDS 70 including its inner disk that polarizes the central star* |keppler|_ *the circumstellar disk and polarized companion of CS Cha* |ginski|_ *and the circumstellar disk of T Cha* |pohl|_ |vanholstein|_

IRDAP yields a multitude of improvements for observations of circumstellar disks: it enables us to accurately study the morphology of disks, measure non-azimuthal polarization and determine scattering phase functions and particle properties. Because IRDAP discerns instrumental polarization from stellar polarization, it is a vital tool for accurate radiative transfer modeling of disks and enables the detection of unresolved (inner) disks and the measurement of the polarization of substellar companions. Finally, IRDAP enables accurate data reduction for targets that cannot be reduced with the conventional data-reduction methods (because there is no bright star in the field of view), such as solar system objects and galaxies.

Contents
---------
.. toctree::
   :maxdepth: 1
   
   Home <self>
   installation
   example
   instructions
   configfile
   citing
   contributing    
   changelog
   mailinglist
   acknowledgements
   GitHub <https://github.com/robvanholstein/IRDAP>

----

Copyright notice:
 
The IRDAP Python module is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software Foundation, 
version 3 of the License.
 
The IRDAP Python module is distributed in the hope that it will be useful, but without 
any warranty; without even the implied warranty of merchantability or fitness for a 
particular purpose. See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with the IRDAP 
Python module. If not, see http://www.gnu.org/licenses/.