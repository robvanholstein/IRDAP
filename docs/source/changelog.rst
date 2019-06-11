
.. |last-commit| image:: https://img.shields.io/github/last-commit/robvanholstein/IRDAP.svg?colorB=e6c000
   :target: https://github.com/robvanholstein/IRDAP/
   
.. |issues| image:: https://img.shields.io/github/issues/robvanholstein/IRDAP.svg?color=b4001e
   :target: https://github.com/robvanholstein/IRDAP/issues

..
   |last-commit| |issues|

Changelog & to-do
=================

Changelog
---------

v0.2.1 June 2019, R.G. van Holstein:
 - Added the possibility to reduce data taken with the narrowband filters using the existing model of the broadband filters

v0.2.0 June 2019, R.G. van Holstein:
 - Changed the definition of the tuple of the input parameters *annulus_star*, *annulus_background* and *flux_annulus_background*: *width* is replaced by *outer_radius* and the definition of *start_angle* and *end_angle* is now the same as in DS9
 - Replaced static flat in Ks by one with P0-90 inserted
	
v0.1.2 June 2019, R.G. van Holstein:
 - Release of IRDAP package on PyPI and GitHub
 - Documentation put online
 
To-do high priority
-------------------

   - Combine data of different OB's with automatic selection of correct center, sky and flux images
   - Add ADI (including PCA) of total intensity images for pupil-tracking observations
   - Change config file option skip_preprocessing to perform_preprocessing: True/False, perform_polarimetry: True/False and perform_adi: True/False
   - Add option to express final images in contrast/arcsec^2 or mJansky/arcsec^2
   - Add plot of center coordinates of center frames vs time
   - Test polarimetric data with dithering applied (e.g. centering methods, subtraction of object file from center file)
   - Add accurate (calibrated) model correction for narrowband filters
   - Add option to apply accurate plate scale and distortion correction for data (especially important for pupil-tracking and bright sources; do we need calibrations?)	
   - Add detection limits (e.g. contrast curves)
   
To-do low priority
------------------
     
   - Add images and plots aimed at analysis of polarimetric data (R^2 scaling, AoLP arrows in polarized intensity or DoLP images, contrast curves etc.)	
   - Add option to subtract DARK,BACKGROUND-frames if sky files are lacking (especially important for Ks; test effect first before completely implementing)
   - Add option for 'stupid ADI' for field-tracking data with derotator offset and option to subtract 180 deg rotated image if no derotator offset
   - Add warning if (part of the) data is taken at a low polarimetric efficiency, and perhaps make plot of efficiency vs file number
   - Add optional ADI (including PCA) for polarization images
   - Make figures with sub-images horizontal, or make multiple lines of left and right images in a single figure   
   - Add options for various methods to shift and rotate images (interpolate, ndimage-fourier, sci-image functions; similar to VIP)
   - When rotating images, make complete frame still visible to not cut out any parts of the data
   - Improve centering of non-coronagraphic data (center found depends a lot on first PSF and affects Qphi and Uphi images)
   - Test finding of satellite spots of center files when waffle pattern is '+'
   - Add weighted least-squares as option for model correction (depending on image quality or polarimetric efficiency)	
   - Add multiprocessing for sigma filtering and centering of frames
   - Make IRDAP runable from a script by importing it and calling main function
   - Add option to scale the master sky frame to subtract from the object frames (especially for Ks; see also Gallicher et al. 2011)
   - Exclude saturated pixels in aperture to determine star polarization (same way as used in function fit_2d_gaussian)
   - Determine star polarization as a function of aperture radius	
   - Add optional RDI for total intensity images	
   - Make docstrings compliant with accepted conventions and create API doc on website