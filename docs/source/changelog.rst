
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

v1.2.3 December 2019, R.G. van Holstein
 - Clarified meaning of standard deviation of fitted center coordinates in log file
 - Very minor correction to computation of HWP angle from headers
 - Corrected minor bug when computing polarimetric efficiency when the number of Q and U cycles is unequal
 - Added warning when spread of star polarization per HWP cycle is large and suggest to try normalized double difference
 
v1.2.2 November 2019, R.G. van Holstein
 - Added polarimetric efficiency to crosstalk/transmission plot, a print statement with the range of the polarimetric efficiency of the observations, and a warning in case of low polarimetric efficiency
 - Corrected a bug in the checks of flux_annulus_star
 - Corrected a bug that when mean combining ADI+PCA-images only the first frame of the cube was written to the FITS-files

v1.2.1 November 2019, R.G. van Holstein
 - Corrected bug introduced in v1.2.0 with computation of DoLP and AoLP of star when there is an unequal number of Q and U measurements
 
v1.2.0 November 2019, R.G. van Holstein
 - Added computation of statistical uncertainty (photon and background noise) on measured star polarization in PDI part
 - Added contrast curves for the detection of polarized extended sources and polarized point source with PDI
 
v1.1.0 October 2019, R.G. van Holstein & J. Milli
 - Added angular differential imaging (ADI) including principal component analysis (PCA) of total intensity images for pupil-tracking observations
 - Added mean combination of ADI final images to command-line option :code:`irdap --meancombine`

v1.0.0 September 2019, R.G. van Holstein & J. Milli:
 - No backward compatibility with previous versions
 - Added analysis of master flux frames to enable the user to express the final images in contrast/arcsec^2 or Jansky/arcsec^2
 - Restructered code and configuration file to allow for the separate execution of pre-processing, polarimetric differential imaging (previously referred to as post-processing) and angular differential imaging 
 - Added plot of center coordinates of center frames vs time
 - Added possibility to open documentation through terminal command :code:`irdap --website`
 - Confirmed that polarimetric data with dithering applied is correctly reduced
 
v0.3.0 July 2019, R.G. van Holstein:
 - Added command-line option :code:`irdap --meancombine` to mean-combine final images of multiple observation blocks

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

   - Make the data reduction work for data taken before approximately May 2015 that do not contain the header keyword indicating the Stokes parameter (using the HWP angle instead)
   - Make IRDAP throw out less files when the order of the Stokes parameters of the files is not Qplus, Qminus, Uplus, Uminus
   - Add possibility to use CENTER-frames with a different exposure time than the OBJECT-frames
   - Add possibility to crop images when doing ADI to save time
   - Create contrast curves of final ADI-images 
   - Add ADI+PCA for polarimetric data reduction of pupil-tracking observations
   - Add accurate (calibrated) model correction for narrowband filters
   - Add option to apply accurate plate scale and distortion correction for data (especially important for pupil-tracking and bright sources; do we need calibrations?)
   - For pupil-tracking data mask bad pixel clusters when computing least squares solution to get rid off sectors of bad data; for field-tracking data set bad pixel clusters to NaN   

   
To-do low priority
------------------

   - Add images and plots aimed at analysis of polarimetric data (R^2 scaling, AoLP arrows in polarized intensity or DoLP images, contrast curves etc.)	
   - Add option to subtract DARK,BACKGROUND-frames if sky files are lacking (especially important for Ks; test effect first before completely implementing)
   - Add option for 'stupid ADI' for field-tracking data with derotator offset and option to subtract 180 deg rotated image if no derotator offset
   - Make figures with sub-images horizontal, or make multiple lines of left and right images in a single figure   
   - Add options for various methods to shift and rotate images (interpolate, ndimage-fourier, sci-image functions; similar to VIP)
   - When rotating images, make complete frame still visible to not cut out any parts of the data
   - Improve centering of non-coronagraphic data (center found depends a lot on first PSF and affects Qphi and Uphi images). Perhaps fit coordinates on each PSF, but do the actual shifts with the mean of these fitted values. This has proven to give a more accurate final result. 
   - Test finding of satellite spots of center files when waffle pattern is '+'
   - Add weighted least-squares as option for model correction (depending on image quality or polarimetric efficiency)	
   - Add multiprocessing for sigma filtering and centering of frames
   - Add option to scale the master sky frame to subtract from the object frames (especially for Ks; see also Gallicher et al. 2011)
   - Exclude saturated pixels in aperture to determine star polarization (same way as used in function fit_2d_gaussian)
   - Determine star polarization as a function of aperture radius	
   - Add optional RDI for total intensity images	
   - Make docstrings compliant with accepted conventions and create API doc on website