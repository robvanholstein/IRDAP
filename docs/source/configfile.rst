
Configuration file
==================

The configuration file with the input parameters for IRDAP is shown below:

.. include:: ../../irdap/config.conf
   :literal:

The input parameters are divided into five groups:

- :ref:`Basic pre-processing options`
- :ref:`Basic PDI options`
- :ref:`Basic ADI options`
- :ref:`Advanced pre-processing options`
- :ref:`Advanced PDI options`

Below we discuss for each parameter of the configuration file what it does and what the valid input options are.

.. important::
   While the user will quite frequently have to adapt the basic options, the advanced options generally do not need to be changed. 
	
Basic pre-processing options
----------------------------

.. _perform_preprocessing:

.. py:function:: perform_preprocessing:

   ``True``, ``False`` (default = ``True``)
   
   If True, perform :ref:`pre-processing <Pre-processing in a nutshell>` on the raw data, i.e. prepare the data for the subsequent :ref:`polarimetric differential imaging (PDI) <PDI in a nutshell>` and/or :ref:`angular differential imaging (ADI) <ADI in a nutshell>` steps. `perform_preprocessing` must be ``True`` the first time IRDAP is run on a particular data set. When re-running IRDAP on the same data set, `perform_preprocessing` can be set to ``False`` to skip the pre-processing and save time when tweaking the input parameters of the PDI and/or ADI steps.


.. py:function:: sigma_filtering:

   ``True``, ``False`` (default = ``True``)
   
   If ``True``, use sigma-filtering to remove bad pixels that are still in the OBJECT-, CENTER- and FLUX-frames after applying the master bad pixel map. After the sigma-filtering the frames should have no bad pixels. If ``False``, only use the master bad pixel map to speed up the reduction. The latter is not recommended for a final reduction however, because a signficant number of bad pixels will remain in the frames.


.. py:function:: object_collapse_ndit:

   ``True``, ``False`` (default = ``False``)

   If ``True``, compute the mean over the (NDIT) frames of the OBJECT-files before sky (or background) subtraction, flat-fielding, bad pixel removal and centering to speed up the pre-processing. If ``False``, perform the above steps for each frame and then compute the mean over the frames. The latter will result in the most accurate reduction and is recommended for a final reduction.


.. _object_centering_method:

.. py:function:: object_centering_method:
   
   ``automatic``, ``center frames``, ``gaussian``, ``cross-correlation``, ``manual`` (default = ``automatic``)

   Method to center the OBJECT-frames. If ``center frames``, center the OBJECT-frames with the center coordinates determined from the CENTER-frames. If ``gaussian``, center the OBJECT-frames by fitting a 2D Gaussian to each frame. If ``cross-correlation``, perform the centering by fitting a 2D Gaussian to the first frame and then using cross-correlation to align (register) the other frames onto the centered first frame. Generally ``gaussian`` yields more accurate results than ``cross-correlation``. For ``gaussian`` and ``cross-correlation`` the determined center coordinates are plotted for each image. If ``manual``, use fixed coordinates as defined with object_center_coordinates_ to center the OBJECT-frames. If ``automatic``, center the OBJECT-frames using the CENTER-files if they exist, and fit a 2D Gaussian to each OBJECT-frame if there are no CENTER-files.

   .. attention::
      In case the centering is not accurate, the advanced parameters center_subtract_object_, center_param_centering_, object_center_coordinates_ and object_param_centering_ can be adapted. In addition, if the shifts among the frames are small when centering using ``gaussian`` or ``cross-correlation``, the quality of the final images can sometimes be improved by repeating the pre-processing with object_centering_method_ set to ``manual`` and object_center_coordinates_ equal to the mean center coordinates found earlier with the Gaussian fitting or cross-correlation (see the log file).

.. py:function:: frames_to_remove:

   `list` (default = ``[]``)

   List of integers and length-2-tuples of integers indicating which files and frames not to use for the data reduction. This parameter allows the user to remove frames of bad quality and to exclude complete files from the reduction without (re)moving them from the raw subdirectory. A complete file can be removed by specifying its integer index, while a frame of a specific file can be removed by specifying a tuple ``(file_index, frame_index)``. For example ``[2, 4, (5, 3)]`` will exclude files 2 and 4 from the reduction and will remove frame 3 of file 5. If no files or frames should be removed, use an empty list ``[]``. The files are sorted in chronological order from oldest to newest. The file indices of each file are shown in the header overview that is created when you run IRDAP.
			
   .. attention::
      The indices of the files and frames are 1-based, i.e. the first file or frame has index 1.

   .. hint::
      To check for frames of bad quality or frames that are not centered correctly, the user can examine the sub-images of the (pre-)processed OBJECT- and FLUX-frames that IRDAP writes to the ``preprocessed`` and ``calibration/flux`` subdirectories. Because the file and frame indices are displayed next to the sub-images, the bad frames can straightforwardly be removed by inserting the corresponding indices in the list of `frames_to_remove`.


Basic PDI options
-----------------------------

.. py:function:: perform_pdi:

   ``True``, ``False`` (default = ``True``)
   
   If True, perform :ref:`polarimetric differential imaging (PDI) <PDI in a nutshell>` using the pre-processed data. To perform this step, perform_preprocessing_ must be ``True`` or the raw data must have been pre-processed before. If ``False``, do not perform PDI. The latter can be useful when tweaking the input parameters of the ADI step only.


.. _annulus_star:

.. py:function:: annulus_star:

   ``automatic``, ``ao residuals``, ``star aperture``, `list`, `tuple` (default = ``automatic``)

   Annulus used to measure the polarization of the central star (or any other source in the field of view) from the *Q*-, *I*\ :sub:`Q`-, *U*-, and *I*\ :sub:`U`-images. This measured polarization signal is subtracted from the polarization images written to the subdirectory ``reduced_star_pol_subtr``. For the most accurate results the annulus should only contain signal from the star, and no signal from for example a circumstellar disk, companion or background star. The measured polarization signal is affected by the subtraction of the background in the images (see annulus_background_).
   
   If ``ao residuals``, the annulus will be centered on the central star and be located over the AO residuals. The inner and outer radius of the annulus will automatically be adapted to the filter used. If ``star aperture``, the star polarization will be determined from a small aperture with a radius of 11 pixels located at the position of the central star. If ``automatic``, an annulus over the AO residuals (as for ``ao residuals``) will be used in case of coronagraphic data, and a small aperture at the position of the central star (as for ``star aperture``) in case of non-coronagraphic data. 
   
   Generally the user would do a first data reduction using ``automatic``. If after this first reduction it appears that more control over the exact shape of the annulus is required, the user can define the annulus with a tuple of 6 parameters:
   
   - *coord_center_x*: x-coordinate of center (pixels)
   - *coord_center_y*: y-coordinate of center (pixels)
   - *inner_radius*: Inner radius (pixels)
   - *outer_radius*: Outer radius (pixels)
   - *start_angle*: Start angle of annulus sector (deg)
   - *end_angle*: End angle of annulus sector (deg)
   
   For example ``(512.5, 512.5, 60, 95, 0, 360)`` will create an annulus that is centered on the central star and is located over the AO residuals in H-band.

   .. attention::
      IRDAP uses the same 1-based indexed coordinates as `DS9 <http://ds9.si.edu/>`_. Therefore coordinates read from DS9 can be used for the configuration file without applying any conversion.
	   
      IRDAP centers the frames at the coordinates (x, y) = (512.5, 512.5). Therefore a centered annulus will have *coord_center_x* and *coord_center_y* equal to 512.5.

   .. attention::
      *start_angle* and *end_angle* are defined with respect to the orientation of the final images with 0 deg oriented to the right and positive angles rotating counterclockwise. Therefore an angle of 0 deg points to the West, except when single_posang_north_up_ is set to ``False`` for observations taken in field-tracking mode with a single derotator position angle.
	   
      The definition of the angles is the same as used in `DS9 <http://ds9.si.edu/>`_. When using DS9 to draw a line starting at the center coordinates, the angle indicated can be used for the configuration file without applying any conversion.

   If even more control over the shape of the annulus is required, the user can define a list of length-6-tuples. This is for instance useful to exclude signal of a circumstellar disk from the annulus. ``[(512.5, 512.5, 60, 95, 105, 255), (512.5, 512.5, 60, 95, 285, 75)]`` will for example create an annulus with a gap at the top and bottom each spanning 30 deg.


.. _annulus_background:

.. py:function:: annulus_background:

   ``large annulus``, `list`, `tuple` (default = ``large annulus``)

   Annulus used to measure and remove the background in the *Q*-, *I*\ :sub:`Q`-, *U*-, and *I*\ :sub:`U`-images. This background subtraction also affects the computation of the polarization of the central star (see annulus_star_). For the most accurate results the annulus should not contain any signal from the star or any other source in the field of view. 
  
   If ``large annulus``, the annulus will be centered on the central star and be located far away from the star with an inner radius of 360 pixels and an outer radius of 420 pixels. In case more control over the exact shape of the annulus is required, the user can define the annulus with a (list of) length-6-tuple(s) in the same way as for annulus_star_.


.. py:function:: normalized_polarization_images:

   ``True``, ``False`` (default = ``False``)

   If ``True``, create final images of degree of linear polarization, normalized Stokes *q* and *u* (see `van Holstein et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..64V/abstract>`_) and degree and angle of linear polarization computed from *q* and *u*: 
   
   - *DoLP* = sqrt(*Q*\ :sup:`2`\ + *U*\ :sup:`2`\) / [0.5 (*I*\ :sub:`Q` + *I*\ :sub:`U`)]
   - *q* = *Q* / *I*\ :sub:`Q`
   - *u* = *U* / *I*\ :sub:`U`
   - *AoLP_norm* = 0.5 arctan(*u* / *q*)
   - *DoLP_norm* = sqrt(*q*\ :sup:`2`\ + *u*\ :sup:`2`\)
   
   These additional final images are only valid if all flux in the images originates from the source of interest. This is generally the case for observations of for example solar system objects or galaxies. The images are generally not valid for observations of circumstellar disks or companions because in that case a large part of the flux in the total intensity images originates from the central star. 
   
   *DoLP_norm* and *AoLP_norm* are potentially more accurate than *DoLP* (above) and *AoLP* = 0.5 arctan(*U* / *Q*) (as always created), especially when there are significant variations in seeing and sky transparency among the measurements. If ``False``, do not create these images.

Basic ADI options
-----------------

.. _perform_adi:

.. py:function:: perform_adi:

   ``True``, ``False`` (default = ``True``)
   
   If True, perform :ref:`angular differential imaging (ADI) <ADI in a nutshell>` using the pre-processed data. ADI is only performed on data taken in pupil-tracking mode. To perform this step, perform_preprocessing_ must be ``True`` or the raw data must have been pre-processed before. If ``False``, do not perform ADI. The latter can be useful when tweaking the input parameters of the PDI step only.

.. _principal_components:

.. py:function:: principal_components:

   ``companion+disk``, ``companion``, ``disk``, `integer`, `list` (default = ``companion+disk``)
  
   The number of principal components to subtract for the reduction combining angular differential imaging (ADI) and principal component analysis (PCA). If ``companion``, create two frames with 10 and 16 principal components subtracted to search for companions. If ``disk``, create two frames with 2 and 4 principal components subtracted to search for disk signal. If ``companion+disk``, create four frames with 2, 4, 10 and 16 principal components subtracted to search for companions and disk signal. 
   
   If more control is required, a strictly positive integer or a list of unique and strictly positive integers can be provided. For example providing ``[10, 16]`` will yield two frames with 10 and 16 principal components subtracted. Because we cannot remove more components than the number of frames in the cube of pre-processed frames, the number of components to be subtracted is reduced if there are not enough frames.
   
.. _pca_radii:

.. py:function:: pca_radii:

   ``automatic``, `list` (default = ``automatic``)

	Annuli used to optimize the ADI+PCA reduction over. If ``automatic``, the reduction is optimized over three annuli with inner and outer radii equal to 10 and 30, 30 and 100, and 100 and 512 pixels. If more control is required, the user can provide a list of at least length 2 containing positive and increasing integers (including 0) indicating the inner and outer radii of the annuli in pixels. For example, if the list has two elements, these elements define the inner and outer radius of a single annulus. If the list has more elements, the second element defines the outer radius of the first annulus and the inner radius of the second annulus. Likewise the third element defines the outer radius of the second annulus and the inner radius of the third annulus, etc. Therefore setting `pca_radii` to ``[10, 30, 100, 512]`` is equivalent to setting it to ``automatic``.

Advanced pre-processing options
-------------------------------

.. _center_subtract_object:

.. py:function:: center_subtract_object:

   ``True``, ``False`` (default = ``True``)
   
   If ``True``, subtract the OBJECT-file(s) taken closest in time from the CENTER-file(s) before determining the center coordinates from the satellite spots. This results in a more accurate determination of the center coordinates because the background and any other sources in the field of view are suppressed. When the difference between the image orientations of the CENTER- and OBJECT-frames is large (i.e. difference in derotator position angle for field-tracking and difference in parallactic angle for pupil-tracking), the OBJECT-frame will be rotated around the initial guess of the centers as defined by object_center_coordinates_ before subtracting it from the CENTER-frame. If ``False``, do not subtract the OBJECT-file(s). Parameter is ignored if no CENTER-files are used for the centering (see object_centering_method_).

.. _center_param_centering:

.. py:function:: center_param_centering: 

   `tuple` (default = ``(12, None, 30000)``)

   Tuple of 3 parameters pertaining to the fitting of the satellite spots of the the CENTER-frames with 2D Gaussians to determine the center coordinates. In case the center coordinates are not accurately determined, the 3 parameters can be adapted in addition to object_center_coordinates_: 
   
   - *crop_radius*: Half the length of the sides of the square cropped sub-images that should each contain a satellite spot and are used to fit the 2D Gaussians to (pixels). Must be integer. object_center_coordinates_ gives the initial guess of the center coordinates to place the sub-images approximately at the locations of the satellite spots. If ``None``, the complete frame is used for the fitting and object_center_coordinates_ is ignored.
   - *sigfactor*: All sub-image pixels with values smaller than sigfactor*standard deviation are replaced by random Gaussian noise to mask them for fitting the 2D Gaussian. This can be useful when there is a lot of noise in the background. If ``None``, no pixels are replaced by Gaussian noise.
   - *saturation_level*: All pixels within the smallest circle encompassing the pixels with a value equal to or higher than saturation_level are ignored when fitting the 2D Gaussian. Ignoring the saturated pixels improves the accuracy of the fit. We use a circle because strongly saturated pixels in the peak of the PSF often have values lower than saturation_level. If ``None``, no pixels are ignored.
   
   `center_param_centering` is ignored if no CENTER-files are used for the centering (see object_centering_method_).


.. _object_center_coordinates:

.. py:function:: object_center_coordinates:

   ``automatic``, `tuple` (default = ``automatic``)

   Tuple of 4 floats with (initial guess of) center coordinates of OBJECT-frames:

   - *x_left*: x-coordinate of center of left frame half (pixels)
   - *y_left*: y-coordinate of center of left frame half (pixels)
   - *x_right*: x-coordinate of center of right frame half (pixels)
   - *y_right*: y-coordinate of center of right frame half (pixels)

   If ``automatic``, `object_center_coordinates` is set to ``(478, 535, 1504, 524)`` when the coronagraph `N_ALC_Ks` is used and to ``(478, 522, 1504, 512)`` otherwise.
     
   .. attention::
      IRDAP uses the same 1-based indexed coordinates as `DS9 <http://ds9.si.edu/>`_. Therefore coordinates read from DS9 can be used for the configuration file without applying any conversion.

      The center coordinates are defined in the complete frame, i.e. with both detector halves. Therefore *x_right* has a value larger than 1024 pixels.

   The center coordinates of `object_center_coordinates` are used at various steps of the pre-processing:

   - During the processing of the CENTER-files, they give the coordinates to rotate the OBJECT-frame around before subtracting it from the CENTER-frame (see center_subtract_object_). 
   - When determining the center coordinates from the CENTER-frames, they give the initial guess of the center coordinates to place the sub-images approximately on the location of the satellite spots (see center_param_centering_). 
   - They give the initial guess of the center coordinates when centering the OBJECT-frames by fitting 2D Gaussians or using cross-correlation (see object_centering_method_ and object_param_centering_). 
   - In case object_centering_method_ is ``manual``, they are the actual center coordinates used to center the OBJECT-frames.


.. _object_param_centering:

.. py:function:: object_param_centering:

   `tuple` (default = ``(60, None, 30000)``)

   Tuple of 3 parameters pertaining to the centering of the OBJECT-frames by fitting a 2D Gaussian to each frame or using cross-correlation. In case the centering is not performed accurately, the 3 parameters can be adapted in addition to object_center_coordinates_ (especially useful for non-coronagraphic data of multiple stars): 

   - *crop_radius*: Half the length of the sides of the square cropped sub-images used to fit the 2D Gaussians to and used for cross-correlating the images (pixels). Must be integer. The sub-image is centered on the coordinates as provided by object_center_coordinates_. If ``None``, the complete frame is used and object_center_coordinates_ is ignored.
   - *sigfactor*: see center_param_centering_.
   - *saturation_level*: see center_param_centering_.
   
   `object_param_centering` is only used when object_centering_method_ is ``gaussian`` or ``cross-correlation``.


.. _flux_centering_method:

.. py:function:: flux_centering_method:

   ``gaussian``, ``manual`` (default = ``gaussian``)
   
   Method to center the FLUX-frames. If ``gaussian``, fit a 2D Gaussian to each frame. The determined center coordinates are plotted for each image. If ``manual``, use fixed coordinates as defined with flux_center_coordinates_. In case the centering is not accurate, the parameters flux_center_coordinates_ and flux_param_centering_ can be adapted.


.. _flux_center_coordinates:

.. py:function:: flux_center_coordinates:

   `tuple` (default = ``(478, 522, 1504, 512)``)

   (Initial guess of) center coordinates of FLUX-frames, defined in the same way as object_center_coordinates_. The center coordinates are the initial guess when centering the FLUX-frames by fitting 2D Gaussians (see flux_centering_method_ and flux_param_centering_). In addition, in case flux_centering_method_ is ``manual``, the center coordinates are the actual center coordinates used to center the FLUX-frames.


.. _flux_param_centering:

.. py:function:: flux_param_centering:

   `tuple` (default = ``(60, None, 30000)``)

   Tuple of 3 parameters pertaining to the centering of the FLUX-frames by fitting a 2D Gaussian to each frame. In case the centering is not performed accurately, the 3 parameters can be adapted in addition to flux_center_coordinates_ (especially useful in case of multiple stars):
   
   - *crop_radius*: Half the length of the sides of the square cropped sub-images used to fit the 2D Gaussians to (pixels). Must be integer. The sub-image is centered on the coordinates as provided by flux_center_coordinates_. If ``None``, the complete frame is used and flux_center_coordinates_ is ignored.
   - *sigfactor*: see center_param_centering_.
   - *saturation_level*: see center_param_centering_.
   
   `flux_param_centering` is only used when flux_centering_method_ is ``gaussian``.  


.. _flux_annulus_background:

.. py:function:: flux_annulus_background:

   ``large annulus``, `list`, `tuple` (default = ``large annulus``)

   Annulus used to remove the background in the FLUX-frames. For the most accurate results the annulus should not contain any signal from the star or any other source in the field of view. If ``large annulus``, the annulus will be centered on the central star and be located far away from the star with an inner radius of 320 pixels and an outer radius of 380 pixels. In case more control over the exact shape of the annulus is required, the user can define the annulus with a (list of) length-6-tuple(s) in the same way as for annulus_star_.


.. _flux_annulus_star:

.. py:function:: flux_annulus_star:

   ``automatic``, `list`, `tuple` (default = ``automatic``)

   Annulus used to measure the total flux of the star in the master flux frames. The total fluxes are converted into reference fluxes that can be used to express the final images (e.g. the *I*\ :sub:`Q`- or *Q*:math:`_\phi`-images) in units of contrast/arcsec\ :sup:`2` or Jansky/arcsec\ :sup:`2`. The annulus should contain all signal from the central star, but no signal from for example a (stellar) companion or background star. 
   
   If ``automatic``, the annulus will be an aparture of radius 120 pixels located at the position of the central star. In case more control over the exact shape of the annulus is required, for example to exclude an object close to the central star, the user can define the annulus with a (list of) length-6-tuple(s) in the same way as for annulus_star_.

Advanced PDI options
--------------------------------

.. py:function:: double_difference_type:

   ``conventional``, ``normalized`` (default = ``conventional``)

   Type of double difference to be computed. In almost all cases one would use ``conventional``. When there are large variations in atmospheric seeing and/or sky transparency among the measurements (often resulting in a large spread of the measured star polarization as a function of HWP cycle number), using ``normalized`` can suppress spurious polarization signals and improve the quality of the final images (see `van Holstein et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..64V/abstract>`_). 

.. _single_posang_north_up:

.. py:function:: single_posang_north_up:
   
   ``True``, ``False`` (default = ``True``)

   For observations taken in field-tracking mode with a single derotator position angle (header keyword ``INS4.DROT2.POSANG`` either 0 or with a fixed offset), the final images are rotated with North up if ``True``, and kept in the orientation of the raw frames if ``False``. In the latter case, the images are more accurate as they suffer less from interpolation errors caused by image rotation. This is useful when for example extracting the polarized surface brightness distribution of a circumstellar disk. Parameter is ignored for field-tracking observations with more than one derotator position angle or observations taken in pupil-tracking mode. In these cases the final images produced are always oriented with North up.









..
   Variables not used anymore

..
   py:function:: combination_method_polarization (``least squares`` or ``trimmed mean`` or ``median``): Method to be used to produce the incident Q- and U-images, i.e. the images that are corrected for the instrumental polarization effects. Valid values are ``least squares``, ``trimmed mean`` or ``median``. The recommended option is ``trimmed mean``. With ``least squares`` the images are obtained by solving for every pixel the system of equations describing the measurements using linear least squares (see Eq. 35 of `van Holstein et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..64V/abstract>`_). With ``trimmed mean`` or ``median`` the images are obtained by solving the system of equations for each pair of double-difference Q- and U-images (each HWP cycle) separately, and then computing the trimmed mean or median over all resulting images. ``least squares`` is the most accurate option, but any unremoved bad pixels will still be visible in the images. Using ``median`` will remove these bad pixels, but is the least accurate option and also yields images with a lower signal-to-noise ratio as is clear from images of circumstellar disks. Using ``trimmed mean`` will yield images that have essentially the same accuracy and signal-to-noise ratio images produced using ``least squares``, but without the bad pixels. Therefore ``trimmed mean`` is the recommended option.

..
   py:function:: combination_method_intensity (``mean`` or ``trimmed mean`` or ``median``): Method to be used to produce the incident I_Q- and I_U-images. These images are computed by combining the I_Q- or I_U-images of all HWP cycles using the ``mean``, ``trimmed mean`` or ``median``. ``mean`` yields the most accurate images, but any unremoved bad pixels will still be visible in the images. Using ``median`` will remove these bad pixels, but is the least accurate option. ``trimmed mean`` produces images similar to ``mean``, but without the bad pixels. When using ``trimmed mean``, also check the input parameter **trimmed_mean_prop_to_cut_intens**. It is generally recommended to use either ``trimmed mean`` or ``mean``.
   
.. 
   py:function:: trimmed_mean_prop_to_cut_polar (`float`): Fraction to cut off of both tails of the distribution if ``trimmed mean`` is used for **combination_method_polarization**. Parameter is ignored in case ``least squares`` or ``median`` is used. Value should be in range 0 <= *trimmed_mean_prop_to_cut_polar* <= 1. In most cases a value of 0.1 or 0.15 removes the bad pixels well while producing images very similar to those obtained with ``least squares``.

.. 
   py:function:: trimmed_mean_prop_to_cut_intens (`float`): Fraction to cut off of both tails of the distribution if ``trimmed mean`` is used for **combination_method_intensity**. Parameter is ignored in case ``mean`` or ``median`` is used. Value should be in range 0 <= *trimmed_mean_prop_to_cut_intens* <= 1. In most cases a value of 0.1 or 0.15 removes the bad pixels well while producing images similar to those obtained with ``mean``.
