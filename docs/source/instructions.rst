
Usage instructions
==================

Running IRDAP
-------------

Reducing data with IRDAP is very straightforward and only requires two commands to be entered in the terminal. First, create a directory (e.g. :file:`/home/T_Cha_2016-02-20`) containing a subdirectory called :file:`raw` (i.e. :file:`/home/T_Cha_2016-02-20/raw`) in which you place the raw FITS-files. 

.. attention::
   The raw FITS-files should have `DPR TYPE` equal to OBJECT, OBJECT,CENTER, OBJECT,FLUX or SKY. If the user wishes to make a new master flat and bad pixel map (see :ref:`Pre-processing in a nutshell`) `DPR TYPE` can also be equal to FLAT,LAMP, DARK,BACKGROUND or DARK. FITS-files with a `DPR TYPE` other than the ones mentioned above will automatically be ignored.

.. hint::
   To only download the raw data that is required for IRDAP, use the `SPHERE data archive <http://archive.eso.org/wdb/wdb/eso/sphere/form>`_ and set the search field `DPR TECH` to `any POLARIMETRY`.
   
In the terminal navigate to the main directory (e.g. ``cd /home/T_Cha_2016-02-20``) and type:
::
 
   irdap --makeconfig

This creates a default configuration file :file:`config.conf` in the main directory. 
You can then adjust the parameters in the configuration file with a text 
editor (see :ref:`Configuration file`). Finally, in the terminal type:
::

   irdap --run

to perform the data reduction. 

.. note::
	In general a first reduction can be performed without making changes to the configuration file. After this first reduction, some of the input parameters can be adjusted and IRDAP can be run again (but this is often not necessary!).
 
Pre-processing in a nutshell
----------------------------

The data reduction is divided into a pre-processing, polarimetric differential imaging (PDI) and angular differential imaging (ADI) part. 

.. warning::
   The data-reduction steps are only described concisely here. The log file that IRDAP creates describes the steps in more detail and also briefly describes the output files created. It is therefore recommended to carefully read the log file the first time you run IRDAP.
 
The pre-processing starts with reading the configuration file and creating the log file that contains all the lines printed on screen. In addition, IRDAP creates overviews of the relevant header keywords for each file in the :file:`raw` subdirectory, sorts the raw data and creates subdirectories to write the output files to.

.. admonition:: Hints

   - The header overviews can also be created without running IRDAP by typing :code:`irdap --headers` in the terminal;
   - When reducing data with IRDAP, keep an eye out for the WARNING messages printed on screen and in the log file; 
   - When running IRDAP a copy of the configuration file is made. When subsequently re-running IRDAP, this copy is moved together with the log file to the ``logs`` subdirectory as a record of the performed reductions. 
   - In case some steps of the reduction are skipped when re-running IRDAP (i.e. setting :ref:`perform_preprocessing <Basic pre-processing options>`, :ref:`perform_pdi <Basic PDI options>` and/or :ref:`perform_adi <Basic ADI options>` equal to ``False`` in the :ref:`configuration file <Configuration file>`), the sections pertaining to these skipped reduction steps are copied from the previous to the new log file and configuration file. This way the latest log file and configuration file always contain  all relevant information and input parameters of the reduction steps taken. 

After these initial steps, IRDAP will pre-process the OBJECT-files. To this end, it loads the static bad pixel map and the static master flat of the right filter, and creates a master sky frame from the provided SKY-files. The OBJECT-files are then sky (or background) subtracted, flat fielded, bad-pixel filtered and centered with a method chosen by the user. The centering would generally be performed using the CENTER-files, which are then processed accordingly.

.. note::
   Rather than using the static bad pixel map and master flat, the user can also create a master flat and bad pixel map by including a sequence of FLAT,LAMP- and DARK,BACKGROUND-files (or DARK-files) in the ``raw`` subdirectory (e.g. a sequence with exposure times 1, 2, 3, 4 and 5 s). The FLAT,LAMP-files preferably have the P0-90 polarizer set inserted, because it causes strong vignetting at the edges of the field of view which will otherwise not be corrected.
  
A cube of master flux frames is created for the left and right detector halves by processing the FLUX-files in a similar fashion as the OBJECT-files. If the data contains SKY-files with the same exposure time and neutral density filter as the FLUX-files, IRDAP processes these to subtract the sky background from the master flux frames. 

.. hint::
   IRDAP automatically determines the reference fluxes from the master flux frames and writes them to a CSV-file. These references fluxes can be used to convert the final images produced by IRDAP (e.g. the *I*\ :sub:`Q`- or *Q*:math:`_\phi`-images) from units of counts (ADU) into units of contrast/arcsec\ :sup:`2`. If the user can determine the stellar flux in Jansky, the final images can be expressed in Jansky/arcsec\ :sup:`2` (see the log file created by IRDAP for more details).

The pre-processed OBJECT-data is written to the subdirectory ``preprocessed`` and the processed SKY-, CENTER- and FLUX-data (and the user-created bad pixel map and master flat) to the subdirectory ``calibration``.

.. important::
   In case the pre-processed OBJECT- or FLUX-data is not correctly centered, IRDAP should be run again after adapting the :ref:`Advanced pre-processing options` of the configuration file.

PDI in a nutshell
-----------------

For the polarimetric differential imaging (PDI) part, IRDAP computes the double sum and double difference from the pre-processed OBJECT-files. It then applies the model-based correction method as described in `van Holstein et al. 2019 <ADS link>`_ to remove the instrumental polarization and cross-talk. The correction method for pupil-tracking observations has some differences compared to that of field-tracking observations (see `van Holstein et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017SPIE10400E..15V>`_). Subsequently, the background in the images is subtracted and the polarization of the star determined. The FITS-files of the final images are then written to two subdirectories: 

- ``reduced_pdi\no_star_pol_subtr``, containing the final images with the polarization of the star still present;
- ``reduced_pdi\star_pol_subtr``, containing the final images with the polarization of the star subtracted.

.. important::
   By default, the polarization of the star is determined with an annulus over the AO residuals. However, for the most accurate results the annulus should only contain signal from the star, and no signal from for example a circumstellar disk, companion or background star. Therefore the user often needs to adjust the input parameter :ref:`annulus_star <Basic PDI options>` in the configuration file.

.. important::
   The possibility to measure the polarization of the central star and the ability to create final images with and without this stellar polarization is a big advantage of IRDAP (see `van Holstein et al. 2019 <ADS link>`_). For images of a star or circumstellar disk for example, the stellar polarization can indicate the presence of an unresolved (inner) disk if it can be proven (or reasonably expected) that the polarization does not originate from interstellar dust. In that case one would use the images with the star polarization still present when making a comparison with radiative transfer models (see e.g. `Keppler et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018A&A...617A..44K>`_). Measuring the polarization of the star is also vital when measuring the polarization of substellar companions.
   
.. warning::
   For targets without a bright star (e.g. solar system objects), one would always use the images in the subdirectory ``reduced_pdi\no_star_pol_subtr``, i.e. those without the polarization subtracted.

.. note::
   The units of the *I*\ :sub:`Q`-, *I*\ :sub:`U`- and *I*\ :sub:`tot`-images are the number of counts (ADU) when summing the left and right frame halves of a single exposure (a single DIT). The images are averaged over the NDIT and the number of FITS-files used. Similarly, the units of the *Q*-, *U*-, *Q*:math:`_\phi`- and *U*:math:`_\phi`-images are the number of counts when subtracting the right from the left frame half of a single exposure. 

ADI in a nutshell
-----------------

.. attention::
   Angular differential imaging is not functional yet. It will be added around mid-October. 

Combining multiple data sets
----------------------------

If a target was observed using multiple observation blocks (OBs), it is recommended to first reduce each OB separately. After that, the final images of the PDI and ADI reductions can be mean-combined by using the terminal to navigate to a directory of your choice and typing:
::

   irdap --meancombine path1 path2 ... pathx
   
The space-separated ``paths`` are absolute paths to the main directories of the reductions, e.g.:
::

   irdap --meancombine /home/T_Cha_2016-02-20 /home/T_Cha_2016-02-21
   
The mean-combined images will be written to the current working directory of the terminal.

To understand the input parameters, continue with the :ref:`configuration file <Configuration file>`. 