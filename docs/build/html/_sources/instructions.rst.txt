
Usage instructions
==================

Running IRDAP
-------------

Reducing data with IRDAP is very straightforward and only requires two commands to be entered in the terminal. First, create a directory (e.g. :file:`/home/T_Cha_2016-02-20`) containing a subdirectory called :file:`raw` (i.e. :file:`/home/T_Cha_2016-02-20/raw`) in which you place the raw FITS-files. 

.. attention::
   The raw FITS-files should have `DPR TYPE` equal to OBJECT, OBJECT,CENTER, OBJECT,FLUX or SKY. If the user wishes to make a new master flat and bad pixel map (see :ref:`Pre-processing in brief`) `DPR TYPE` can also be equal to FLAT,LAMP, DARK,BACKGROUND or DARK.

.. hint::
   It is recommended to download the raw data from the `SPHERE data archive <http://archive.eso.org/wdb/wdb/eso/sphere/form>`_ and when doing so setting the search field `DPR TECH` to `any POLARIMETRY`. This way only data that is required for IRDAP will be downloaded.
   
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
 
Pre-processing in brief
-----------------------

The data reduction is divided in a pre-processing and a post-processing part. 

.. warning::
   The data-reduction steps are only described concisely here. The log file that IRDAP creates describes the steps in more detail and also briefly describes the output files created. It is therefore recommended to carefully read the log file the first time you run IRDAP.
 
The pre-processing starts with reading the configuration file and creating the log file that contains all the lines printed on screen. In addition, IRDAP creates overviews of the relevant header keywords for each file in the raw subdirectory, sorts the raw data and creates subdirectories to write the output files to.

.. hint::
   The header overviews can also be created without running IRDAP by typing :code:`irdap --headers` in the terminal.
   
.. note::
   When running IRDAP a copy of the configuration file is made. When subsequently re-running IRDAP, this copy is moved together with the log file to the ``logs`` subdirectory as a record of the performed reductions. 

After these initial steps, IRDAP will pre-process the OBJECT-files. To this end, it loads the static bad pixel map and the static master flat of the right filter, and creates a master sky frame from the provided SKY-files. The OBJECT-files are then sky (or background) subtracted, flat fielded, bad-pixel filtered and centered with a method chosen by the user. The centering would generally be performed using the CENTER-files, which are then processed accordingly.

.. note::
   Rather than using the static bad pixel map and master flat, the user can also create a master flat and bad pixel map by including a sequence of FLAT,LAMP- and DARK,BACKGROUND-files (or DARK-files) in the ``raw`` subdirectory (e.g. a sequence with exposure times 1, 2, 3, 4 and 5 s). The FLAT,LAMP-files need to have the P0-90 polarizer set inserted.
  
A master flux frame is created by processing the FLUX-files in a similar fashion as the OBJECT-files. If the data contains SKY-files with the same exposure time and neutral density filter as the FLUX-files, IRDAP processes these to subtract the sky background from the master flat frame. 

The pre-processed OBJECT-data is written to the subdirectory ``preprocessed`` and the processed SKY-, CENTER- and FLUX-data (and the user-created bad pixel map and master flat) to the subdirectory ``calibration``.

Post-processing in brief
------------------------

For the post-processing part, IRDAP computes the double sum and double difference from the pre-processed OBJECT-files. It then applies the model-based correction method as described in `van Holstein et al. 2019 <ADS link>`_ to remove the instrumental polarization and cross-talk. The correction method for pupil-tracking observations has some differences compared to that of field-tracking observations (see `van Holstein et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017SPIE10400E..15V>`_). Subsequently, the background in the images is subtracted and the polarization of the star determined. The FITS-files of the final images are then written to two subdirectories: 

- ``reduced``, containing the final images with the polarization of the star still present;
- ``reduced_star_pol_subtr``, containing the final images with the polarization of the star subtracted.

.. important::
   The possibility to measure the polarization of the central star and the ability to create final images with and without this stellar polarization is a big advantage of IRDAP (see `van Holstein et al. 2019 <ADS link>`_). For images of a star or circumstellar disk for example, the stellar polarization can indicate the presence of an unresolved (inner) disk if it can be proven (or reasonably expected) that the polarization does not originate from interstellar dust. In that case one would use the images with the star polarization still present when making a comparison with radiative transfer models. Measuring the polarization of the star is also vital when measuring the polarization of substellar companions.
   
.. warning::
   For targets without a bright star (e.g. solar system objects), one would always use the images in the subdirectory ``reduced``, i.e. those without the polarization subtracted.

To understand all the input parameters, continue with the :ref:`configuration file <Configuration file>`. 