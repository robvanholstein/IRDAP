
Usage instructions
==================

Reducing data with IRDAP is very straightforward and only requires
two commands to be entered in the terminal. First, create a directory 
(e.g. :file:`/home/T_Cha_2016-02-20`) containing a subdirectory called 
:file:`raw` in which you place the raw FITS-files. Then in the terminal 
navigate to the directory (e.g. ``cd /home/T_Cha_2016-02-20``) and type:
::
 
   irdap --makeconfig

This creates a default configuration file :file:`config.conf` in the main directory. 
You can then adjust the parameters in the configuration file with a text 
editor (see :ref:`Configuration file` for an overview of what each parameter does). 
However, at first often no changes to the configuration file are required.
Finally in the terminal type:
::

   irdap --run

to perform the data reduction.
 
Output of IRDAP
---------------
 
 
	1. Usage instructions 1) only two commands, makeconfig  and run commands from current working directory, adapt config, but for majority will work with default config file, explain output directories. Explain other commands if they exist.
	2. 
Which raw file types, what pipeline does briefly [keep short!] (rest in log file) in bulleted list pre-processing (+ where output): read config file, create log on screen and as file, sort data, reads a static BPM and flat, but can make own, create a master flat frame, center frames, object frames, sky for flux, flux frames. Save files to enable skipping pre-processing. Post-processing: double difference, apply model based correction method to remove IP, derotate images and correct for cross-talk and rotation. Different for ft and PT. Subtract background and determine star polarization. With and without star polarization, output fits images.


After downloading the data from the `SPHERE data archive <http://archive.eso.org/wdb/wdb/eso/sphere/form>`_,

sort data
directories
log
which steps
file output
 
What will it do (directories, steps, output)


**WORK: make hyperlink of header overview to frames_to_remove on configuration file section**

It is also possible to create an overview of the relevant headers of the 
FITS-files in the raw subdirectory without running the pipeline by typing:
::  

   irdap --headers
   
Refer to configuration file section