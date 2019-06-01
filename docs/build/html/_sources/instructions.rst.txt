
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