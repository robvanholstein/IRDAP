
Example reduction
=================

As an introduction to the user, but also to check if IRDAP is installed correctly, a quick 
end-to-end example reduction can be performed. This demo uses a small portion of the H-band 
polarimetric data of the circumstellar disk of T Cha from ESO Program 096.C-0248(C). This
data set was the first to be reduced with IRDAP's correction method and is published in
`Pohl et al. (2017) <http://adsabs.harvard.edu/abs/2017A%26A...605A..34P>`_. It is also 
used in `van Holstein et al. (2019) <ADS link>`_ to exemplify the correction method.

To run the demo, use the terminal to navigate to a directory of your choice and type:
::

   irdap --demo

IRDAP will then perform a complete data reduction and will write the resulting files
to the terminal's current directory. The produced files are explained in 
:ref:`Usage Instructions`.

.. note::
   An internet connection is required to run the demo, because 56.6 MB of raw data will
   be downloaded from the 
   `IRDAP GitHub repository <https://github.com/robvanholstein/IRDAP/tree/master/irdap/example_data>`_.

After running the demo, you can navigate to the :file:`reduced` directory and open the 
files :file:`T_Cha_2016-02-20_Q_phi.fits` and :file:`T_Cha_2016-02-20_U_phi.fits` that
show the final :math:`Q_\phi`- and :math:`U_\phi`-images.
Opening these files with `DS9 <http://ds9.si.edu/>`_, setting color to cool, scale to linear 
and the limits to -500 to 1500 and -300 to 150, respectively, you will see the following:

.. figure:: ./figs/t_cha_quphi_combined.svg
    :width: 750px
    :align: center
	
These images correspond to Figure 14 (left column) of `van Holstein et al. (2019) <ADS link>`_, 
but with lower signal-to-noise because for the demo the data of only 1 instead of 30 HWP cycles 
was used.

To start reducing your own data with IRDAP, please continue to the :ref:`Usage instructions`.