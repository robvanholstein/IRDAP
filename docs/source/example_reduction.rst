
Example Reduction
=================

As an introduction to the user, but also to test if IRDAP really works, a quick 
example reduction can be performed. This demo uses (a limited amount of )raw data of 
the circumstellar disk of T Cha from program `PROGRAM ID <link>`_ as published in
`Pohl et al. (2017) <http://adsabs.harvard.edu/abs/2017A%26A...605A..34P>`_ and
as used in `van Holstein et al. (2019) <ADS link>`_ to exemplify the correction
method.

To run the demo, navigate to a directory of your choice in the terminal and enter:
::

   irdap --demo

IRDAP will then download 48.5 MB of raw data from its `GitHub repository <>`_ and perform
a complete data reduction. The resulting files will be written to the working
directory of the terminal.

If you open the subdirectory reduced_star_pol_subs Qphi, Uphi, (linear scale, -300 to 150, cool)
you will see the following images. These images look a lot like those of Figure 14 (left column)
of van Holstein et al. (2017), but with a lower signal-to-noise.

.. image:: ./figs/t_cha_quphi_combined.svg
    :width: 750px
    :align: center

Congratulations, you have reduced your first data set with IRDAP!

To reduce your own data set, please continue to the :ref:`User Instructions`.


- Make image of Qphi and Uphi
- Example include flux too
- In text mention that center, sky, flux used
- Finish text