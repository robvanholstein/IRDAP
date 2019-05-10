# irdis_polarimetry_pipeline
Pipeline for the accurate reduction of SPHERE-IRDIS polarimetric data

'''
path_main_dir:
    
String of path to location of data directory, e.g. 'C:\Data\T Cha'. The raw
data, i.e. FITS-files containing the OBJECT, CENTER, SKY and/or FLUX images, 
need to be present in a subdirectory 'Raw' within this directory, so for the 
example above in 'C:\Data\T Cha\Raw'. One can use the directory separator 
('\' above) specific to the operating system in use. The reduced calibration 
and object data will be written to subdirectories called 'Calibrations' and 
'Reduced', respectively.

path_static_flat_badpixelmap:

String of path to location of directory containing the FITS-files of the static
master flats and static bad pixel map, e.g. C:\Data\CentralFiles'. It might be 
useful to keep these static files in a central directory so that they are 
easily accessible when reducing various data sets.

save_preprocessed_data:
    
If True, save preprocessed cubes of single-sum and single-difference images in 
the 'preprocessed' folder so that the preprocessing can be skipped when 
re-running the pipeline.

double_difference_type:

Type of double difference to be computed, either 'standard' or 'normalized' 
(see van Holstein et al. 2019). In almost all cases one would use 'standard'.
When there are large variations in seeing and sky transparency among the
measurements, using 'normalized' can suppress spurious polarization signals and
improve the quality of the final images.

remove_vertical_band_detector_artefact:

If True remove the vertical band detector artefact seen in the double-
difference Q- and U-images by subtracting the median of the top and bottom 60
pixel of each pixel column in the images. If False don't remove the artefact.
In almost all cases one would use True. One would only use False in case 
there is a lot of astrophysical signal within in the region where the median 
is computed.

annulus_star: 
    
Parameter(s) defining with which annulus/annuli the star polarization will be
determined. If 'ao residuals' the annulus will be star-centered and located 
over the AO residuals. The inner radius and width of the annulus will depend on
the filter used. If 'star aperture' a small aparture located at the position of
the central star will be used. One can manually define the annulus/annuli by
specifying a (list of) length-6-tuple(s) of floats with parameters:
- x-coordinate of center (pixels, 1-based, center of image is at 512.5)
- y-coordinate of center (pixels, 1-based, center of image is at 512.5)
- inner radius (pixels)
- width (pixels)
- start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
- end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
For example, when using 'star aperture' the annulus used will be: 
(512.5, 512.5, 0, 11, 0, 360). The coordinates and angles of the annulus are
defined with respect to the final orientation of the image. Generally this is 
North up, except when single_posang_north_up = False (see below) for observations in 
field-tracking mode with a single derotator POSANG. The annulus used can be
checked with annulus_star.fits that can be found in the directories with the
reduced images.

annulus_background: 
    
Parameter(s) defining with which annulus/annuli the background will be
determined. If 'large annulus' the annulus will be star-centered and located 
far away from the star with an inner radius of 360 pixels and a width of 60
pixels. One can manually define the annulus/annuli by specifying a (list of) 
length-6-tuple(s) of floats with parameters:
- x-coordinate of center (pixels, 1-based, center of image is at 512.5)
- y-coordinate of center (pixels, 1-based, center of image is at 512.5)
- inner radius (pixels)
- width (pixels)
- start angle of annulus sector (deg; 0 due right and rotating counterclockwise)
- end angle of annulus sector (deg; 0 due right and rotating counterclockwise)
For example, when using 'large annulus' the annulus used will be: 
(512.5, 512.5, 360, 60, 0, 360). The coordinates and angles of the annulus are
defined with respect to the final orientation of the image. Generally this is 
North up, except when single_posang_north_up = False (see below) for observations in 
field-tracking mode with a single derotator POSANG. The annulus used can be
checked with annulus_background.fits that can be found in the directories with 
the reduced images.

combination_method_polarization:

Method to be used to produce the incident Q- and U-images, i.e. the images that
are corrected for the instrumental polarization effects. Valid values are 
'least squares', 'trimmed mean' or 'median'. The recommended option is 'trimmed
mean.' With 'least squares' the images are obtained by solving for every pixel 
the system of equations describing the measurements using linear least squares 
(see Eq. 35 of van Holstein et al. 2019). With 'trimmed mean' or 'median' the 
images are obtained by solving the system of equations for each pair of double-
difference Q- and U-images (each HWP cycle) separately, and then computing the 
trimmed mean or median over all resulting images. 'least squares' is the most 
accurate option, but any unremoved bad pixels will still be visible in the 
images. Using 'median' will remove these bad pixels, but is the least accurate 
option and also yields images with a lower signal-to-noise ratio as is clear
from images of circumstellar disks. Using 'trimmed mean' will yield images that
have essentially the same accuracy and signal-to-noise ratio images produced 
using 'least squares', but without the bad pixels. Therefore 'trimmed mean' is 
the recommended option.

trimmed_mean_prop_to_cut_polar:

Fraction to cut off of both tails of the distribution if 'trimmed mean' is used
for combination_method_polarization. Parameter is ignored in case 
'least squares' or 'median' is used. Value should be in range
0 <= trimmed_mean_prop_to_cut_polar <= 1. In most cases a 
value of 0.1 or 0.15 removes the bad pixels well while producing images very 
similar to those obtained with 'least squares'.
    
combination_method_intensity:

Method to be used to produce the incident I_Q- and I_U-images. These images are 
computed by combining the I_Q- or I_U-images of all HWP cycles using the 
'mean', 'trimmed mean' or 'median'. 'mean' yields the most accurate images, but
any unremoved bad pixels will still be visible in the images. Using 'median' 
will remove these bad pixels, but is the least accurate option. 'trimmed mean'
produces images similar to 'mean', but without the bad pixels. It is generally
recommended to use either 'trimmed mean' or 'mean'.

trimmed_mean_prop_to_cut_intens:

Fraction to cut off of both tails of the distribution if 'trimmed mean' is used
for combination_method_intensity. Parameter is ignored in case 
'mean' or 'median' is used. Value should be in range
0 <= trimmed_mean_prop_to_cut_intens <= 1. In most cases a 
value of 0.1 or 0.15 removes the bad pixels well while producing images similar
to those obtained with 'mean'.

single_posang_north_up:

For observations taken in field-tracking mode with a single derotator position 
angle (INS4.DROT2.POSANG either 0 or with a fixed offset), the final images 
are rotated with North up if True, and kept in the orientation of the raw 
frames if False. In the latter case, the images are more accurate as they 
suffer less from interpolation errors. This is useful when for example 
extracting the polarized surface brightness distribution of a circumstellar
disk. Parameter is ignored for field-tracking observations with more than one 
derotator position angle or observations taken in pupil-tracking mode, because 
in these cases the final images produced always have North up.

normalized_polarization_images: 
    
if True create final images of degree of linear polarization, normalized Stokes
q and u and degree and angle of linear polarization computed from q and u:
- DoLP = sqrt(Q^2 + U^2) / (0.5*(I_Q + I_U))
- q = Q / I_Q
- u = U / I_U
- AoLP_norm = 0.5 * arctan(u / q)
- DoLP_norm = sqrt(q^2 + u^2)
These images are only valid if all flux in the images originates from the 
astrophysical source of interest. This is generally the case for observations 
of for example solar system objects or galaxies. The images are generally not 
valid for observations of circumstellar disks or companions because in that 
case a large part of the flux in the total intensity images originates from the 
central star. AoLP_norm and DoLP_norm are potentially more accurate than 
AoLP = 0.5 * arctan(U / Q) and DoLP, especially when there are significant 
variations in seeing and sky transparency among the measurements.


'''