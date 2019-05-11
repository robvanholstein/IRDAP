
Installation
============

fcmaker is available on pypi, which makes its installation easier than ever. 
In a terminal, type:
::

   pip install fcmaker

And that should take care of things.

The most recent release of fcmaker is also available for download from its `Github repository <https://github.com/fpavogt/fcmaker/releases/latest/>`_. 
Interested users can fork the fcmaker repository if they want to get access to the 
latest updates not yet released. *Push requests* for bug fixes and new features are 
welcome and will be examined in detail. 
      
Requirements
------------
fcmaker is written in Python 3.6. The following packages are required for it to work 
properly:

* aplpy (1.1.1 or above)
* astroplan
* astropy (3.0 or above)
* astroquery (0.3.4 or above)
* matplotlib (2.0.2 or above)
* numpy (1.13.1 or or above)
* p2api (0.92 or above)
* pillow (4.2.1 or above)
* pytz (2018 or above)
* PyYAML (3.12 or above)
* scipy (0.19.0 or above)

Optional: 

* Montage and montage-wrapper (0.9.9 or above, only until the next release of aplpy)
* Proper system-wide LateX installation (allows for prettier plots)

The Montage package is required to rotate the finding charts North, even if the underlying
FITS file isn't. It won't be required once aplpy gets upgraded.

Testing the installation
------------------------

In a terminal shell, try to access the basic help of fcmaker::
 
   python -m fcmaker --help
 
If that works, chances are, you will probably be fine. Note that fcmaker requires a 
connection to the internet to work (even in local mode).

 