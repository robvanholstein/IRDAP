
Usage instructions
==================

There are two main ways to use fcmaker. You can either execute the module as a script, i.e.::

   python -m fcmaker
   
or you can import the module and execute it within a script, i.e.::

   >>> import fcmaker
   >>> fcmaker.make_fc( ... ) # or fcmaker.make_fc_local ( ... )
   
When running fcmaker as a script, any argument you feed it is pretty much sent to the 
underlying functions ``make_fc()`` and ``make_fc_local()``. For simplicity, this page
only discusses how to run the entire module as a script, which ought to be slightly more 
friendly to users not (yet!) familiar with Python. Still, the hope is that after reading 
this page, the use of the functions ``make_fc()`` and ``make_fc_local()`` shouldn't be
too mysterious. Also, don't forget the built-in help::

   >>> import fcmaker
   >>> help(fcmaker.make_fc)
 