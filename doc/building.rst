.. _pySpline_building:

Building
--------

For speed purposes, ``pySpline`` uses a small compiled Fortran library for doing the time consuming computational operations. 
It is therefore necessary to build this library before using ``pySpline``.

pySpline follows the standard MDO Lab build procedure.
To start, find a configuration file close to your current setup in::

    $ config/defaults

and copy it to ''config/config.mk''. For example::

    $ cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk

If you are a beginner user installing the packages on a linux desktop, 
you should use the ``config.LINUX_GFORTRAN.mk`` versions of the configuration 
files. The ``config.LINUX_INTEL.mk`` versions are usually used on clusters.

Once you have copied the config file, compile `pySpline`` by running::

    $ make

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module pyspline can be imported...
   Module pyspline was successfully imported.

If you don't see this, it will be necessary to configure the build manually.
To configure manually, open ``config/config.mk`` and modify options as necessary.

Lastly, to build and install the Python interface, type::

    pip install .

