.. _pySpline_install:

Installation
============

Building
--------

For speed purposes, pySpline uses a small compiled Fortran library for doing the time consuming computational operations.
It is therefore necessary to build this library before using pySpline.

pySpline follows the standard MDO Lab build procedure.
To start, find a configuration file close to your current setup in::

    config/defaults

and copy it to ``config/config.mk``. For example:

.. prompt:: bash

    cp config/defaults/config.LINUX_GFORTRAN.mk config/config.mk

If you are a beginner user installing the packages on a linux desktop,
you should use the ``config.LINUX_GFORTRAN.mk`` versions of the configuration
files. The ``config.LINUX_INTEL.mk`` versions are usually used on clusters.
Once you have copied the config file, compile pySpline by running:

.. prompt:: bash

    make

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module libspline can be imported...
   Module libspline was successfully imported.

If you don't see this, it will be necessary to configure the build manually.
To configure manually, open ``config/config.mk`` and modify options as necessary.

Lastly, to build and install the Python interface, type:

.. prompt:: bash

    pip install .

Verification
------------
To verify the library, pySpline contains a set of tests that can be run automatically to ensure it reproduces the expected reference results.
To do so, testing dependencies need to be installed first, by typing:

.. prompt:: bash

    pip install .[testing]

Once testing dependencies are installed, then to execute all tests, run the following in the root directory,

.. prompt:: bash

    testflo .
