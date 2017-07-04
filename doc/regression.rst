.. _pySpline_regression:

Regression Tests
----------------

:ref:`pySpline` includes a set of built-in regression tests that are
used to ensure that the code reproduces a consistent set of
reproducible results. It is *HIGHLY* recommended that the regression
tests are run after installation on a new system to ensure consistency
of the code build.

To run the regression tests go to the following directory::

  python/reg_tests

and run the automatic testing script::

  python run_reg_tests.py

If successful, you should see the following line::
 
  pyspline: Success!

if you don't, then *something* went wrong in the test. By default an
``xxdiff`` window will appear that should quickly highlight where the
discrepency is. 

After running the regression tests, there will be *three* important
files to look at::

   pyspline_reg.ref
   pyspline_reg
   pyspline_ref.orig

* ``pyspline_reg.ref`` is the reference file that contains the
   'truth' values that we should reproduce

* ``pyspline_ref`` is the generated output that is *identical* to the
  reference output *except* for @value lines where a discrepency has
  been detected. A graphical diff of ``pyspline_reg`` and
  ``pyspline_reg.ref`` should quickly show where the discrepancies
  are.
   
* ``pyspline_ref.org`` is the *actual* output of the regression
  run. Generally speaking there will be small discrepancies in the 15
  or 16th digit which depends on the compiler. This will trigger a
  text-based diff program which would make it difficult to determine
  which values are incorrect in a file with several thousand lines. 

.. note::

  Only lines that contain the @value magic string are actually
  compared by the regression tester. 


Regression Test Script Help 
+++++++++++++++++++++++++++
To see a list of options for the regression test script run::

  python run_reg_test.py --help

which will show::

  usage: run_reg_tests.py [-h] [-mode {train,compare}] [-diff_cmd DIFF_CMD]
                        [-nodiff]

  optional arguments:
    -h, --help            show this help message and exit
    -mode {train,compare}
                        Train generates reference file. Compare runs test
    -diff_cmd DIFF_CMD    Command to run for displaying diff. Default: xxdiff
    -nodiff               Suppress displaying the comparison if not successful
