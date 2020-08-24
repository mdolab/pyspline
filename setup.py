from setuptools import setup
import re

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('pyspline/__init__.py').read(),
)[0]

setup(name='pyspline',
      version=__version__,


      description="pySpline is a package containing code for working with b-spline curve, surface, and volume objects",
      long_description="""Please see the [documentation](http://mdolab.engin.umich.edu/docs/packages/pyspline/doc/index.html) for installation details and API documentation.

      To locally build the documentation, enter the `doc` folder and enter `make html` in terminal.
      You can then view the built documentation in the `_build` folder.


      Citation
      --------

      Please cite pySpline in any publication for which you find it useful.
      For more background, theory, and figures, see the [pyGeo journal article](http://mdolab.engin.umich.edu/sites/default/files/mao2010_final.pdf).

      G. K. W. Kenway, Kennedy, G. J., and Martins, J. R. R. A., “A CAD-Free Approach to High-Fidelity Aerostructural Optimization”, in Proceedings of the 13th AIAA/ISSMO Multidisciplinary Analysis Optimization Conference, Fort Worth, TX, 2010.

      @conference {Kenway:2010:C,
        title = {A {CAD}-Free Approach to High-Fidelity Aerostructural Optimization},
        booktitle = {Proceedings of the 13th AIAA/ISSMO Multidisciplinary Analysis Optimization Conference},
        year = {2010},
        note = {AIAA 2010-9231},
        month = {September},
        address = {Fort Worth,~TX},
        author = {Gaetan K. W. Kenway and Graeme J. Kennedy and Joaquim R. R. A. Martins}
      }
      """,
      long_description_content_type="text/markdown",
      keywords='spline b-spline optimization',
      author='',
      author_email='',
      url='https://github.com/mdolab/pyspline',
      license='Apache License Version 2.0',
      packages=[
          'pyspline',
      ],
      package_data={
          'pyspline': ['*.so']
      },
      install_requires=[
            'numpy>=1.16',
            'scipy>=1.2',
      ],
      classifiers=[
        "Operating System :: Linux",
        "Programming Language :: Python, Fortran"]
      )

