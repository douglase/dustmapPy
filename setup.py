"""
edouglas@mit.edu

Modified from scipy-lectures.org:
http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id10
CC Attribution license: https://github.com/scipy-lectures/scipy-lecture-notes/blob/master/LICENSE.rst

to build: 
>>> python setup.py build_ext --inplace
"""

from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext


setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("dustmap",
                 sources=["_dustmap.pyx", "dustmap/dustmap.c"],
                 include_dirs=[numpy.get_include()])],
)

'''
from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("cos_doubles",
                 sources=["_cos_doubles.pyx", "cos_doubles.c"],
                 include_dirs=[numpy.get_include()])],
)'''


'''#!/usr/bin/env python

from distutils.core import setup, Extension

import numpy
from Cython.Distutils import build_ext

setup(name='jsv_rt_cython',
      version='0.1',
      description='',
      requires = ['numpy', 'cython'],
      author='E.S. Douglas, G.Geddes and J.S. Vickers',
      url='https://bitbucket.org/douglase/jsv_rt_cython',
      py_modules=['OO_inversion', 'raids_functions','InversionFunctions','RAIDS_auric_funcs','fit_chapman_alpha','tec'],
      cmdclass={'build_ext': build_ext},
      ext_modules=[
          Extension("gsrp",
              sources=["RT/_gsrp.pyx", "RT/gsrp.c","RT/Intensities.c","Matrix/matrix.c"],
              include_dirs=[numpy.get_include()]),
          Extension("multi",
              sources=["RT/_multi.pyx", "RT/multi.c","RT/Intensities.c","Matrix/matrix.c"],
              include_dirs=[numpy.get_include()])
          ],
      ) #equivalent to  "gcc -g -o gsrp gsrp.c Intensities.c ../Matrix/matrix.c -lm"
'''
