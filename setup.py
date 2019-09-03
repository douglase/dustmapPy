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
