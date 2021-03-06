from distutils.core import setup, Extension
from glob import glob
from os.path import dirname, realpath
import numpy as np

sources = [x for x in glob(dirname(realpath(__file__)) + '/src/*.c', recursive=True)
           if not x.endswith('main.c')]
py_sources = glob(dirname(realpath(__file__)) + '/python/*.c')
extras = ['-Wall', '-Wextra', '-Wconversion', '-Wno-unused-variable',
          '-Wno-unused', '-pedantic', '-Wmissing-prototypes',
          '-Wstrict-prototypes', '-O3']

module = Extension('qop', sources=sources + py_sources,
                   include_dirs=[
                       np.get_include(),
                       dirname(realpath(__file__))
                   ],
                   extra_compile_args=extras)

setup(name='qop',
      version='1.0',
      description='A variational quantum circuit simulator.',
      ext_modules=[module])
