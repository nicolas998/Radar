#!/usr/bin/env python
import os
from numpy.distutils.core import setup, Extension

ext1 = Extension(name = 'radar_f90',
	sources = ['radar/fort_func.f90'])

setup(
    name='radar',
    version='0.1.0',
    author='Nicolas Velasquez G',
    author_email='nicolas.velasquezgiron@gmail.com',    
    packages=['radar'],
    package_data={'radar':['radar_f90.so','ajuste_multicapaall_77.pkl']},
    data_files=[('', ['LICENSE.txt'])],
    url='http://pypi.python.org/pypi/Radar/',
    license='LICENSE.txt',
    description='.',
    long_description=open('README.txt').read(),
    ext_modules=[ext1],
	)

