#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function, unicode_literals
from setuptools import find_packages, setup


version = '0.9.5'


setup_args = dict(
    name             = 'GWINC',
    version          = version,
    url              = 'https://git.ligo.org/gwinc/pygwinc',
    author           = 'LIGO Laboratory',
    author_email     = 'jrollins@ligo.caltech.edu ',
    description      = "Gravitation Wave Interferometer Noise Calculator",
    license          = 'Copyright 2017 LIGO Laboratory',
    install_requires = [
        'numpy',
        'scipy',
        'matplotlib',
    ],
    packages = find_packages(
        exclude = ['docs',],
    ),
    include_package_data = True,
    zip_safe = True,
    keywords = 'Noise, LIGO, Gravitational Wave,',
    classifiers=[
        'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
)

if __name__ == "__main__":
    setup(**setup_args)
