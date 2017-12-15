#!/usr/bin/python3
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='confdist',
    version='0.1.0',
    description='Compute transcriptional edit distance, or confusion distance, between pairs of words.',
    url='https://github.com/balazsdukai/confdist',
    author='Balázs Dukai, Vincent van Altena',
    author_email='balazs.dukai@gmail.com, vpvanaltena@live.com',
#    license='GPLv3',
    packages=['confdist'],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',
        
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: GIS',
        
        # Pick your license as you wish (should match "license" above)
#         'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        
        'Operating System :: POSIX :: Linux'
    ],
    python_requires='>=3',
#    keywords='GIS 3DGIS CityGML LiDAR',
    entry_points={
        'console_scripts': ['confdist = confdist.confdistapp:main']
    },
    include_package_data=True,
    zip_safe=False
    )
