from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup

setup(
       name = 'pcmap',
       version = '0.1',
       license='GPL',
       description = 'Computing contact map for protein structures',
       author = 'Guillaume Launay',
       author_email = 'guillaume.launay@ibcp.fr',
       url = 'https://github.com/MMSB-MOBI/pcmap', # use the URL to the github repo
       package_dir={'': 'src'},
       packages=['pcmap'],
       zip_safe=False,
       py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
       download_url = 'https://github.com/glaunay/pcmap/tarball/0.1', # I'll explain this in a second
       keywords = ['protein', 'structure', 'atomic contact'], # arbitrary keywords       
       install_requires=[
              'ccmap', 'docopt', 'pypstruct'
       ],
       classifiers=[
              "Intended Audience :: Science/Research",
              "License :: OSI Approved :: MIT License",
              "Programming Language :: Python :: 3",
              "Topic :: Scientific/Engineering :: Bio-Informatics",
              "Topic :: Software Development :: Libraries :: Python Modules",
              "Development Status :: 5 - Production/Stable" 
       ]
)