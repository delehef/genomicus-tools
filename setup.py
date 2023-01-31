#!/usr/bin/env python3
from distutils.core import setup

setup(name='GenomicusTools',
      version='1.0.0',
      description='Tools to create Genomicus databases',
      author='DYOGEN lab',
      author_email='agora@bio.ens.psl.eu',
      url='https://github.com/DyogenIBENS/Agora/',
      packages=['src', 'utils'],
      package_dir= {
          'agora': 'src/',
          'utils': 'src/utils',
      },
      scripts=[
          'src/createTable.BlocksStats.py',
          'src/createTable.CNE-CNEitems.py',
          'src/createTable.Gene-Search.py',
          'src/createTable.Species.py',
          'src/createTable.Species_new.py',
          'src/createTable.SpeciesTree.py',
          'src/createTable.Synteny.py',
          'src/createTable.Tree-Orthologs.py',
      ]
     )
