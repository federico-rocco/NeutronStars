# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 11:46:42 2021

@author: cosmo
"""

from setuptools import setup, find_packages

with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = f.read()
    
setup(
      name='myproject', 
      version='1.0', 
      packages=find_packages(),
      install_requires=requirements
      )
