#!/usr/bin/env python
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
# -------------------------------------
from setuptools import setup

setup(name='feon',
      version='1.0.0',
      description='FEA python-based framework',
      author='YAOYAO PEI',
      author_email='yaoyao.bae@foxmail.com',
      url = 'https://github.com/YaoyaoBae/Feon.git',
      license = "HBUT",
      packages=['feon','feon.sa'],
      install_requires=['numpy'],
      
      )
