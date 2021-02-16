#!/usr/bin/env python3

import setuptools

setuptools.setup(
name='RiboridDesign',
version='0.1.0',
author='Saugat Poudel',
author_email='sapoudel@ucsd.edu',
description='Python package for designing oligos for riborid protocol',
maintainer='Saugat Poudel',
url='https://github.com/SBRG/RiboRid_Design',
packages=setuptools.find_packages(),
python_requires=">3.7.9",
include_package_data=True,
install_requires=open('requirements.txt').read(),
platforms="GNU/Linux, Mac OS X > 10.7"
)
