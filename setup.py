import sysconfig
import sys
import os
import subprocess
import re
from setuptools import setup, find_packages

setup(
    # package information
    name='ocms_shotgun',
    version="1.0.0",
    description='OCMS_Shotgun : Oxford Centre for Microbiome Studies pipelines for shotgun metagenome processing',
    author='Sandi Yen, Nicholas Ilott, Jethro Johnson',
    license="MIT",
    platforms=["any"],
    keywords="shotgun, metagenomics",
    url="https://github.com/OxfordCMS/OCMS_Shotgun",
    #packages=find_packages("./") + find_packages("./ocmsshotgun/"),
    packages=find_packages(),
    entry_points={
        'console_scripts': ['ocms_shotgun = ocmsshotgun.ocms_shotgun:main']
    },
    include_package_data=True,
    python_requires='>=3.8.2'                                            
)
