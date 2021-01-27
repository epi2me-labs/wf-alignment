"""
Mapping Stats setup script

Copyright (c) 2021 by Oxford Nanopore Technologies Ltd.
"""
import os
from setuptools import setup, find_packages

os.environ['GIT_SSL_NO_VERIFY'] = 'true'

__version__ = '1.0.0'

setup(
    name='mapping_stats',
    version=__version__,
    author='epi2melabs',
    setup_requires=['pytest-runner'],
    description='Mapping Stats',
    zip_safe=False,
    install_requires=[
        'pysam',
        'pandas',
        'aplanat'
    ],
    packages=find_packages(exclude=("tests",)),
    package_data={
        'mapping_stats': ['mapping_stats/tests/data/*'],
    },
    entry_points={
        "console_scripts": [
            'gather_mapping_stats = mapping_stats.gather:run_main',
            'combine_mapping_stats = mapping_stats.combine:run_main',
            'visualise_mapping_stats = mapping_stats.visualise:run_main',
        ]
    }
)
