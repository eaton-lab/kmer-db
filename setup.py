#!/usr/bin/env python

import re
from setuptools import setup

# parse version from init.py
with open("kmunity/__init__.py") as init:
    CUR_VERSION = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]",
        init.read(),
        re.M,
    ).group(1)

# setup installation
setup(
    name="kmunity",
    packages=["kmunity"],
    version=CUR_VERSION,
    author="Deren Eaton",
    author_email="de2356@columbia.edu",
    install_requires=[
        "pandas",
        "requests",
        "loguru",
    ],
    entry_points={
        'console_scripts': ['kmunity = kmunity.__main__:main']},
    license='GPL',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ],
)
