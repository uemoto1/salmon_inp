#!/usr/bin/env python3
import os
import sys
from setuptools import setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    os.system("python setup.py bdist_wheel upload")
    print("You probably want to also tag the version now:")
    print("  git tag -a VERSION -m 'version VERSION'")
    print("  git push --tags")
    sys.exit()

setup(
    name="salmon_inp",
    version="0.0.1",
    author="Mitsuharu Uemoto",
    description="Convert CIF to SALMON inputfile",
    url="https://github.com/uemoto1/salmon_inp",
    py_modules=["salmon_inp"],
    install_requires=[
        "numpy>=1.1.0",
        "scipy",
    ],
    entry_points="""
        [console_scripts]
        cif2salmon=cif2salmon:main
    """,
)
