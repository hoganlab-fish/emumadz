#!/usr/bin/python
from setuptools import setup, find_packages

setup(
    name="emumadz",
    version="0.2.0",
    packages=find_packages(),
    package_dir={"": "emumadz"},
    entry_points={
        'console_scripts': [
            'emumadz=emumadz.parse_vcf:main',
        ],
    },
    python_requires=">=3.9",
    author="Tyrone Chen",
    author_email="tyrone.chen@petermac.org",
    description="Enhanced MUtation MApping and Detection in Zebrafish",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/hoganlab-fish/whole_genome_sequencing",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)