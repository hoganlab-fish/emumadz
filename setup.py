#!/usr/bin/python
from setuptools import setup

setup(
    name="emumadz",
    version="0.2.0",
    packages=["emumadz"],
    package_dir={"../": "emumadz"},
    entry_points={
        'console_scripts': [
            'emumadz=emumadz.parse_vcf:main',
            'emumadz=emumadz.plot_homozygosity:main',
        ],
    },
    python_requires=">=3.9",
    install_requires=[
        "cryptography",
        "flask",
        "numpy",
        "pandas",
        "pysam",
        "tqdm",
    ],
    author="Tyrone Chen",
    author_email="tyrone.chen@petermac.org",
    description="Enhanced MUtation MApping and Detection in Zebrafish",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/hoganlab-fish/emumadz",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)