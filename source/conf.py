# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath("../"))


project = 'EMUMADZ'
# Enhanced MUtation MApping and Detection in Zebrafish
# '2025, Tyrone Chen, Richard Lupat, Michelle Meier, Maia Zethoven, Greg Baillie, Scott Paterson, Oguzhan Baltaci, Cas Simons, Jason Li, Benjamin Hogan'

release = '0.2.0'

authors = [
    ("Tyrone Chen",     "0000-0002-9207-0385"),
    ("Richard Lupat",   "0000-0002-6435-7100"),
    ("Michelle Meier",  "0009-0005-5595-3882"),
    ("Maia Zethoven",   "0000-0002-6528-891X"),
    ("Greg Baillie",    "0000-0002-6130-250X"),
    ("Scott Paterson",  "0000-0000-0000-0000"),
    ("Oguzhan Baltaci", "0009-0001-5651-1331"),
    ("Cas Simons",      "0000-0003-3147-8042"),
    ("Jason Li",        "0000-0002-1150-3549"),
    ("Benjamin Hogan",  "0000-0002-0651-7065"),
]

def add_orcid(author_name, orcid):
    return f'<a href="https://orcid.org/{orcid}">{author_name} <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, '

# Set author and copyright fields
author = ",\n".join([x for x, y in authors])
copyright = ", ".join(["2025", author])

# HTML-specific configuration
html_theme_options = {
    'style_external_links': True,
}

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode'
    ]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
