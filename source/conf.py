# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'EMUMADZ'
# Enhanced MUtation MApping and Detection in Zebrafish
copyright = '2025, Tyrone Chen, Richard Lupat, Michelle Meier, Maia Zethoven, Greg Baillie, Scott Paterson, Oguzhan Baltaci, Cas Simons, Jason Li, Benjamin Hogan'
author = 'Tyrone Chen, Richard Lupat, Michelle Meier, Maia Zethoven, Greg Baillie, Scott Paterson, Oguzhan Baltaci, Cas Simons, Jason Li, Benjamin Hogan'
release = '0.1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'
    ]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
