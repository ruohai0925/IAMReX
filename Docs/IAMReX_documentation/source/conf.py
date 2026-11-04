# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
import subprocess
import sphinx_rtd_theme

# build the doxygen documentation
subprocess.call('cd ../doxygen; make html', shell=True)

project = 'IAMReX'
copyright = '2026, IAMReX Team'
author = 'IAMReX Team'
release = '0.1.0'

bibtex_bibfiles = ["refs.bib"]
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}
main_doc = 'index'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.mathjax',
              'sphinxcontrib.bibtex',
              'sphinx.ext.githubpages',
              'sphinx.ext.viewcode',
              'sphinx.ext.intersphinx',
              'sphinx.ext.autosectionlabel',
              'breathe'
              ]

breathe_projects = { 'IAMReX': '../doxygen/xml' }

templates_path = ['_templates']
exclude_patterns = []

# -- Options for equation numbering ------------------------------------------
numfig = True
math_numfig = True
numfig_secnum_depth = 2
numfig_format = {'table': 'Table %s', 'figure': 'Figure %s', 'code-block': 'Listing %s', 'section': 'Section'}

mathjax3_config = {
    'tex': {
        'tags': 'ams',
        'tag_side': 'right',
    }
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
# html_static_path = ['_static']
