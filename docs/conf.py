# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.append(os.path.abspath('../src'))
from gdt.missions.burstcube import __version__


# -- Project information -----------------------------------------------------

project = 'The Gamma-ray Data Tools'
license = 'MIT'
author = 'The BurstCube Team'

# The full version, including alpha/beta/rc tags
version = __version__
release = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.viewcode', 'sphinx.ext.inheritance_diagram',
              'sphinx.ext.mathjax', 'sphinx.ext.napoleon', 'nbsphinx',
              'IPython.sphinxext.ipython_console_highlighting',
              'sphinx_automodapi.automodapi', 'sphinx.ext.intersphinx']
# for external cross-references that are on the local machine, the relative
# path must be converted to an absolute path otherwise the intersphinx generates
# a broken link.
intersphinx_mapping = {'gdt-core':
                       ('https://astro-gdt.readthedocs.io/en/latest/', None)}

napoleon_google_docstring = True
napoleon_use_ivar = False
napoleon_use_rtype = False
numpydoc_show_class_members = False

# nbsphinx: the notebooks under notebooks/ are committed with their outputs
# already executed against the live HEASARC archive (see notebooks.rst), so
# never re-run them at doc-build time -- a build machine has no guarantee of
# network access to the archive, and re-running would defeat the point of
# committing real, verified outputs.
nbsphinx_execute = 'never'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
mathjax_path = "https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = '.rst'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'bootstrap-astropy'
html_theme_options = {
    'logotext1': 'gdt',  # white,  semi-bold
    'logotext2': '-burstcube',  # orange, light
    'logotext3': ':docs',   # white,  light
    'astropy_project_menubar': False
    }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- PDF options -------------------------------------------------------------
latex_documents = [('index', 'gammaray-data-tools.tex', project, author, 'manual')]
