# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import os
import sys
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

here = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(here, '..')))

# -- Project information -----------------------------------------------------

project = u"MOMAP"
copyright = u"2021, Xiaobo Sun"
author = u"Xiaobo Sun"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = "0.1.0"
# The full version, including alpha/beta/rc tags.
release = version
needs_sphinx = '2.1'
# -- General configuration ------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "autoapi.extension",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

source_suffix = ['.rst', '.md']

master_doc = 'index'

html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}



# -- Use autoapi.extension to run sphinx-apidoc -------

autoapi_dirs = ['../SPEC','../TRANSPORT','../source']

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# -- Options for Intersphinx

intersphinx_mapping = {'python': ('https://docs.python.org/3', None),
                       # Commonly used libraries, uncomment when used in package
                       # 'numpy': ('http://docs.scipy.org/doc/numpy/', None),
                       # 'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None),
                       # 'scikit-learn': ('https://scikit-learn.org/stable/', None),
                       # 'matplotlib': ('https://matplotlib.org/stable/', None),
                       # 'pandas': ('http://pandas.pydata.org/docs/', None),
                       }
