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
from datetime import (
    date,
)

# from deepmd.common import (
#     ACTIVATION_FN_DICT,
#     PRECISION_DICT,
# )
# from deepmd.utils.argcheck import (
#     list_to_doc,
# )

sys.path.append(os.path.dirname(__file__))
# import sphinx_contrib_exhale_multiproject  # noqa: F401

# -- Project information -----------------------------------------------------


def run_apidoc(_):
    import sys

    from sphinx.ext.apidoc import (
        main,
    )

    sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    module = os.path.join(cur_dir, "..")
    main(
        [
            "-M",
            "--tocfile",
            "api_py",
            "-H",
            "Python API",
            "-o",
            os.path.join(cur_dir, "api_py"),
            module,
            "source/*",
            "--force",
        ]
    )


def setup(app):
    # Add hook for building doxygen xml when needed
    app.connect("builder-inited", run_apidoc)


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = [
#     'recommonmark',
#     "sphinx_rtd_theme",
#     'myst_parser',
#     'sphinx_markdown_tables',
#     'sphinx.ext.autosummary'
# ]

extensions = [
    # "deepmodeling_sphinx",
    # "dargs.sphinx",
    "sphinx_rtd_theme",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    # "sphinxarg.ext",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "autoapi.extension",
    # "sphinxcontrib.bibtex",
    # "myst_nb",
    # "numpydoc",
    # "breathe",
    "exhale",
]

# breathe_domain_by_extension = {
#         "h" : "cpp",
# }
breathe_projects = {
    "cc": "_build/cc/xml/",
    "c": "_build/c/xml/",
    "core": "_build/core/xml/",
}
breathe_default_project = "cc"

# exhale_args = {
#     "doxygenStripFromPath": "..",
#     # Suggested optional arguments
#     # "createTreeView":        True,
#     # TIP: if using the sphinx-bootstrap-theme, you need
#     # "treeViewIsBootstrap": True,
#     "exhaleExecutesDoxygen": True,
#     # "unabridgedOrphanKinds": {"namespace"}
#     # "listingExclude": [r"namespace_*"]
# }
# exhale_projects_args = {
#     "cc": {
#         "containmentFolder": "./API_CC",
#         "exhaleDoxygenStdin": "INPUT = ../source/api_cc/include/",
#         "rootFileTitle": "C++ API",
#         "rootFileName": "api_cc.rst",
#     },
#     "c": {
#         "containmentFolder": "./api_c",
#         "exhaleDoxygenStdin": "INPUT = ../source/api_c/include/",
#         "rootFileTitle": "C API",
#         "rootFileName": "api_c.rst",
#     },
#     "core": {
#         "containmentFolder": "./api_core",
#         "exhaleDoxygenStdin": """INPUT = ../source/lib/include/
#                                  PREDEFINED += GOOGLE_CUDA
#                                               TENSORFLOW_USE_ROCM
#         """,
#         "rootFileTitle": "Core API",
#         "rootFileName": "api_core.rst",
#     },
# }

# Tell sphinx what the primary language being documented is.
# primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
# highlight_language = 'cpp'

#
myst_heading_anchors = 4
nb_execution_mode = "off"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

intersphinx_mapping = {
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "python": ("https://docs.python.org/", None),
    "tensorflow": (
        "https://www.tensorflow.org/api_docs/python",
        "https://github.com/mr-ubik/tensorflow-intersphinx/raw/master/tf2_py_objects.inv",
    ),
    "ase": ("https://wiki.fysik.dtu.dk/ase/", None),
}
numpydoc_xref_param_type = True


numpydoc_xref_aliases = {}
import typing

for typing_type in typing.__all__:
    numpydoc_xref_aliases[typing_type] = "typing.%s" % typing_type

# rst_epilog = f"""
# .. |ACTIVATION_FN| replace:: {list_to_doc(ACTIVATION_FN_DICT.keys())}
# .. |PRECISION| replace:: {list_to_doc(PRECISION_DICT.keys())}
# """

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_logo = "_static/logo.svg"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

autodoc_default_flags = ["members"]
autosummary_generate = True
master_doc = "index"
mathjax_path = (
    "https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.0/es5/tex-mml-chtml.min.js"
)
myst_enable_extensions = [
    "dollarmath",
    "colon_fence",
]
myst_fence_as_directive = ("math",)
# fix emoji issue in pdf
latex_engine = "xelatex"
latex_elements = {
    "fontpkg": r"""
\usepackage{fontspec}
\setmainfont{Symbola}
""",
    "preamble": r"""
\usepackage{enumitem}
\setlistdepth{99}
""",
}

# For TF automatic generated OP docs
napoleon_google_docstring = True
napoleon_numpy_docstring = False

bibtex_bibfiles = ["../CITATIONS.bib"]





on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

here = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(here, '..')))

# -- Project information -----------------------------------------------------

project = u"PyWorld"
copyright = u"2024, Xiaobo Sun"
author = u"Xiaobo Sun"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = "0.1.0"
# The full version, including alpha/beta/rc tags.
release = version
# needs_sphinx = '2.1'
# -- General configuration ------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
# extensions = [
#     "sphinx.ext.autodoc",
#     "sphinx.ext.coverage",
#     "sphinx.ext.doctest",
#     "sphinx.ext.intersphinx",
#     "sphinx.ext.mathjax",
#     "sphinx.ext.napoleon",
#     "sphinx.ext.todo",
#     "sphinx.ext.viewcode",
#     "autoapi.extension",
# ]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
# exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

source_suffix = ['.rst', '.md']

# master_doc = 'index'

html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}



# -- Use autoapi.extension to run sphinx-apidoc -------

autoapi_dirs = ['../source']

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
