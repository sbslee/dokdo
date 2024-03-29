# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))

# -- Project information -----------------------------------------------------

project = 'dokdo'
copyright = '2020, Seung-been "Steven" Lee'
author = 'Seung-been "Steven" Lee'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.linkcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.autosectionlabel',
    'sphinx_rtd_theme',
    'sphinx_issues'
]

autodoc_mock_imports = [
    'numpy',
    'pandas',
    'qiime2',
    'scipy',
    'matplotlib',
    'seaborn',
    'skbio',
    'statsmodels',
]

issues_github_path = 'sbslee/dokdo'

napoleon_use_param = False

# Include the example source for plots in API docs.
plot_include_source = True
plot_formats = [('png', 90)]
plot_html_show_formats = False
plot_html_show_source_link = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autosectionlabel_prefix_document = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# -- Add external links to source code with sphinx.ext.linkcode --------------

import inspect

def linkcode_resolve(domain, info):
    if domain != 'py':
        return None

    modname = info['module']

    if not modname:
        return None

    submod = sys.modules.get(modname)

    if submod is None:
        return None

    fullname = info['fullname']
    obj = submod

    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except AttributeError:
            return None

    try:
        fn = inspect.getsourcefile(inspect.unwrap(obj))
    except TypeError:
        fn = None
    if not fn:
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except OSError:
        lineno = None

    if lineno:
        linespec = f'#L{lineno}-L{lineno + len(source) - 1}'
    else:
        linespec = ''

    fn = fn.split('/api/')[1]

    return f'https://github.com/sbslee/dokdo/tree/master/dokdo/api/{fn}/{linespec}'
