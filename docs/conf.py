#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file was created by sphinx-quickstart on Sun Nov  6 17:05:42 2016
# and has been modified.
#
# Since Sherpa is tricky to build, in that it needs both a C and FORTRAN
# compiler, which has been hard to set up on Read The Docs (although it
# is getting close as of May 2018), the current approach is to mock
# all the modules needed to create the documentation. That is, there is
# no need to build Sherpa in order to build the documentation. This means
# that we can get documentation for optional modules, but it does
# require that the docstrings for functions in compiled modules are
# added to the Python module rather than within the C API (as an example,
# look at sherpa.utils which provides docstrings for the code in
# sherpa.utils._utils.
#
# The minimum supported version of Sphinx is 1.8.
#
# The minimum requirements are:
#    numpy  - since setup.py enforces this
#    sphinx_rtd_theme
#    sphinx_astropy
#
# The documentation can be built
#   a) from the top level with 'python setup.py build_sphinx'
#   b) from the docs/ directory with 'sphinx-build -b html . build/html'
#
import datetime
import glob
import os
import shutil
import sys

from sphinx_astropy.conf.v1 import intersphinx_mapping, default_role

import sphinx_rtd_theme


# Based on http://read-the-docs.readthedocs.io/en/latest/faq.html#i-get-import-errors-on-libraries-that-depend-on-c-modules
#
from unittest.mock import MagicMock as BaseMock


on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

runtime_dir = os.path.split(os.getcwd())[1]


# I found this somewhere (probably the rtd link given above). I was
# hoping it would support building with Python 2.7 but it doesn't seem
# to. I have left it in for now, even though Python 2.7 is no longer
# supported, in case it is being used somewhere.
#
class Mock(BaseMock):
    @classmethod
    def __getattr__(cls, name):
        # Hope that none of the attributes need a recursive Mock
        return BaseMock()


class PylabMock(Mock):
    """Provide defaults for pybab.errorbar"""

    # The pylab backend inspects the errorbar function to extract the
    # ecolor, capsize, and barsabove settings. So this will need to
    # be kept up-to-date with pylab.errorbar.
    #
    def errorbar(self, ecolor=None, capsize=None, barsabove=False):
        pass


class XSPECMock(Mock):
    """Need to override get_xsversion to return a "useful" string"""

    def get_xsversion(self):
        "Presumably XSPEC will never reach this"
        return "999.999.999"


# It's not clear if pylab is needed here. At present non of the
# CIAO-related modules (e.g. the crates I/O backend) are mocked;
# should they be?
#
sys.modules['astropy.io.fits'] = Mock()
sys.modules['pylab'] = PylabMock()

# This list was found by trying to import various Sherpa modules, adding
# in mocks as necessary, and repeating. Some of these may be unnescessary.
#
for mod_name in ['sherpa.utils._utils',
                 'sherpa.utils._psf',
                 'sherpa.models._modelfcts',
                 'sherpa.estmethods._est_funcs',
                 'sherpa.optmethods._saoopt',
                 'sherpa.stats._statfcts',
                 'sherpa.astro.utils._utils',
                 'sherpa.astro.utils._pileup',
                 'sherpa.astro.utils._region',
                 'sherpa.astro.utils._wcs',
                 'sherpa.astro.models._modelfcts',
                 'sherpa.astro.io.pyfits_backend',
                 'sherpa.image.DS9',
                 # Extra modules
                 'group'
                ]:
    sys.modules[mod_name] = Mock()

# Specialized mocks
#
sys.modules['sherpa.astro.xspec._xspec'] = XSPECMock()

# If in the docs directory then add the parent to the system path
# so that Sherpa can be found.
#
if runtime_dir == 'docs':
    sys.path.insert(0, os.path.abspath('..'))


import sherpa

# For now include the '+...' part of the version string
# in the full version, but drop the ".dirty" suffix.
#
sherpa_release = sherpa._version.get_versions()['version']
if on_rtd and sherpa_release.endswith('.dirty'):
    sherpa_release = sherpa_release[:-6]

sherpa_version = sherpa_release[:sherpa_release.find('+')]

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# The use of app.add_js_file requires Sphinx 1.8.
#
needs_sphinx = '1.8'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    # Use napoleon over numpydoc for now since it stops a large number
    # of warning messages (about missing links) that I don't have time
    # to investigate.
    'sphinx.ext.napoleon',
    # 'numpydoc.numpydoc',
    'sphinx.ext.intersphinx',
    'sphinx_astropy.ext.intersphinx_toggle',
    'sphinx_astropy.ext.edit_on_github',
    # notebooks
    'nbsphinx',
    'sphinx_rtd_theme',
    'pytest_doctestplus.sphinx.doctestplus',
]
# Imported from sphinx_astropy so we don't have to maintain the list
# of servers
# intersphinx_mapping = {'python': ('https://docs.python.org/3', None)}

# notebook support
# - for now never execute a notebook
#
nbsphinx_execute = 'never'

# Copy over the notebooks to the documentation area
# if they don't exist or are older. This assumes we
# are running in the top-level directory!
#
nbmapping = {}
infile = 'notebooks/nbmapping.dat'
if runtime_dir == 'docs':
    infile = f"../{infile}"

try:
    with open(infile, 'r') as fh:
        for l in fh.readlines():
            l = l.strip()
            if l == '' or l.startswith('#'):
                continue

            toks = l.split('\t')
            if len(toks) != 2:
                sys.stderr.write(f"Invalid nbmapping line: {l}\n")
                continue

            nbmapping[toks[0]] = toks[1]
except OSError:
    pass


def getmtime(infile):

    try:
        t = os.path.getmtime(infile)
    except FileNotFoundError:
        return None
    return datetime.datetime.fromtimestamp(t)


if len(nbmapping) == 0:
    sys.stderr.write("No notebook mapping found\n")

else:
    path = 'notebooks/*ipynb'
    if runtime_dir == 'docs':
        path = f"../{path}"

    nbfiles = glob.glob(path)
    if len(nbfiles) == 0:
        sys.stderr.write("No notebooks found!\n")

    for nbfile in nbfiles:
        head = os.path.basename(nbfile)
        try:
            outpath = nbmapping[head]
        except KeyError:
            sys.stderr.write(f"Dropping nbfile: {nbfile}\n")
            continue

        if runtime_dir == 'docs':
            assert outpath.startswith('docs/')
            outpath = outpath[5:]

        t1 = getmtime(nbfile)
        t2 = getmtime(outpath)

        if t2 is None or t2 < t1:
            shutil.copy2(nbfile, outpath)


# prolog/epilog based on
# https://nbsphinx.readthedocs.io/en/0.7.1/conf.py
#
# Note: we assume the notebooks are in the notebooks/
# top-level directory and are copied into the docs/
# directory at build time.
#
nbsphinx_prolog = r"""
{% set docname = 'notebooks/' + env.docname.split('/')[-1] %}

.. raw:: html

    <div class="admonition note">

      This page was generated from the Jupyter notebook
      <a class="reference external" href="https://github.com/sherpa/sherpa/blob/{{ env.config.githash|e }}/{{ docname|e }}.ipynb">{{ docname|e }}</a>.
      <script>
        if (document.location.host) {
          $(document.currentScript).replaceWith(
            '<a class="reference external" ' +
            'href="https://nbviewer.jupyter.org/url' +
            (window.location.protocol == 'https:' ? 's/' : '/') +
            window.location.host +
            window.location.pathname.slice(0, -4) +
            'ipynb">View in <em>nbviewer</em></a>.'
          );
        }
      </script>
    </div>

"""

# docstrings
#
napoleon_google_docstring = False

autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = '.rst'

# The encoding of source files.
#
# source_encoding = 'utf-8-sig'

# The main toctree document.
main_doc = 'index'

# General information about the project.
project = 'Sherpa'
copyright = '2019-2022, Chandra X-ray Center, Smithsonian Astrophysical Observatory.'
author = 'Chandra X-ray Center, Smithsonian Astrophysical Observatory'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = sherpa_version
# The full version, including alpha/beta/rc tags.
release = sherpa_release

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#
# today = ''
#
# Else, today_fmt is used as the format for a strftime call.
#
# today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The reST default role (used for this markup: `text`) to use for all
# documents.
# Imported from sphinx_astropy
# default_role = 'obj'

# Setting from sphinx_astropy for numpydoc xref settings
numpydoc_xref_param_type = True

# If true, '()' will be appended to :func: etc. cross-reference text.
#
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# Graphviz based on values from AstroPy - see
# https://github.com/astropy/sphinx-astropy/blob/main/sphinx_astropy/conf/v1.py
#
graphviz_output_format = "svg"

graphviz_dot_args = [
    '-Nfontsize=10',
    '-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif',
    '-Efontsize=10',
    '-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif',
    '-Gfontsize=10',
    '-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif'
]

# Ensure sphinx_astropy.ext.edit_on_github knows where to send
# the edit links.
#
edit_on_github_project = 'sherpa/sherpa'
# edit_on_github_branch = '4.10.1'
edit_on_github_source_root = ''
edit_on_github_doc_root = 'docs'

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# html_theme_options = {}
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# The name for this set of Sphinx documents.
# "<project> v<release> documentation" by default.
#
# html_title = 'Sherpa v4.8.2'

# A shorter title for the navigation bar.  Default is the same as html_title.
#
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#
# html_logo = None

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
# html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
#
# html_extra_path = []

# If not None, a 'Last updated on:' timestamp is inserted at every page
# bottom, using the given strftime format.
# The empty string is equivalent to '%b %d, %Y'.
#
# html_last_updated_fmt = None

# Follow AstroPy's convention for the date format.
#
html_last_updated_fmt = '%d %b %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#
# html_additional_pages = {}

# If false, no module index is generated.
#
# html_domain_indices = True

# If false, no index is generated.
#
# html_use_index = True

# If true, the index is split into individual pages for each letter.
#
# html_split_index = False

# If true, links to the reST sources are added to the pages.
#
# html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#
# html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#
# html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Language to be used for generating the HTML full-text search index.
# Sphinx supports the following languages:
#   'da', 'de', 'en', 'es', 'fi', 'fr', 'h', 'it', 'ja'
#   'nl', 'no', 'pt', 'ro', 'r', 'sv', 'tr', 'zh'
#
# html_search_language = 'en'

# A dictionary with options for the search language support, empty by default.
# 'ja' uses this config value.
# 'zh' user can custom change `jieba` dictionary path.
#
# html_search_options = {'type': 'default'}

# The name of a javascript file (relative to the configuration directory) that
# implements a search results scorer. If empty, the default will be used.
#
# html_search_scorer = 'scorer.js'

# Output file base name for HTML help builder.
htmlhelp_basename = 'Sherpadoc'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
     # The paper size ('letterpaper' or 'a4paper').
     #
     # 'papersize': 'letterpaper',

     # The font size ('10pt', '11pt' or '12pt').
     #
     # 'pointsize': '10pt',

     # Additional stuff for the LaTeX preamble.
     #
     # 'preamble': '',

     # Latex figure (float) alignment
     #
     # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (main_doc, 'Sherpa.tex', 'Sherpa Documentation',
     'Chandra X-ray Center', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#
# latex_use_parts = False

# If true, show page references after internal links.
#
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
#
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
#
# latex_appendices = []

# It false, will not define \strong, \code, 	itleref, \crossref ... but only
# \sphinxstrong, ..., \sphinxtitleref, ... To help avoid clash with user added
# packages.
#
# latex_keep_old_macro_names = True

# If false, no module index is generated.
#
# latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (main_doc, 'sherpa', 'Sherpa Documentation',
     [author], 1)
]

# If true, show URL addresses after external links.
#
# man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (main_doc, 'Sherpa', 'Sherpa Documentation',
     author, 'Sherpa',
     'Sherpa is a modeling and fitting environment for Python.',
     'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#
# texinfo_appendices = []

# If false, no module index is generated.
#
# texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#
# texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#
# texinfo_no_detailmenu = False

# -- try and get copybuttons to work --
#
# docs/_static/copybutton.js is copied from
# https://raw.githubusercontent.com/scipy/scipy-sphinx-theme/master/_theme/scipy/static/js/copybutton.js
# version is from
# https://github.com/scipy/scipy-sphinx-theme/commit/a8aa8a6aad1524c9577a861fc4faa82d6c167138
#

def setup(app):
    app.add_js_file('copybutton.js')

    # Add in the git commit id so that the nbsphinx_prolog
    # can pick it up. This may well be repeating logic that
    # is already known by Sphinx.
    #
    # Guessing at rebuild and types arguments.
    #
    githash = sherpa._version.get_versions()['full-revisionid']
    if githash.endswith('.dirty'):
        githash = githash[:-6]

    app.config.add('githash', githash, rebuild=False, types=None)
