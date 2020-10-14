#
#  Copyright (C) 2007, 2014, 2015, 2016, 2019, 2020
#     Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""

Modeling and fitting package for scientific data analysis

Sherpa is a modeling and fitting package for scientific data analysis.
It includes general tools of interest to all users as well as
specialized components for particular disciplines (e.g. astronomy).

Note that the top level sherpa package does not import any
subpackages.  This allows the user to import a particular component
(e.g. `sherpa.optmethods`) without pulling in any others.  To import all
the standard subpackages, use ``import sherpa.all`` or
``from sherpa.all import *``.

"""

import datetime
import logging
import os
import os.path
import subprocess
import sys


__all__ = ('citation', 'get_config', 'get_include', 'smoke')

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


class Formatter(logging.Formatter):
    def format(self, record):
        if record.levelno > logging.INFO:
            msg = '%s: %s' % (record.levelname, record.msg)
        else:
            msg = record.msg
        return msg

log = logging.getLogger('sherpa')
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(Formatter())
log.addHandler(handler)
log.setLevel(logging.INFO)

dbg = log.debug

del Formatter, log, handler



DEFAULT_CITATION = """Please review the Zenodo Sherpa page at
       https://doi.org/10.5281/zenodo.593753
to identify the latest release of Sherpa. The Zenodo page
https://help.zenodo.org/#versioning provides information on how to
cite a specific version.

If you want a general-purpose citation then please use either of the
following, kindly provided by the SAO/NASA Astrophysics Data System
service:

@INPROCEEDINGS{2001SPIE.4477...76F,
       author = {{Freeman}, Peter and {Doe}, Stephen and {Siemiginowska}, Aneta},
        title = "{Sherpa: a mission-independent data analysis application}",
     keywords = {Astrophysics},
    booktitle = {Astronomical Data Analysis},
         year = 2001,
       editor = {{Starck}, Jean-Luc and {Murtagh}, Fionn D.},
       series = {Society of Photo-Optical Instrumentation Engineers (SPIE) Conference Series},
       volume = {4477},
        month = nov,
        pages = {76-87},
          doi = {10.1117/12.447161},
archivePrefix = {arXiv},
       eprint = {astro-ph/0108426},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2001SPIE.4477...76F},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@INPROCEEDINGS{2007ASPC..376..543D,
   author = {{Doe}, S. and {Nguyen}, D. and {Stawarz}, C. and {Refsdal}, B. and
	{Siemiginowska}, A. and {Burke}, D. and {Evans}, I. and {Evans}, J. and
	{McDowell}, J. and {Houck}, J. and {Nowak}, M.},
    title = "{Developing Sherpa with Python}",
booktitle = {Astronomical Data Analysis Software and Systems XVI},
     year = 2007,
   series = {Astronomical Society of the Pacific Conference Series},
   volume = 376,
   editor = {{Shaw}, R.~A. and {Hill}, F. and {Bell}, D.~J.},
    month = oct,
    pages = {543},
   adsurl = {http://adsabs.harvard.edu/abs/2007ASPC..376..543D},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""


def _make_citation(version, title, date, authors, idval, latest=True):

    year = date.year
    month = date.strftime('%b').lower()
    nicedate = date.strftime('%B %d, %Y')

    authnames = ' and\n                  '.join(authors)

    if latest:
        out = "The latest release of Sherpa is {}\n".format(title)
        out += "released on {}.".format(nicedate)
    else:
        out = "Sherpa {} was released on {}.".format(version, nicedate)

    out += """

@software{{sherpa_{year}_{idval},
  author       = {{{authnames}}},
  title        = {{{title}}},
  month        = {month},
  year         = {year},
  publisher    = {{Zenodo}},
  version      = {{{version}}},
  doi          = {{10.5281/zenodo.{idval}}},
  url          = {{https://doi.org/10.5281/zenodo.{idval}}}
}}

""".format(authnames=authnames, idval=idval, title=title,
           version=version, year=year, month=month)

    out += DEFAULT_CITATION
    return out


def _get_citation_hardcoded(version):
    """Retrieve the citation information.

    Parameters
    ----------
    version : str
        The version to retrieve the citation before. It is expected
        to be '4.8.0' or higher.

    Returns
    -------
    citation : str or None
        Citation information if known, otherwise None.

    """

    def todate(year, mnum, dnum):
        return datetime.datetime(year, mnum, dnum)

    db = 'DougBurke'
    ol = 'Omar Laurino'
    dn = 'dtnguyen2'
    ta = 'Tom Aldcroft'
    an = 'Aneta Siemiginowska'

    jb = 'Jamie Budynkiewicz'
    cd = 'Christoph Deil'
    bs = 'Brigitta Sipocz'

    wm = 'wmclaugh'

    kl = 'Katrin Leinweber'

    mt = 'Marie-Terrell'
    bs2 = 'Brigitta Sipőcz'
    hm = 'Hans Moritz Günther'
    td = 'Todd'

    cite = {}

    def add(version, **kwargs):
        assert version not in cite
        cite[version] = dict(**kwargs)
        cite[version]['version'] = version

    add(version='4.8.0', title='sherpa: Sherpa 4.8.0',
        date=todate(2016, 1, 27),
        authors=[db, ol, dn, ta, an],
        idval='45243')
    add(version='4.8.1', title='sherpa: Sherpa 4.8.1',
        date=todate(2016, 4, 15),
        authors=[db, ol, dn, ta, jb, an, cd, bs],
        idval='49832')
    add(version='4.8.2', title='sherpa: Sherpa 4.8.2',
        date=todate(2016, 9, 23),
        authors=[db, ol, dn, jb, ta, an, cd, wm, bs],
        idval='154744')

    add(version='4.9.0', title='sherpa/sherpa: Sherpa 4.9.0',
        date=todate(2017, 1, 27),
        authors=[db, ol, dn, jb, ta, an, wm, cd, bs],
        idval='260416')
    add(version='4.9.1', title='sherpa/sherpa: Sherpa 4.9.1',
        date=todate(2017, 8, 3),
        authors=[db, ol, dn, jb, ta, an, wm, cd, bs],
        idval='838686')

    add(version='4.10.0', title='sherpa/sherpa: Sherpa 4.10.0',
        date=todate(2018, 5, 1),
        authors=[db, ol, dn, jb, ta, an, wm, cd, bs],
        idval='1245678')
    add(version='4.10.1', title='sherpa/sherpa: Sherpa 4.10.1',
        date=todate(2018, 10, 16),
        authors=[db, ol, dn, jb, ta, an, wm, bs, cd, kl],
        idval='1463962')
    add(version='4.10.2', title='sherpa/sherpa: Sherpa 4.10.2',
        date=todate(2018, 12, 14),
        authors=[db, ol, dn, jb, ta, an, cd, wm, bs, kl],
        idval='2275738')

    add(version='4.11.0', title='sherpa/sherpa: Sherpa 4.11.0',
        date=todate(2019, 2, 20),
        authors=[db, ol, dn, jb, ta, an, cd, wm, bs, kl],
        idval='2573885')
    add(version='4.11.1', title='sherpa/sherpa: Sherpa 4.11.1',
        date=todate(2019, 8, 1),
        authors=[db, ol, dn, jb, ta, an, cd, wm, bs, kl],
        idval='3358134')

    add(version='4.12.0', title='sherpa/sherpa: Sherpa 4.12.0',
        date=todate(2020, 1, 30),
        authors=[db, ol, dn, wm, jb, ta, an, cd, mt, bs2, hm, td, kl],
        idval='3631574')
    add(version='4.12.1', title='sherpa/sherpa: Sherpa 4.12.1',
        date=todate(2020, 7, 14),
        authors=[db, ol, wm, dn, mt, hm, jb, an, ta, cd, bs2, kl, td],
        idval='3944985')

    kwargs = cite.get(version, None)
    if kwargs is None:
        return None

    return _make_citation(**kwargs, latest=False)


def _download_json(url):
    """Return the JSON data or None.

    Parameters
    ----------
    url : str
        The URL to process.

    Returns
    -------
    response : dict
        The data in JSON as 'success' or an error message as 'failed'.

    """

    import json
    import ssl
    from urllib.request import Request, urlopen
    from urllib.error import HTTPError

    # Need to set the Content-Type so it looks like I need
    # a request. Darn it. Even though https://developers.zenodo.org/#list36
    # claims you can set the content-type, it just seems to
    # return JSON.
    #
    # req = Request(url,
    #               headers={'Content-Type': 'application/x-bibtex'})
    #
    req = Request(url)

    # This is not ideal, but given the way that we package
    # CIAO we can end up with issues evaluating SSL certificates,
    # so as we are not sending secure information we drop
    # the SSL checks.
    #
    context = ssl._create_unverified_context()

    try:
        res = urlopen(req, context=context)
    except HTTPError as he:
        dbg("Unable to access {}: {}".format(url, he))
        return {'failed': 'Unable to access the Zenodo site.'}

    try:
        jsdata = json.load(res)
    except UnicodeDecodeError as ue:
        dbg("Unable to decode JSON from Zenodo: {}".format(ue))
        return {'failed': 'Unable to understand the response from Zenodo.'}

    return {'success': jsdata}


def _make_zenodo_citation(jsdata, latest=True):
    """Convert Zenodo record to a citation.

    Parameters
    ----------
    jsdata : object
        A Zenodo record.

    Returns
    -------
    response : dict
        Success in 'success' (dict) or failure message in 'failed'.

    """

    try:
        created = jsdata['created']
        mdata = jsdata['metadata']
        idval = jsdata['id']
    except KeyError as ke:
        dbg("Unable to find metadata: {}".format(ke))
        return {'failed': 'Unable to parse the Zenodo response.'}

    created = created.split('T')[0]
    isoformat = '%Y-%m-%d'
    try:
        date = datetime.datetime.strptime(created, isoformat)
    except ValueError:
        dbg("Unable to convert created: '{}'".format(created))
        return {'failed': 'Unable to parse the Zenodo response.'}

    try:
        version = mdata['version']
        title = mdata['title']
        creators = mdata['creators']
    except KeyError as ke:
        dbg("Unable to find metadata: {}".format(ke))
        return {'failed': 'Unable to parse the Zenodo response.'}

    authors = [c['name'] for c in creators]
    out = _make_citation(version=version, title=title, date=date,
                         authors=authors, idval=idval,
                         latest=latest)
    return {'success': out}


def _get_citation_version():
    """What version of Sherpa are we using?

    Returns
    -------
    version : str
        The Sherpa version.
    """

    out = 'You are using Sherpa {}'.format(__version__)
    if '+' in __version__:
        out += " (it is not a released version)"

    out += ".\n\n"
    return out


def _get_citation_zenodo_failure(failed):
    """Standard response when there's a problem.

    Returns
    -------
    text : str
        The failure message.

    """
    out = 'There was a problem retireving the data from Zenodo:\n'
    out += failed
    out += '\n\n'
    out += DEFAULT_CITATION
    return out


def _get_citation_zenodo_latest():
    """Query Zenodo for the latest release.

    Returns
    -------
    citation : str
        Citation information. It will include information
        on any failure.
    """

    # Can we retrieve the information from Zenodo?
    #
    # We do not access the DOI but the Zenodo API
    # url = 'https://doi.org/10.5281/zenodo.593753'
    url = 'https://zenodo.org/api/records/593753'
    jsdata = _download_json(url)
    if 'failed' in jsdata:
        return _get_citation_zenodo_failure(jsdata['failed'])

    out = _make_zenodo_citation(jsdata['success'])
    if 'failed' in out:
        return _get_citation_zenodo_failure(out['failed'])

    return out['success']


def _parse_zenodo_data(jsdata, version):
    """Extract data for the given version from a Zenodo query.

    Parameters
    ----------
    jsdata : json
        The response from Zenodo.
    version : str
        The version we are looking for.

    Returns
    -------
    response : dict
        Success in 'success' (dict) or failure message in 'failed'.
    """

    try:
        hits = jsdata['hits']['hits']
    except KeyError:
        dbg("Unable to find hits/hits")
        return {'failed': 'Unable to parse the Zenodo response.'}

    data = None
    try:
        for hit in hits:
            if hit['metadata']['version'] == version:
                data = hit
                break

    except KeyError:
        dbg("Record missing version")
        return {'failed': 'Unable to parse the Zenodo response.'}

    except TypeError:
        dbg('hits/hits is not iterable!')
        return {'failed': 'Unable to parse the Zenodo response.'}

    if data is None:
        dbg('Version {} not found'.format(version))
        return {'failed': "Zenodo has no information for version {}.".format(version)}

    return {'success': data}


def _get_citation_zenodo_version(version):
    """Query Zenodo for the specific release.

    As this has to return all Sherpa records it is slow.

    Returns
    -------
    citation : str
        Citation information.
    """

    # Is there a better way to do this?
    #
    url = 'https://zenodo.org/api/records/?q=conceptrecid:"593753"&all_versions=True'
    jsdata = _download_json(url)
    if 'failed' in jsdata:
        return _get_citation_zenodo_failure(jsdata['failed'])

    data = _parse_zenodo_data(jsdata['success'], version)
    if 'failed' in data:
        return _get_citation_zenodo_failure(data['failed'])

    out = _make_zenodo_citation(jsdata['success'], latest=False)
    if 'failed' in out:
        return _get_citation_zenodo_failure(out['failed'])

    return out['success']


def _get_citation(version='current'):
    """Retrieve the citation information.

    Parameters
    ----------
    version : str, optional
        The version to retrieve the citation for. The supported values
        are limited to 'current', to return the citation for the
        installed version of Sherpa, 'latest' which will return the
        latest release, and the current set of releases available on
        Zenodo (this goes back to '4.8.0').

    Returns
    -------
    citation : str
        Citation information.

    """

    vstr = _get_citation_version()

    if version == 'latest':
        return vstr + _get_citation_zenodo_latest()

    if version == 'current':
        # Replace with the version (excluding any extraneous
        # information). Note that development versions close to a
        # release will have a version which has no release, but
        # I believe it is okay to say "Hey, Zenodo has no info
        # on this", rather than to try and fall back to the
        # previous version (since there's no good Zenodo reference
        # in this case).
        #
        version = __version__.split('+')[0]

    # We could include a check on the version number (e.g. a.b.c
    # where all elements are an integer), but do we want to
    # require this naming scheme?

    # If we know this version there's no need to call Zenodo
    #
    out = _get_citation_hardcoded(version)
    if out is not None:
        return vstr + out

    # In case the hardcoded list hasn't been updated.
    #
    return vstr + _get_citation_zenodo_version(version)


def citation(version='current', filename=None, clobber=False):
    """Return citatation information for Sherpa.

    The citation information is taken from Zenodo [1]_, using the
    Sherpa "latest release" identifier [2]_, and so requires an
    internet connection. The message is displayed on screen, using
    pagination, or can be written to a file.

    Parameters
    ----------
    version : str, optional
        The version to retrieve the citation for. The supported values
        are limited to 'current', to return the citation for the
        installed version of Sherpa, 'latest' which will return the
        latest release, and the current set of releases available on
        Zenodo (this goes back to '4.8.0').
    filename : str or StringIO or None, optional
        If not None, write the output to the given file or filelike
        object.
    clobber : bool, optional
        If filename is a string, then - when clobber is set - refuse
        to overwrite the file if it already exists.

    Notes
    -----
    If there is no internet connection, or there was a problem in
    downloading the data, or the Zenodo API has started to return
    different informatioon than expected, then the code will return
    information on why the call failed and other citation options.

    Zenodo only lets you perform a limited number of calls in
    a set time span, so if you call this routine too many times
    then it may start to fail.

    If a specific version is given then a hard-coded list of versions
    is checked, and if it matches then this information is used,
    rather than requiring a call to Zenodo.

    References
    ----------

    .. [1] https://zenodo.org/

    .. [2] https://doi.org/10.5281/zenodo.593753

    Examples
    --------

    Display the citation information for the current release on
    Zenodo. The information is paged to the display:

    >>> import sherpa
    >>> sherpa.citation()

    Write out the citation information for Sherpa 4.12.1 to the file
    ``cite.txt``:

    >>> sherpa.citation('4.12.1', outfile='cite.txt')

    Display the information for the latest release:

    >>> sherpa.citation('latest')

    """

    from sherpa.utils import send_to_pager
    cite = _get_citation(version=version)
    send_to_pager(cite, filename=filename, clobber=clobber)


def get_include():
    "Get the root path for installed Sherpa header files"

    return os.path.join(os.path.dirname(__file__), 'include')


def get_config():
    "Get the path for the installed Sherpa configuration file"

    filename = "sherpa-standalone.rc"

    home_dir = None
    config = None

    # If NOSHERPARC is set, read in system config file
    # ignore any user config file
    if (('NOSHERPARC' in os.environ) == True):
        return os.path.join(os.path.dirname(__file__), filename)

    # If SHERPARC is set, read in config file from there,
    # and ignore default location
    if (('SHERPARC' in os.environ) == True):
        config = os.environ.get('SHERPARC')
        if os.path.isfile(config):
            return config

    # SHERPARC was not set, so look for .sherpa.rc in default
    # location, which is user's home directory.
    home_dir = os.environ.get('HOME')
    config = os.path.join(home_dir, '.'+filename)

    if os.path.isfile(config):
        return config

    # If no user config file is set, fall back to system config file
    return os.path.join(os.path.dirname(__file__), filename)


def smoke(verbosity=0, require_failure=False, fits=None, xspec=False, ds9=False):
    """Run Sherpa's "smoke" test.

    The smoke test is a simple test that ensures the Sherpa
    installation is functioning. It is not a complete test suite, but
    it fails if obvious issues are found.

    Parameters
    ----------
    verbosity : int, optional
        The level of verbosity of this test
    require_failure : boolean, optional
        For debugging purposes, the smoke test may be required to
        always fail. Defaults to False.
    fits : str or None, optional
        Require a fits module with this name to be present before
        running the smoke test. This option makes sure that when the
        smoke test is run the required modules are present.  Note that
        tests requiring fits may still run if any fits backend is
        available, and they might still fail on their own.
    xspec : boolean, optional
        Require xspec module when running tests. Tests requiring xspec
        may still run if the xspec module is present.
    ds9 : boolean, optional
        Requires DS9 when running tests.

    Raises
    ------
    SystemExit
        Raised if any errors are found during the tests.

    """
    from sherpa.astro.utils import smoke
    smoke.run(verbosity=verbosity, require_failure=require_failure, fits=fits, xspec=xspec, ds9=ds9)


def _smoke_cli(verbosity=0, require_failure=False, fits=None, xspec=False, ds9=False):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-v", "--verbosity", dest="verbosity",
                      help="verbosity level")
    parser.add_option("-0", "--require-failure", dest="require_failure", action="store_true",
                      help="require smoke test to fail (for debugging)")
    parser.add_option("-f", "--fits-module", dest="fits", action="store",
                      help="require a specific fits module to be present")
    parser.add_option("-x", "--require-xspec", dest="xspec", action="store_true",
                      help="require xspec module")
    parser.add_option("-d", "--require-ds9", dest="ds9", action="store_true",
                      help="require DS9")

    options, _ = parser.parse_args()

    xspec = options.xspec or xspec
    verbosity = options.verbosity or verbosity
    require_failure = options.require_failure or require_failure
    fits = options.fits or fits
    ds9 = options.ds9 or ds9

    smoke(verbosity=verbosity, require_failure=require_failure, fits=fits, xspec=xspec, ds9=ds9)


def _install_test_deps():
    def install(package_name):
        try:
            subprocess.call([sys.executable, '-m', 'pip', 'install', package_name],
                               stdout=sys.stdout, stderr=sys.stderr)
        except:
            print("""Cannot import pip or install packages with it.
            You need pytest, and possibly pytest-cov, in order to run the tests.
            If you downloaded the source code, please run 'pip install -r test_requirements.txt'
            from the source directory first.
            """)
            raise

    deps = ['pytest>=3.3,!=5.2.3']  # List of packages to be installed
    pytest_plugins = []  # List of pytest plugins to be installed

    # If the plugins are installed "now", pytest won't load them because they are not registered as python packages
    # by the time this code runs. So we need to keep track of the plugins and explicitly pass them to `pytest.main()`.
    # Note that explicitly passing plugins that were already registered would result in an error, hence this
    # complexity seems to be required.
    installed_plugins = []

    for dep in deps:
        try:
            __import__(dep)
        except ImportError:
            install(dep)

    for plugin_name in pytest_plugins:
        module = plugin_name.replace("-", "_")
        try:
            __import__(module)
        except ImportError:
            install(plugin_name)
            installed_plugins.append(module)

    return installed_plugins


def clitest():
    plugins = _install_test_deps()
    import pytest
    import os
    sherpa_dir = os.path.dirname(__file__)
    errno = pytest.main([sherpa_dir, '-rs'], plugins=plugins)  # passing the plugins that have been installed "now".
    sys.exit(errno)
