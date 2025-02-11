#
#  Copyright (C) 2007, 2014 - 2016, 2019 - 2024
#  Smithsonian Astrophysical Observatory
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
import importlib
import importlib.resources
import logging
import os
import os.path
import subprocess
import sys
from typing import Any, Optional

from . import _version
__version__ = _version.get_versions()['version']


__all__ = ('citation', 'get_config', 'get_include', 'smoke')


class Formatter(logging.Formatter):
    def format(self, record):
        # Get the full message, #1688 shows we can not rely
        # in record.msg.
        #
        msg = record.getMessage()
        if record.levelno > logging.INFO:
            msg = f"{record.levelname}: {msg}"
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

The most recent article will appear in the Astrophysical Journal Supplement Series and it is
available on the archives:

@ARTICLE{2024arXiv240910400S,
       author = {{Siemiginowska}, Aneta and {Burke}, Douglas and {G{\"u}nther}, Hans Moritz and {Lee}, Nicholas P. and {McLaughlin}, Warren and {Principe}, David A. and {Cheer}, Harlan and {Fruscione}, Antonella and {Laurino}, Omar and {McDowell}, Jonathan and {Terrell}, Marie},
        title = "{Sherpa: An Open Source Python Fitting Package}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2024,
        month = sep,
          eid = {arXiv:2409.10400},
        pages = {arXiv:2409.10400},
          doi = {10.48550/arXiv.2409.10400},
archivePrefix = {arXiv},
       eprint = {2409.10400},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv240910400S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


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


def _make_citation(version: str,
                   title: str,
                   date: datetime.datetime,
                   authors: list[str],
                   idval: str,
                   latest: bool = True) -> str:

    year = date.year
    month = date.strftime('%b').lower()
    nicedate = date.strftime('%B %d, %Y')

    authnames = ' and\n                  '.join(authors)

    if latest:
        out = f"The latest release of Sherpa is {title}\n"
        out += f"released on {nicedate}."
    else:
        out = f"Sherpa {version} was released on {nicedate}."

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


def _get_citation_hardcoded(version: str) -> Optional[str]:
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

    Notes
    -----
    The entries can be created with the script:
    scripts/make_zenodo_release.py

    """

    def todate(year, mnum, dnum):
        return datetime.datetime(year, mnum, dnum)

    cite = {}

    def add(version, **kwargs):
        assert version not in cite
        cite[version] = dict(**kwargs)
        cite[version]['version'] = version

    add(version='4.16.1', title='sherpa/sherpa: Sherpa 4.16.1',
        date=todate(2024, 5, 21),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'Hans Moritz Günther', 'Marie-Terrell', 'dtnguyen2', 'Aneta Siemiginowska', 'Harlan Cheer', 'Jamie Budynkiewicz', 'Tom Aldcroft', 'luzpaz', 'Christoph Deil', 'Brigitta Sipőcz', 'Johannes Buchner', 'nplee', 'Axel Donath', 'Iva Laginja', 'Katrin Leinweber', 'Todd'],
        idval='11236879')
    add(version='4.16.0', title='sherpa/sherpa: Sherpa 4.16.0',
        date=todate(2023, 10, 26),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'Hans Moritz Günther', 'Marie-Terrell', 'dtnguyen2', 'Aneta Siemiginowska', 'Harlan Cheer', 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Christoph Deil', 'Brigitta Sipőcz', 'Johannes Buchner', 'nplee', 'Axel Donath', 'Iva Laginja', 'Katrin Leinweber', 'Todd'],
        idval='825839')

    add(version='4.15.1', title='sherpa/sherpa: Sherpa 4.15.1',
        date=todate(2023, 5, 18),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'Marie-Terrell',
                 'dtnguyen2', 'Hans Moritz Günther', 'Aneta Siemiginowska',
                 'Jamie Budynkiewicz', 'Harlan Cheer', 'Tom Aldcroft',
                 'Christoph Deil', 'Brigitta Sipőcz', 'Johannes Buchner',
                 'Axel Donath', 'Iva Laginja', 'Katrin Leinweber',
                 'nplee', 'Todd'],
        idval='7948720')
    add(version='4.15.0', title='sherpa/sherpa: Sherpa 4.15.0',
        date=todate(2022, 10, 11),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'dtnguyen2',
                 'Hans Moritz Günther', 'Marie-Terrell', 'Aneta Siemiginowska',
                 'Jamie Budynkiewicz', 'Harlan Cheer', 'Tom Aldcroft',
                 'Christoph Deil', 'Brigitta Sipőcz', 'Johannes Buchner',
                 'Axel Donath', 'Iva Laginja', 'Katrin Leinweber',
                 'nplee', 'Todd'],
        idval='7186379')
    add(version='4.14.1', title='sherpa/sherpa: Sherpa 4.14.1',
        date=todate(2022, 5, 20),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'dtnguyen2',
                 'Hans Moritz Günther', 'Marie-Terrell', 'Aneta Siemiginowska',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Christoph Deil',
                 'Harlan Cheer', 'Brigitta Sipőcz', 'Johannes Buchner',
                 'Axel Donath', 'Iva Laginja', 'Katrin Leinweber',
                 'nplee', 'Todd'],
        idval='6567264')
    add(version='4.14.0', title='sherpa/sherpa: Sherpa 4.14.0',
        date=todate(2021, 10, 7),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'dtnguyen2',
                 'Hans Moritz Günther', 'Marie-Terrell', 'Aneta Siemiginowska',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Christoph Deil',
                 'Brigitta Sipőcz', 'Johannes Buchner', 'Iva Laginja',
                 'Katrin Leinweber', 'nplee', 'Todd'],
        idval='5554957')

    add(version='4.13.1', title='sherpa/sherpa: Sherpa 4.13.1',
        date=todate(2021, 5, 18),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'dtnguyen2',
                 'Marie-Terrell', 'Hans Moritz Günther', 'Aneta Siemiginowska',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Christoph Deil',
                 'Brigitta Sipőcz', 'Johannes Buchner', 'Iva Laginja',
                 'Katrin Leinweber', 'nplee', 'Todd'],
        idval='4770623')
    add(version='4.13.0', title='sherpa/sherpa: Sherpa 4.13.0',
        date=todate(2021, 1, 8),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'dtnguyen2',
                 'Marie-Terrell', 'Hans Moritz Günther', 'Aneta Siemiginowska',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Christoph Deil',
                 'Brigitta Sipőcz', 'Johannes Buchner', 'Iva Laginja',
                 'Katrin Leinweber', 'nplee', 'Todd'],
        idval='4428938')

    add(version='4.12.2', title='sherpa/sherpa: Sherpa 4.12.2',
        date=todate(2020, 10, 27),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'dtnguyen2',
                 'Hans Moritz Günther', 'Marie-Terrell', 'Aneta Siemiginowska',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Christoph Deil',
                 'Brigitta Sipőcz', 'Johannes Buchner', 'Iva Laginja',
                 'Katrin Leinweber', 'nplee', 'Todd'],
        idval='4141888')
    add(version='4.12.1', title='sherpa/sherpa: Sherpa 4.12.1',
        date=todate(2020, 7, 14),
        authors=['Doug Burke', 'Omar Laurino', 'wmclaugh', 'dtnguyen2',
                 'Marie-Terrell', 'Hans Moritz Günther', 'Jamie Budynkiewicz',
                 'Aneta Siemiginowska', 'Tom Aldcroft', 'Christoph Deil',
                 'Brigitta Sipőcz', 'Katrin Leinweber', 'Todd'],
        idval='3944985')
    add(version='4.12.0', title='sherpa/sherpa: Sherpa 4.12.0',
        date=todate(2020, 1, 30),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2', 'wmclaugh',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'Christoph Deil', 'Marie-Terrell', 'Brigitta Sipőcz',
                 'Hans Moritz Günther', 'Todd', 'Katrin Leinweber'],
        idval='3631574')

    add(version='4.11.1', title='sherpa/sherpa: Sherpa 4.11.1',
        date=todate(2019, 8, 1),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'Christoph Deil', 'wmclaugh', 'Brigitta Sipocz',
                 'Katrin Leinweber'],
        idval='3358134')
    add(version='4.11.0', title='sherpa/sherpa: Sherpa 4.11.0',
        date=todate(2019, 2, 20),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'Christoph Deil', 'wmclaugh', 'Brigitta Sipocz',
                 'Katrin Leinweber'],
        idval='2573885')

    add(version='4.10.2', title='sherpa/sherpa: Sherpa 4.10.2',
        date=todate(2018, 12, 14),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'Christoph Deil', 'wmclaugh', 'Brigitta Sipocz',
                 'Katrin Leinweber'],
        idval='2275738')
    add(version='4.10.1', title='sherpa/sherpa: Sherpa 4.10.1',
        date=todate(2018, 10, 16),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'wmclaugh', 'Brigitta Sipocz', 'Christoph Deil',
                 'Katrin Leinweber'],
        idval='1463962')
    add(version='4.10.0', title='sherpa/sherpa: Sherpa 4.10.0',
        date=todate(2018, 5, 11),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'wmclaugh', 'Christoph Deil', 'Brigitta Sipocz'],
        idval='1245678')

    add(version='4.9.1', title='sherpa/sherpa: Sherpa 4.9.1',
        date=todate(2017, 8, 3),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'wmclaugh', 'Christoph Deil', 'Brigitta Sipocz'],
        idval='838686')
    add(version='4.9.0', title='sherpa/sherpa: Sherpa 4.9.0',
        date=todate(2017, 1, 27),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'wmclaugh', 'Christoph Deil', 'Brigitta Sipocz'],
        idval='260416')

    add(version='4.8.2', title='sherpa/sherpa: Sherpa 4.8.2',
        date=todate(2016, 9, 23),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2',
                 'Jamie Budynkiewicz', 'Tom Aldcroft', 'Aneta Siemiginowska',
                 'Christoph Deil', 'wmclaugh', 'Brigitta Sipocz'],
        idval='154744')
    add(version='4.8.1', title='sherpa: Sherpa 4.8.1',
        date=todate(2016, 4, 15),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2', 'Tom Aldcroft',
                 'Jamie Budynkiewicz', 'Aneta Siemiginowska', 'Christoph Deil',
                 'Brigitta Sipocz'],
        idval='49832')
    add(version='4.8.0', title='sherpa: Sherpa 4.8.0',
        date=todate(2016, 1, 27),
        authors=['Doug Burke', 'Omar Laurino', 'dtnguyen2', 'Tom Aldcroft',
                 'Aneta Siemiginowska'],
        idval='45243')

    kwargs = cite.get(version, None)
    if kwargs is None:
        return None

    return _make_citation(**kwargs, latest=False)


def _download_json(url: str) -> dict[str, Any]:
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
        dbg(f"Unable to access {url}: {he}")
        return {'failed': 'Unable to access the Zenodo site.'}

    try:
        jsdata = json.load(res)
    except UnicodeDecodeError as ue:
        dbg(f"Unable to decode JSON from Zenodo: {ue}")
        return {'failed': 'Unable to understand the response from Zenodo.'}

    return {'success': jsdata}


def _make_zenodo_citation(jsdata,
                          latest: bool = True) -> dict[str, str]:
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
        dbg(f"Unable to find metadata: {ke}")
        return {'failed': 'Unable to parse the Zenodo response.'}

    created = created.split('T')[0]
    isoformat = '%Y-%m-%d'
    try:
        date = datetime.datetime.strptime(created, isoformat)
    except ValueError:
        dbg(f"Unable to convert created: '{created}'")
        return {'failed': 'Unable to parse the Zenodo response.'}

    try:
        version = mdata['version']
        title = mdata['title']
        creators = mdata['creators']
    except KeyError as ke:
        dbg(f"Unable to find metadata: {ke}")
        return {'failed': 'Unable to parse the Zenodo response.'}

    authors = [c['name'] for c in creators]
    out = _make_citation(version=version, title=title, date=date,
                         authors=authors, idval=idval,
                         latest=latest)
    return {'success': out}


def _get_citation_version() -> str:
    """What version of Sherpa are we using?

    Returns
    -------
    version : str
        The Sherpa version.
    """

    out = f'You are using Sherpa {__version__}'
    if '+' in __version__:
        out += " (it is not a released version)"

    out += ".\n\n"
    return out


def _get_citation_zenodo_failure(failed: str) -> str:
    """Standard response when there's a problem.

    Returns
    -------
    text : str
        The failure message.

    """
    out = 'There was a problem retrieving the data from Zenodo:\n'
    out += failed
    out += '\n\n'
    out += DEFAULT_CITATION
    return out


def _get_citation_zenodo_latest() -> str:
    """Query Zenodo for the latest release.

    Returns
    -------
    citation : str
        Citation information. It will include information
        on any failure.
    """

    # Can we retrieve the information from Zenodo?
    #
    # This used to access the information via
    #     https://zenodo.org/api/records/593753
    # but this no-longer works and I do not trust Zenodo not
    # to change things again, so this is no a simplified
    # version of _download_zenodo_data().
    #
    from urllib.parse import urlencode

    params = {"q": "parent.id:593753",
              "all_versions": 0,
              "sort": "mostrecent"}
    paramstr = urlencode(params)
    url = f'https://zenodo.org/api/records?{paramstr}'

    dbg(f"Zenodo query: {url}")
    jsdata = _download_json(url)
    if 'failed' in jsdata:
        return _get_citation_zenodo_failure(jsdata['failed'])

    try:
        hit = jsdata['success']['hits']['hits'][0]
    except (KeyError, IndexError):
        dbg("Unable to find hits/hits[0]")
        return _get_citation_zenodo_failure('Unable to parse the Zenodo response')

    out = _make_zenodo_citation(hit)
    if 'failed' in out:
        return _get_citation_zenodo_failure(out['failed'])

    return out['success']


def _zenodo_missing(version: str) -> str:
    """

    Parameters
    ----------
    version : str
        The version we are looking for.

    Returns
    -------
    message : str
        The "unable to find" message.
    """

    return f'Zenodo has no information for version {version}.'


def _parse_zenodo_data(jsdata, version: str) -> dict[str, Any]:
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
        dbg(f'Version {version} not found')
        return {'failed': _zenodo_missing(version)}

    return {'success': data}


def _download_zenodo_data(version: str) -> dict[str, Any]:
    """Query Zenodo for the specific release.

    We have to deal with pagination in the Zenodo response.

    Parameters
    ----------
    version : str
        The release number (e.g. '4.12.2').

    Returns
    -------
    response : dict
        The data in JSON as 'success' or an error message as 'failed'.

    """

    # We could set the size parameter to something very large, to
    # get all responses with one call, but instead we use pagination.
    # Zenodo helpfully provides a links/next record with the
    # next URL, but it seems to be missing the all_versions=True
    # option, which makes it less-than-useful, hence the addition
    # of it below. The alternative would be to manually track the
    # page counter and add '&page=n' to the call.
    #
    from urllib.parse import urlencode

    params = {"q": "parent.id:593753",
              "all_versions": 1,
              "sort": "mostrecent"}
    paramstr = urlencode(params)
    url = f'https://zenodo.org/api/records?{paramstr}'

    missing = _zenodo_missing(version)

    while True:
        dbg(f"Zenodo query: {url}")
        jsdata = _download_json(url)

        # If the query fails then we error out
        #
        if 'failed' in jsdata:
            return jsdata

        jsdata = jsdata['success']
        data = _parse_zenodo_data(jsdata, version)
        if 'success' in data:
            return data

        # There are two failures we care about:
        #  - we can parse the information but have not been able to
        #    find the version
        #  - any other reason
        #
        # If the former then we look for the links/next entry to
        # look at the next page of the response. If it doesn't
        # exist we assume we are on the last page and so can error
        # out.
        #
        # If the latter then we error out rather than trying anything
        # else.
        #
        if data['failed'] != missing:
            return data

        try:
            url = jsdata['links']['next']

            # Add in the necessary all_versions tag: see
            # https://github.com/zenodo/zenodo/issues/1662
            #
            if 'all_versions=True' not in url:
                url += '&all_versions=True'

        except KeyError:
            # There is no links/next field so assume we've checked all
            # pages
            return data


def _get_citation_zenodo_version(version: str) -> str:
    """Query Zenodo for the specific release.

    As this has to return all Sherpa records it is slow.

    Parameters
    ----------
    version : str
        The release number (e.g. '4.12.2').

    Returns
    -------
    citation : str
        Citation information.
    """

    jsdata = _download_zenodo_data(version)
    if 'failed' in jsdata:
        return _get_citation_zenodo_failure(jsdata['failed'])

    out = _make_zenodo_citation(jsdata['success'], latest=False)
    if 'failed' in out:
        return _get_citation_zenodo_failure(out['failed'])

    return out['success']


def _get_citation(version: str = 'current') -> str:
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


def citation(version: str = 'current',
             filename=None,
             clobber: bool = False) -> None:
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


def get_include() -> str:
    "Get the root path for installed Sherpa header files"

    base = importlib.resources.files("sherpa")
    return str(base / 'include')


def get_config() -> str:
    "Get the path for the installed Sherpa configuration file"

    filename = "sherpa-standalone.rc"

    # The behavior depends on whether the NOSHERPARC
    # environment variable is set.
    #
    if 'NOSHERPARC' not in os.environ:

        # If SHERPARC is set, try that
        #
        config = os.environ.get('SHERPARC')
        if config is not None:
            if os.path.isfile(config):
                return config

        # Can we use the HOME directory?
        #
        home_dir = os.environ.get('HOME')
        if home_dir is not None:
            config = os.path.join(home_dir, f'.{filename}')
            if os.path.isfile(config):
                return config

    # Fall back to the system config file
    base = importlib.resources.files("sherpa")
    return str(base / filename)


def smoke(verbosity: int = 0,
          require_failure: bool = False,
          fits: Optional[str] = None,
          xspec: bool = False,
          ds9: bool = False) -> None:
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


def _install_test_deps() -> list[str]:
    """Check the needed modules for running the tests are installed.

    Since pytest has historically had issues with being able to use
    packages installed during this check, the list of such modules is
    returned so it can be passed to pytest.main().

    Returns
    -------
    plugins : list of module
        The pytest modules that were added by this call
        (may be empty).

    """

    def install(package_name):
        try:
            subprocess.call([sys.executable, '-m', 'pip', 'install', package_name],
                            stdout=sys.stdout, stderr=sys.stderr)
        except:
            print("""Cannot import pip or install packages with it.

            You need pytest in order to run the tests.
            If you downloaded the source code, please run 'pip install -r test_requirements.txt'
            from the source directory first.
            """)
            raise

    # packages are stored as a dictionary with keys
    #   name:  package name, as used by pip
    #   check: module to check it is installed
    #          (defaults to name if not given)
    #   constraint: the constraint to use with pip
    #               (defaults to name if not given)
    #
    # So we can have
    #
    #  {'name': 'pycheck'}
    #  {'name': 'pycheck', 'constraint': 'pytest>=5.0,!=5.2.3'}
    #  {'name': 'pytest-doctestplus',
    #   'check': 'pytest_doctestplus.output_checker'}
    #
    # Note that for pytest plugins, the module name is assumed to be
    # the name field with hyphens replaced by underscores. If necessary
    # this could be updated to be another field in the package
    # dictionary.
    #
    # The reason for the "check" field is that mid 2023 it was found
    # that
    #
    #    importlib.import_module("pytest_doctestplus")
    #
    # would work even when not installed, hence the need to import a
    # module within the package as a check. This has been left in even
    # though the current list of required plugins is empty.
    #
    deps: list[dict[str, str]] = [{'name': 'pytest',
                                   'constraint': 'pytest>=8.0'}]
    pytest_plugins: list[dict[str, str]] = []

    def get(dep: dict[str, str]) -> tuple[str, str, str]:
        name = dep['name']
        constraint = dep.get('constraint', name)
        check = dep.get('check', name)
        return name, constraint, check

    installed_plugins = []
    for dep in deps:
        name, constraint, check = get(dep)
        try:
            importlib.import_module(check)
        except ImportError:
            install(constraint)

    for dep in pytest_plugins:
        name, constraint, check = get(dep)
        try:
            importlib.import_module(check)
        except ImportError:
            install(constraint)
            installed_plugins.append(name.replace('-', '_'))

    return installed_plugins


def clitest() -> None:
    """The sherpa_test endpoint."""

    plugins = _install_test_deps()
    import pytest

    # Add in command-line arguments to allow configuring the Sherpa tests
    base = importlib.resources.files("sherpa")
    args = [str(base), '-rs'] + sys.argv[1:]

    # passing the plugins that have been installed "now".
    errno = pytest.main(args, plugins=plugins)
    sys.exit(errno)
