#!/usr/bin/env python

"""
Usage:

  ./make_zenodo_release.py [version]

Aim:

Create the release information for the sherpa._get_citation_hardcoded
routine. The output is to the screen, for manual inclusion.

If version is not given, all entries are included (which may
include CIAO releases), otherwise a version number like '4.12.2'
is expected.

"""

import datetime
import json
import ssl
from urllib.parse import urlencode
from urllib.request import Request, urlopen


def download_json(url):
    """Return the JSON data.

    Parameters
    ----------
    url : str
        The URL to process.

    Returns
    -------
    response : dict
        The data in JSON.

    """

    req = Request(url)
    context = ssl._create_unverified_context()

    # Do not bother trapping errors
    res = urlopen(req, context=context)
    return json.load(res)


def make_zenodo_citation(jsdata):
    """Convert Zenodo record to a citation.

    Prints the record to stdout.

    Parameters
    ----------
    jsdata : object
        A Zenodo record.

    """

    created = jsdata['created']
    mdata = jsdata['metadata']
    idval = jsdata['id']

    created = created.split('T')[0]
    isoformat = '%Y-%m-%d'
    date = datetime.datetime.strptime(created, isoformat)

    version = mdata['version']
    title = mdata['title']
    creators = mdata['creators']

    authors = [c['name'] for c in creators]

    print(f"    add(version='{version}', title='{title}',")
    print(f"        date=todate({date.year}, {date.month}, {date.day}),")
    print(f"        authors={authors},")
    print(f"        idval='{idval}')")


def dump_zenodo(version):
    """Access Zenodo for the version information.

    Parameters
    ----------
    version : str or None
        The version string (e.g. '4.12.2'). If None then all versions
        are reported.
    """

    params = {"q": "parent.id:593753",
              "all_versions": 1,
              "sort": "mostrecent"}
    paramstr = urlencode(params)
    url = f'https://zenodo.org/api/records?{paramstr}'
    while True:
        jsdata = download_json(url)
        for hit in jsdata['hits']['hits']:
            if version is None or hit['metadata']['version'] == version:
                make_zenodo_citation(hit)
                if version is not None:
                    return

        try:
            url = jsdata['links']['next']

            # Add in the necessary all_versions tag: see
            # https://github.com/zenodo/zenodo/issues/1662
            #
            if 'all_versions=True' not in url:
                url += '&all_versions=True'

        except KeyError:
            return


if __name__ == "__main__":

    import sys
    if len(sys.argv) > 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} [version]\n")
        sys.exit(1)

    if len(sys.argv) == 2:
        version = sys.argv[1]
    else:
        version = None

    dump_zenodo(version)
