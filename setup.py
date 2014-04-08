#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from setup_helpers import setup


# COMMON METADATA
meta = {
        'name' : 'sherpa',
        'version' : '4.7b1',
        'author' : 'Smithsonian Astrophysical Observatory / Chandra X-Ray Center',
        'author_email' : 'cxchelp@head.cfa.harvard.edu',
        'url' : 'http://cxc.harvard.edu/sherpa/',
        'description' : 'Modeling and fitting package for scientific data analysis',
        'license' : 'GNU GPL v3',
        'long_description' : 'Modeling and fitting package for scientific data analysis',
        'platforms' : 'Linux, Mac OS X',
        'install_requires' : ['numpy', 'pyfits', 'matplotlib'],

        }

meta['packages'] = packages=['sherpa',
                'sherpa.estmethods',
                'sherpa.image',
                'sherpa.models',
                'sherpa.optmethods',
                'sherpa.plot',
                'sherpa.sim',
                'sherpa.stats',
                'sherpa.ui',
                'sherpa.utils',
                'sherpa.astro',
                'sherpa.astro.io',
                'sherpa.astro.models',
                'sherpa.astro.optical',
                'sherpa.astro.sim',
                'sherpa.astro.ui',
                'sherpa.astro.utils',
                ]

meta['package_data'] = {'sherpa': ['include/sherpa/*.hh',
                               'include/sherpa/astro/*.hh',
                               'tests/test_*.py'],
                    'sherpa.estmethods': ['tests/test_*.py'],
                    'sherpa.image': ['tests/test_*.py'],
                    'sherpa.models': ['tests/test_*.py'],
                    'sherpa.optmethods': ['tests/test_*.py'],
                    'sherpa.plot': ['tests/test_*.py'],
                    'sherpa.sim': ['tests/test_*.py'],
                    'sherpa.stats': ['tests/test_*.py'],
                    'sherpa.ui': ['tests/test_*.py'],
                    'sherpa.utils': ['tests/test_*.py'],
                    'sherpa.astro': ['tests/test_*.py'],
                    'sherpa.astro.io': ['tests/test_*.py'],
                    'sherpa.astro.models': ['tests/test_*.py'],
                    'sherpa.astro.optical': ['tests/test_*.py'],
                    'sherpa.astro.sim': ['tests/test_*.py'],
                    'sherpa.astro.ui': ['tests/test_*.py'],
                    'sherpa.astro.utils': ['tests/test_*.py'],
                    }

meta['data_files'] = data_files = [('sherpa', ['sherpa/sherpa.rc']),
#            ('', ['build/lib/group.so']),
            ]

setup(**meta)