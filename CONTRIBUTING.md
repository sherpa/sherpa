Contributing to Sherpa
======================

Thanks for reporting an issue, fixing a bug, or adding a new feature to
Sherpa. The following guidelines are not meant to be too onerous; please open
an issue if you find this to be wrong! Thanks are due to the [AstroPy
project](http://www.astropy.org/), from which we based this document.

Reporting Issues
----------------

When [opening an issue](https://github.com/sherpa/sherpa/issues) to report a
problem, please try and provide a minimal code example that reproduces the
issue, and also include details of the operating system, and the versions of
the code you are using - in particular for Sherpa, NumPy, and Python. The
Sherpa version can be found using

    $ python -c 'import sherpa; print(sherpa.__version__)'
    40701

When using a development version of Sherpa, the value will include the
Git commit. An example of this is ``4.7+524.g2d013d1``.

Contributing code
-----------------

Contributions to Sherpa are generally done via a [pull
request](https://github.com/sherpa/sherpa/pulls) from your GitHub fork of the
[Sherpa repository][https://github.com/sherpa/sherpa), although you will see
core developers using the Sherpa repository for development. At present we do
not have a document describing our suggested development workflow, so we
suggest you follow the [AstroPy
workflow](http://docs.astropy.org/en/latest/development/workflow/development_workflow.html)
for now. We also do not have a written up set of coding guidelines, so we also
suggest using the [AstroPy
guidelines](http://docs.astropy.org/en/latest/development/codeguide.html); the
current Sherpa code base is slowly being converted to a single style.

Once you open a pull request (which should be opened against the ``master``
branch, not against any of the other branches), please make sure that you
include the following:

- **Code**: the code you are adding.

- **Tests**: these are usually tests that ensures that code that previously
  failed now works (regression tests) or tests that cover as much as possible
  of the new functionality to make sure it doesn't break in future, and also
  returns consistent results on all platforms (since we run these tests on many
  platforms/configurations). If you have questions about how the tests in
  Sherpa are arranged, please include them in your pull request.

- **Documentation**: at present, the documentation is limited to Python
  docstrings, but it is expected that more-extensive documentation will
  be added soon (in the ``docs/`` directory).

- **Changelog entry**: unlike AstroPy, we do not have a "change log" file.
  Instead, we ask that the Pull Request include the relevant information,
  in a ``Release Note`` section at the start of the request (please see
  the [accepted requests](https://github.com/sherpa/sherpa/pulls?q=is%3Apr+is%3Aclosed)
  for examples of style and format.
  
Other Tips
----------

When contributing trivial documentation fixes (i.e. fixes to typos, spelling,
grammar) that do not contain any special markup and are not associated with
code changes, include the string "[skip ci]" at the end of your commit
message.  For example:

    $ git commit -m "Fixed typo [skip ci]"

This will prevent automated tests for running against your change, freeing
up resources for testing non-trivial changes.

If you already made the commit without including this string, you can edit
your existing commit message by running:

    $ git commit --amend

Checklist for Contributed Code
------------------------------

A pull request for a new feature will be reviewed to see if it meets the
following requirements.  For any pull request, a Sherpa maintainer can help to
make sure that the pull request meets the requirements for inclusion in the
package. A good first test is to ensure that you can check out your new branch
and run

    $ python setup.py develop
    $ python setup.py test

without a failure.

**Git branch**

The changes should be made in a branch which is based off of the ``master``
branch of the [Sherpa repository](https://github.com/sherpa/sherpa), and
should have a descriptive name. For bug fixes, the name should include the
bug ID and a short description of the bug fix. For new - or updated - 
functionality, the name should be descriptive.

Example: fix for [bug #64](https://github.com/sherpa/sherpa/issues/64):

    % git remote add upstream https://github.com/sherpa/sherpa
    % git fetch upstream
    % git checkout -b bug-#64-comparison-to-None upstream/master

Example: new functionality, such as
[Setuptools not required](https://github.com/sherpa/sherpa/pull/65):

    % git remote add upstream https://github.com/sherpa/sherpa
    % git fetch upstream
    % git checkout -b setuptools-not-required upstream/master

**Software versions**

At present Sherpa is built using Python 2.7, but there are plans to support
Python 3 in the near future (this will be tracked in [issue
#76](https://github.com/sherpa/sherpa/issues/76)).

Ideally, NumPy support should be 1.6 or greater, but please include a comment
if you need to restrict (or relax) this further.

Currently Sherpa requires both a C and Fortran compiler, along with
related tools, to build. Additional dependencies must be highlighted,
and should be made optional if at all possible.

**Testing**

The Travis Continuous Integration service is used to test all pull requests
(examples can be found at the [Sherpa
page](https://travis-ci.org/sherpa/sherpa/)), and these tests must pass
before a pull request can be considered for inclusion in Sherpa.

Ideally each pull request will come with new, or updated, tests to check
the change. This may be relaxed when fixes to existing functionality are
made, but please indicate that this is the case in the pull request.

**Documentation**

Does the Pull Request include the following (for simple fixes, such as
a typo in the documentation, not all are required):

 - a complete description of the feature
 - the design outline
 - the test cases exercised by the test suite

Do the docstrings for any  new, or updated, code describe what the
code does, the format of the inputs and outputs, references to
further reading or the original algorithm, exceptions raised, and
examples?

**License**

The Sherpa code is licensed under the Gnu Public License, version 3
or later, so any changes must not conflict with this.
