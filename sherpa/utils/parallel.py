#
#  Copyright (C) 2007, 2015, 2016, 2018, 2019, 2020, 2021, 2022, 2023
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
Routines for running code in parallel.
"""

from configparser import ConfigParser, NoSectionError
import logging

import numpy as np

from sherpa import get_config

warning = logging.getLogger("sherpa").warning

config = ConfigParser()
config.read(get_config())

_ncpu_val = "NONE"
try:
    _ncpu_val = config.get('parallel', 'numcores').strip().upper()
except NoSectionError:
    pass

_ncpus = None
"""The number of CPU cores to use when running jobs in parallel.

This is taken from the parallel.numcores setting from the Sherpa
configuration file (returned by `sherpa.get_config()`), where the
default setting of ``None`` will use all available cores.
"""

if not _ncpu_val.startswith('NONE'):
    _ncpus = int(_ncpu_val)

_multi = False
"""Can jobs be run in parallel?

The ability to run jobs in parallel depends on whether the Python
`multiprocessing` module can be configured to use the
multiprocessing.multiprocessing_start_method setting from the
Sherpa configuration file (returned by `sherpa.get_config()`).
"""

try:
    import multiprocessing

    multiprocessing_start_method = config.get('multiprocessing', 'multiprocessing_start_method', fallback='fork')

    if multiprocessing_start_method not in ('fork', 'spawn', 'default'):
        raise ValueError('multiprocessing_start method must be one of "fork", "spawn", or "default"')

    if multiprocessing_start_method != 'default':
        multiprocessing.set_start_method(multiprocessing_start_method, force=True)

    _multi = True

    if _ncpus is None:
        _ncpus = multiprocessing.cpu_count()

except Exception as e:
    warning("parallel processing is unavailable,\n"
            f"multiprocessing module failed with \n'{e}'")
    _ncpus = 1
    _multi = False

del _ncpu_val, config, get_config, ConfigParser, NoSectionError


__all__ = ("_multi", "_ncpus",
           "parallel_map", "parallel_map_funcs",
           "run_tasks")


# Can this be replaced by itertools.batched once Python 3.12 is the
# minimum supported version? The fact that it may group elements
# differently should not be a problem for downstream users.
#
def split_array(arr, m):
    """Split array ``arr`` into ``m`` roughly equal chunks
    >>> split_array(range(27), 6)
    [[0, 1, 2, 3, 4],
     [5, 6, 7, 8],
     [9, 10, 11, 12, 13],
     [14, 15, 16, 17],
     [18, 19, 20, 21, 22],
     [23, 24, 25, 26]]

    >>> import numpy as np
    >>> split_array(np.arange(25), 6)
    [array([0, 1, 2, 3]),
     array([4, 5, 6, 7]),
     array([ 8,  9, 10, 11, 12]),
     array([13, 14, 15, 16]),
     array([17, 18, 19, 20]),
     array([21, 22, 23, 24])]

    >>> split_array(np.arange(30).reshape(5,-1), 3)
    [array([[ 0,  1,  2,  3,  4,  5],
           [ 6,  7,  8,  9, 10, 11]]),
    array([[12, 13, 14, 15, 16, 17]]),
    array([[18, 19, 20, 21, 22, 23],
          [24, 25, 26, 27, 28, 29]])]

    Author: Tom Aldcroft
      split_array() - originated from Python users working group
    """
    n = len(arr)
    idx = [int(round(i * n / float(m))) for i in range(m + 1)]
    return [arr[idx[i]:idx[i + 1]] for i in range(m)]


def worker(f, ii, chunk, out_q, err_q):
    """Evaluate a function for each element, add response to queue.

    Parameters
    ----------
    f : callable
        The function. It accepts a single argument, matching the
        elements of chunk.
    ii : int
        The identifier for this worker.
    chunk : sequence
        The elements to be passed to func.
    out_q, err_q : manager.Queue
        The success channel - which will be sent (idx, retval) on
        success - and the error channel, which will be sent any
        exception.

    """
    try:
        vals = [f(c) for c in chunk]
    except Exception as e:
        err_q.put(e)
        return

    # output the result and task ID to output queue
    out_q.put((ii, vals))


def process_tasks(procs, err_q):
    """Ensure all the processes are run error-ing out if needed.

    Parameters
    ----------
    procs : list of multiprocessing.Process tasks
        The processes to run.
    err_q : manager.Queue
        The error channel used by the processes.

    """

    def die():
        """Clean up processes that are still running"""
        for proc in procs:
            if proc.exitcode is None:
                proc.terminate()

    try:
        for proc in procs:
            proc.start()

        for proc in procs:
            proc.join()

    except KeyboardInterrupt as e:
        # kill all slave processes on ctrl-C
        die()
        raise e

    if not err_q.empty():
        die()
        raise err_q.get()


def run_tasks(procs, err_q, out_q, num=None):
    """Run the processes, exiting early if necessary, and return the results.

    .. versionchanged:: 4.16.0
       The num argument is not needed and has been marked optional.
       It will be removed in a future release.

    Parameters
    ----------
    procs : list of multiprocessing.Process tasks
        The processes to run.
    err_q, out_q : manager.Queue
        The error and success channels used by the processes.
    num : optional
        This argument is unused and will be removed.

    Returns
    -------
    result : list
        The result from the processes. This may contain more elements
        than procs, as each process may return multiple results.

    Notes
    -----
    Each process sends its output - a pair with index and a list of
    results - to the out_q queue, and any error encounteded to the
    err_q queue. There should be len(procs) messages sent to the out_q
    queue for a successful run.

    """

    process_tasks(procs, err_q)

    # Loop through and insert the results to match the original order.
    #
    results = [None] * len(procs)
    while not out_q.empty():
        idx, result = out_q.get()
        results[idx] = result

    # Since each process may contain multiple results, flatten the
    # results.
    #
    vals = []
    for r in results:
        vals.extend(r)

    return vals


def parallel_map(function, sequence, numcores=None):
    """Run a function on a sequence of inputs in parallel.

    A parallelized version of the native Python map function that
    utilizes the Python multiprocessing module to divide and
    conquer sequence.

    Parameters
    ----------
    function : function
       This function accepts a single argument (an element of
       ``sequence``) and returns a value.
    sequence : array_like
       The data to be passed to ``function``.
    numcores : int or None, optional
       The number of calls to ``function`` to run in parallel. When
       set to ``None``, all the available CPUs on the machine - as
       set either by the 'numcores' setting of the 'parallel' section
       of Sherpa's preferences or by multiprocessing.cpu_count - are
       used.

    Returns
    -------
    ans : array
       The return values from the calls, in the same order as the
       ``sequence`` array.

    Notes
    -----
    A tuple or dictionary should be used to pass multiple values to
    the function.

    The input list is split into ``numcores`` chunks, and then each
    chunk is run in parallel. There is no guarantee to the ordering
    of the tasks.

    Examples
    --------

    In the following examples a simple set of computations are used;
    in reality the function is expected to be run on computations
    that take a significant amount of time to run.

    Run the computation (summing up each element of the input array)
    on a separate core and return the results (unless the machine only
    has a single core or the parallel.numcores setting is set to 1).

    >>> import numpy as np
    >>> args = [np.arange(5), np.arange(3), np.arange(7)]
    >>> parallel_map(np.sum, args)
    [10, 3, 21]

    Use two jobs to evaluate the results: one job will sum up two arrays
    while the other will only sum one array since there are 3 jobs to
    run.

    >>> parallel_map(np.sum, args, numcores=2)
    [10, 3, 21]

    An example of sending in multiple arguments to a function (``comp``)
    via a dictionary (although in this case there is only one task to
    execute):

    >>> parallel_map(comp, [{'idx1': 23, 'idx2': 47}])

    Here the ``tcomp`` function accepts a single parameter which it
    can deconstruct to extract the two values it needs:

    >>> parallel_map(tcomp, [(23, 47), (2, 20), (5, 10)])

    """
    if not callable(function):
        raise TypeError(f"input function '{repr(function)}' is not callable")

    if not np.iterable(sequence):
        raise TypeError(f"input '{repr(sequence)}' is not iterable")

    size = len(sequence)

    if not _multi or size == 1 or (numcores is not None and numcores < 2):
        return list(map(function, sequence))

    if numcores is None:
        numcores = _ncpus

    # Returns a started SyncManager object which can be used for sharing
    # objects between processes. The returned manager object corresponds
    # to a spawned child process and has methods which will create shared
    # objects and return corresponding proxies.
    manager = multiprocessing.Manager()

    # Create FIFO queue and lock shared objects and return proxies to them.
    # The managers handles a server process that manages shared objects that
    # each slave process has access to.  Bottom line -- thread-safe.
    out_q = manager.Queue()
    err_q = manager.Queue()
    # lock = manager.Lock() - currently unused

    # if sequence is less than numcores, only use len sequence number of
    # processes
    if size < numcores:
        numcores = size

    # group sequence into numcores-worth of chunks
    sequence = split_array(sequence, numcores)

    procs = [multiprocessing.Process(target=worker,
                                     args=(function, ii, chunk, out_q, err_q))
             for ii, chunk in enumerate(sequence)]

    return run_tasks(procs, err_q, out_q)


# TODO: this routine needs a review
def parallel_map_funcs(funcs, datasets, numcores=None):
    """Run a sequence of function on a sequence of inputs in parallel.

    Sherpa's parallel_map runs a single function to an iterable set of
    sequence.  parallel_map_funcs is generalized parallelized version
    of sherpa's parallel_map function since each element of the ordered
    iterable funcs shall operate on the each element of the datasets.

    Parameters
    ----------
    funcs : a list or tuple of functions
       An ordered iterable sequence of functions which accepts an element
       of the datasets and returns a value.  The number of elements in
       funcs must match the number of elements of the datasets.
    datasets : a list or tuple of array_like
       The data to be passed to ``func``. The number of elements in
       datasets must match the number of elements of funcs.
    numcores : int or None, optional
       The number of calls to ``funcs`` to run in parallel. When
       set to ``None``, all the available CPUs on the machine - as
       set either by the 'numcores' setting of the 'parallel' section
       of Sherpa's preferences or by multiprocessing.cpu_count - are
       used.

    Returns
    -------
    ans : array
       The return values from the calls, in the same order as the
       ``sequence`` array.

    Notes
    -----
    Due to the overhead involved in passing the functions and datasets
    to the different cores, the functions should be very time consuming
    to compute (of order 0.1-1s).  This is similar to the ``parallel_map``
    function.

    An ordered iterable (i.e. tuple or list) should be used to pass multiple
    values to the multiple functions. The lengths of the iterable funcs and
    datasets must be equal. The corresponding funcs and datasets are passed
    to the different cores to distribute the work in parallel. There is no
    guarantee to the ordering of the tasks.

    Examples
    --------

    In the following examples a simple set of computations, sum and std
    deviations, are used; in reality the function is expected to be run
    on computations that take a significant amount of time to run.

    Run the computation (summing up each element of the first input array
    and calculate the standard deviation of the second input array)
    on a separate core and return the results (unless the machine only
    has a single core or the parallel.numcores setting is set to 1).

    >>> import numpy as np
    >>> funcs = [np.sum, np.std]
    >>> datasets = [np.arange(3), np.arange(4)]
    >>> parallel_map_funcs(funcs, datasets, numcores=2)
    [0, 1, 2, 0.0, 0.0, 0.0, 0.0]

    """
    if not np.iterable(funcs):
        raise TypeError(f"input '{repr(funcs)}' is not iterable")

    if not np.iterable(datasets):
        raise TypeError(f"input '{repr(datasets)}' is not iterable")

    for func in funcs:
        if not callable(func):
            raise TypeError(f"input func '{repr(func)}' is not callable")

    funcs_size = len(funcs)
    datasets_size = len(datasets)
    if funcs_size != datasets_size:
        raise TypeError(f"input funcs ({funcs_size}) and datsets "
                        "({datasets_size}) size must be same")

    if not _multi or datasets_size == 1 or \
            (numcores is not None and numcores < 2):
        # TODO: see issue #1743
        #
        return list(map(funcs[0], datasets))

    if numcores is None:
        numcores = _ncpus

    # Returns a started SyncManager object which can be used for sharing
    # objects between processes. The returned manager object corresponds
    # to a spawned child process and has methods which will create shared
    # objects and return corresponding proxies.
    manager = multiprocessing.Manager()

    # Create FIFO queue and lock shared objects and return proxies to them.
    # The manager handles a server process that manages shared objects that
    # each slave process has access to.  Bottom line -- thread-safe.
    out_q = manager.Queue()
    err_q = manager.Queue()
    # lock = manager.Lock() - currently unused

    procs = [multiprocessing.Process(target=worker,
                                     args=(funcs[ii], ii, chunk, out_q, err_q))
             for ii, chunk in enumerate(datasets)]

    return run_tasks(procs, err_q, out_q)
