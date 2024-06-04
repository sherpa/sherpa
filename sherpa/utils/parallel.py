#
#  Copyright (C) 2007, 2015, 2016, 2018 - 2024
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

"""Routines for running code in parallel.

.. versionchanged:: 4.16.1
   All `multiprocessing` calls are now done using an explicit context,
   available as the `context` field, rather than using the global
   version.

.. versionadded:: 4.16.0
   Prior to this, these symbols were provided by the
   :py:mod:`sherpa.utils` module.

"""

from abc import abstractmethod
from configparser import ConfigParser
import inspect
import logging
from typing import Any, Final, Optional, Protocol, Sequence, TypeVar

import numpy as np

from sherpa import get_config
from .random import RandomType

# A number of symbols have been added to this module in release 4.17.0
# to allow typing statements to be made. This is partly because the
# multiprocessing module, which provides a number of symbols that
# would be used in such statements, is optional (and so is the
# threading module which multiprocessing often uses to define
# symbols).  Therefore a number of protocols have been added to allow
# statements to make statements like "this requires a queue which
# handles exceptions".
#
I_contra = TypeVar("I_contra", contravariant=True)
O_co = TypeVar("O_co", covariant=True)
T = TypeVar("T")


class Callback(Protocol[I_contra, O_co]):
    """Simple callback of a single argument."""

    def __call__(self, arg: I_contra) -> O_co:
        ...


class CallbackWithRNG(Protocol[I_contra, O_co]):
    """Simple callback of a single argument plus RNG argument."""

    def __call__(self,
                 arg: I_contra,
                 rng: Optional[RandomType]) -> O_co:
        ...


# Since importing multiprocessing and threading are allowed to fail it
# is hard to type things with a class defined in these modules.
#
class SupportsProcess(Protocol):
    """Label those methods from multiprocessing.Process we need."""

    exitcode: Optional[int]

    @abstractmethod
    def start(self) -> None:
        ...

    # Technically join accepts an optional timeout argument but we do
    # not need that.
    #
    @abstractmethod
    def join(self) -> None:
        ...

    @abstractmethod
    def terminate(self) -> None:
        ...

# This typing rule is actually less generic than the actual queue
# type, since there's no guarantee that it is only used to
# send/receive a single type. However, the use here does have the
# queues having a fixed type.
#
class SupportsQueue(Protocol[T]):
    """Represent the multiprocessing Queue class"""

    @abstractmethod
    def empty(self) -> bool:
        ...

    # we do not need the block or timeout options
    @abstractmethod
    def put(self, obj: T) -> None:
        ...

    # we do not need the block or timeout options
    @abstractmethod
    def get(self) -> T:
        ...


class SupportsLock(Protocol):
    """Represent the multiprocessing Lock class"""

    @abstractmethod
    def acquire(self,
                block: bool = True,
                timeout: Optional[float] = None) -> bool:
        ...

    @abstractmethod
    def release(self) -> None:
        ...


class SupportsManager(Protocol):
    """Represent the multiprocessing Manager class"""

    @abstractmethod
    def Queue(self) -> SupportsQueue:
        ...

    @abstractmethod
    def Lock(self) -> SupportsLock:
        ...


class SupportsContext(Protocol):
    """The multiprocessing context.

    Note that the BaseContext does not provide the Process symbol but
    all the "runnable" variants do.

    """

    @property
    @abstractmethod
    def Process(self):
        ...

    @abstractmethod
    def Manager(self) -> SupportsManager:
        ...

    @abstractmethod
    def cpu_count(self) -> int:
        ...


debug = logging.getLogger(__name__).debug
warning = logging.getLogger(__name__).warning

config = ConfigParser()
config.read(get_config())

_ncpus = None

_ncpu_val = config.get('parallel', 'numcores', fallback="NONE").upper()
if not _ncpu_val.startswith('NONE'):
    _ncpus = int(_ncpu_val)

_multi = False

# This should be Optional[multiprocessing.context.BaseContext] but
# we do not require multiprocessing to be available.
#
_context : Optional[SupportsContext]

try:
    import multiprocessing

    multiprocessing_start_method = config.get('multiprocessing',
                                              'multiprocessing_start_method',
                                              fallback='fork')

    if multiprocessing_start_method not in ('fork', 'forkserver',
                                            'spawn', 'default'):
        raise ValueError('multiprocessing_start method must be one of "fork", "forkserver", "spawn", or "default"')


    if multiprocessing_start_method == "default":
        _context = multiprocessing.get_context()  # type: ignore[assignment]
    else:
        _context = multiprocessing.get_context(multiprocessing_start_method)  # type: ignore[assignment]

    _multi = True

    # The '_context is not None' check is added to simplify the type
    # checking.
    #
    if _ncpus is None and _context is not None:
        _ncpus = _context.cpu_count()

except Exception as e:
    warning("parallel processing is unavailable,\n"
            "multiprocessing module failed with \n'%s'", str(e))
    _ncpus = 1
    _multi = False
    _context = None

del _ncpu_val, config, get_config, ConfigParser


multi: Final[bool] = _multi
"""Can jobs be run in parallel?

The ability to run jobs in parallel depends on whether the Python
`multiprocessing` module can be configured to use the
``multiprocessing.multiprocessing_start_method`` setting from the
Sherpa configuration file (returned by `sherpa.get_config()`).

See Also
--------
context, ncpus

"""

assert _ncpus is not None  # to please typing
ncpus: Final[int] = _ncpus
"""The number of CPU cores to use when running jobs in parallel.

This is taken from the ``parallel.numcores`` setting from the Sherpa
configuration file (returned by `sherpa.get_config()`), where the
default setting of ``None`` will use all available cores.
"""

# This is hard to type given that multiprocessing is not a required
# module.
#
context = _context
"""The multiprocessing context used to run the processes.

This will be ``None`` when multiprocessing support is not available
(that is, `multi` is ``False``). It is set by the
``multiprocessing_start_method`` setting from the ``multiprocessing``
block in the configuration file (returned by `sherpa.get_config`).

.. versionadded:: 4.16.1

See Also
--------
multi

"""


__all__ = ("multi", "ncpus", "context",
           "parallel_map", "parallel_map_funcs", "parallel_map_rng",
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


def worker(f: Callback[I_contra, O_co],
           idx: int,
           chunk: Sequence[I_contra],
           out_q: SupportsQueue[tuple[int, list[O_co]]],
           err_q: SupportsQueue[Exception]) -> None:
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
    out_q.put((idx, vals))


def worker_rng(func: CallbackWithRNG[I_contra, O_co],
               idx: int,
               chunk: Sequence[I_contra],
               out_q: SupportsQueue[tuple[int, list[O_co]]],
               err_q: SupportsQueue[Exception],
               rng: Optional[RandomType]
               ) -> None:
    """Evaluate a function for each element, add response to queue.

    Unlike worker(), this also sends in a RNG.

    Parameters
    ----------
    func : callable
       The function. It accepts two arguments, the first being an
       element of chunk and the second the rng argument.
    idx : int
       The identifier for this worker.
    chunk : sequence
       The elements to be passed to func.
    out_q, err_q : manager.Queue
       The success channel - which will be sent (idx, retval) on
       success - and the error channel, which will be sent any
       exception.
    rng : numpy.random.Generator, numpy.random.RandomState, or None, optional
       Determines how random numbers are created. If set to None then
       the routines in `numpy.random` are used, and so can be
       controlled by calling `numpy.random.seed`. It is up to the
       caller to ensure that the generator will work correctly if
       called in parallel.

    """

    try:
        vals = [func(c, rng) for c in chunk]
    except Exception as e:
        err_q.put(e)
        return

    out_q.put((idx, vals))


def cleanup_tasks(procs: Sequence[SupportsProcess]) -> None:
    """Clean up processes that are still running"""
    for proc in procs:
        if proc.exitcode is None:
            proc.terminate()


def process_tasks(procs: Sequence[SupportsProcess],
                  err_q: SupportsQueue[Exception]
                  ) -> None:
    """Ensure all the processes are run error-ing out if needed.

    Parameters
    ----------
    procs : list of multiprocessing.Process tasks
        The processes to run.
    err_q : manager.Queue
        The error channel used by the processes.

    """

    try:
        for proc in procs:
            proc.start()

        for proc in procs:
            proc.join()

    except KeyboardInterrupt as e:
        # kill all slave processes on ctrl-C
        cleanup_tasks(procs)
        raise e

    if not err_q.empty():
        cleanup_tasks(procs)
        raise err_q.get()


def run_tasks(procs: Sequence[SupportsProcess],
              err_q: SupportsQueue[Exception],
              out_q: SupportsQueue[tuple[int, list[O_co]]],
              num: Optional[Any] = None) -> list[O_co]:
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
    results: list = [None] * len(procs)
    while not out_q.empty():
        idx, result = out_q.get()
        results[idx] = result

    # Since each process may contain multiple results, flatten the
    # results.
    #
    vals: list = []
    for r in results:
        # The assumption is that this task worked (i.e. returned a value).
        if r is None:
            raise RuntimeError("task failed")
        vals.extend(r)

    return vals




def parallel_map(function: Callback[I_contra, O_co],
                 sequence: Sequence[I_contra],
                 numcores: Optional[int] = None
                 ) -> list[O_co]:
    """Run a function on a sequence of inputs in parallel.

    A parallelized version of the native Python map function that
    utilizes the Python multiprocessing module to divide and conquer
    sequence. If ``function`` uses random numbers then
    `parallel_map_rng` should be used instead.

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

    See Also
    --------
    parallel_map_rng

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

    # Using np.iterable does not imply to mypy that you can use len,
    # so add the ignore call.
    #
    size = len(sequence)  # type: ignore[arg-type]

    ncores = ncpus if numcores is None else numcores
    if not _multi or size == 1 or ncores < 2:
        return list(map(function, sequence))

    # At this point we know context is not None but the typing code
    # does not.
    assert context is not None

    # Returns a started SyncManager object which can be used for sharing
    # objects between processes. The returned manager object corresponds
    # to a spawned child process and has methods which will create shared
    # objects and return corresponding proxies.
    manager = context.Manager()

    # Create FIFO queue and lock shared objects and return proxies to them.
    # The managers handles a server process that manages shared objects that
    # each slave process has access to.  Bottom line -- thread-safe.
    out_q = manager.Queue()
    err_q = manager.Queue()

    # if sequence is less than numcores, only use len sequence number of
    # processes
    if size < ncores:
        ncores = size

    # group sequence into numcores-worth of chunks
    sequence = split_array(sequence, ncores)

    assert context.Process is not None
    procs = [context.Process(target=worker,
                             args=(function, idx, chunk, out_q, err_q))
             for idx, chunk in enumerate(sequence)]

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
        raise TypeError(f"input funcs ({funcs_size}) and datasets "
                        "({datasets_size}) size must be same")

    if not _multi or datasets_size == 1 or \
            (numcores is not None and numcores < 2):
        # TODO: see issue #1743
        #
        return list(map(funcs[0], datasets))

    # At this point we know context is not None but the typing code
    # does not.
    assert context is not None

    if numcores is None:
        numcores = _ncpus

    # Returns a started SyncManager object which can be used for sharing
    # objects between processes. The returned manager object corresponds
    # to a spawned child process and has methods which will create shared
    # objects and return corresponding proxies.
    manager = context.Manager()

    # Create FIFO queue and lock shared objects and return proxies to them.
    # The manager handles a server process that manages shared objects that
    # each slave process has access to.  Bottom line -- thread-safe.
    out_q = manager.Queue()
    err_q = manager.Queue()

    assert context.Process is not None
    procs = [context.Process(target=worker,
                             args=(funcs[idx], idx, chunk, out_q, err_q))
             for idx, chunk in enumerate(datasets)]

    return run_tasks(procs, err_q, out_q)


# The typing is not quite right for function
#
def parallel_map_rng(function: CallbackWithRNG[I_contra, O_co],
                     sequence: Sequence[I_contra],
                     numcores: Optional[int] = None,
                     rng: Optional[RandomType] = None
                     ) -> list[O_co]:
    """Run a function on a sequence of inputs in parallel with a RNG.

    Similar to parallel_map, but the function takes two arguments,
    with the second one being ``rng``, the random generator to use. This
    is for those functions which need random numbers, and this routine
    takes care to create a separate generator for each process run in
    parallel.

    .. versionadded:: 4.16.0

    Parameters
    ----------
    function : function
       This function accepts two arguments - the first being an
       element of ``sequence`` and second called ``rng`` - and returns a
       value.
    sequence : array_like
       The data to be passed to ``function`` as the first argument.
    numcores : int or None, optional
       The number of calls to ``function`` to run in parallel. When
       set to ``None``, all the available CPUs on the machine - as
       set either by the 'numcores' setting of the 'parallel' section
       of Sherpa's preferences or by multiprocessing.cpu_count - are
       used.
    rng : numpy.random.Generator, numpy.random.RandomState, or None, optional
       Controls how random numbers are generated. When code is run in
       parallel, each worker is sent a separate generator, to ensure
       that the sequences are different, and the rng parameter is used
       to create the seed number passed to `numpy.random.SeedSequence`
       for this case.

    Returns
    -------
    ans : array
       The return values from the calls, in the same order as the
       ``sequence`` array.

    See Also
    --------
    parallel_map

    Notes
    -----
    The input ``rng`` argument is used to create a seed number, which
    is passed to the `numpy.random.SeedSequence` object to create a
    separate generator for each worker. The generator used for these
    is always created by a call to `numpy.random.default_rng`, even when
    the ``rng`` argument is None (indicating that the legacy random API is
    being used), or a different generator that that used by
    ``default_rng``.

    """
    if not callable(function):
        raise TypeError(f"input function '{repr(function)}' is not callable")

    # Check the function takes two arguments, with the second one
    # called rng.
    #
    sig = inspect.signature(function)
    if len(sig.parameters) != 2:
        raise TypeError(f"input function '{repr(function)}' does not take two arguments")

    names = list(sig.parameters)
    if names[1] != "rng":
        raise TypeError(f"input function '{repr(function)}' second argument is not called rng")

    if not np.iterable(sequence):
        raise TypeError(f"input '{repr(sequence)}' is not iterable")

    size = len(sequence)  # type: ignore[arg-type]

    if not _multi or size == 1 or (numcores is not None and numcores < 2):
        # As this is not in parallel the supplied generator can be
        # used.
        #
        debug("parallel_map_rng: running %d items in serial with rng=%s",
              size, rng)
        return [function(s, rng=rng) for s in sequence]

    # At this point we know context is not None but the typing code
    # does not.
    assert context is not None

    ncores = ncpus if numcores is None else numcores

    # Returns a started SyncManager object which can be used for sharing
    # objects between processes. The returned manager object corresponds
    # to a spawned child process and has methods which will create shared
    # objects and return corresponding proxies.
    manager = context.Manager()

    # Create FIFO queue and lock shared objects and return proxies to them.
    # The managers handles a server process that manages shared objects that
    # each slave process has access to.  Bottom line -- thread-safe.
    out_q = manager.Queue()
    err_q = manager.Queue()

    # if sequence is less than numcores, only use len sequence number of
    # processes
    if size < ncores:
        ncores = size

    # group sequence into numcores-worth of chunks
    sequence = split_array(sequence, ncores)

    debug("parallel_map_rng: running %d items in parallel (%d processes) with rng=%s",
          size, len(sequence), rng)

    # Create the RNG for each chunk. To be repeatable an explicit seed
    # is created given the input rng.
    #
    # TODO: Should the seed be created for a larger (i.e. more bits)
    # data type?
    #
    maxval = np.iinfo(np.uint64).max
    if rng is None:
        root_seed = np.random.randint(maxval, dtype=np.uint64)
    else:
        try:
            root_seed = rng.integers(maxval,  # type: ignore[union-attr]
                                     endpoint=False,
                                     dtype=np.uint64)
        except AttributeError:
            root_seed = rng.randint(maxval,  # type: ignore[union-attr]
                                    dtype=np.uint64)

    seeds = np.random.SeedSequence(root_seed).spawn(len(sequence))

    # This always uses default_rng; I can not see a sensible way to
    # allow the user to over-ride this (if they want to choose a
    # particular generator) or to get it to use np.random.RandomState
    # if required.
    #
    assert context.Process is not None
    procs = [context.Process(target=worker_rng,
                             args=(function, idx, chunk, out_q, err_q,
                                   np.random.default_rng(seed)))
             for idx, (chunk, seed) in enumerate(zip(sequence, seeds))]

    return run_tasks(procs, err_q, out_q)
