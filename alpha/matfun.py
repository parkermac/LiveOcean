"""
Functions cribbed from J Paul Rinheimer for reading in matlab files.
"""

import scipy.io as spio
import numpy as np

def print_mat_nested(d, indent=0, nkeys=0):
    """
    Pretty print nested structures from .mat files
    Inspired by: `StackOverflow <http://stackoverflow.com/questions/ ...
    3229419/pretty-printing-nested-dictionaries-in-python>`_
    """

    # Subset dictionary to limit keys to print.  Only works on first level
    if nkeys>0:
        d = {k: d[k] for k in d.keys()[:nkeys]}
        # Dictionary comprehension: limit to first nkeys keys.

    if isinstance(d, dict):
        for key, value in d.iteritems():
        # iteritems loops through key, value pairs
            print('\t' * indent + 'Key: ' + str(key))
            print_mat_nested(value, indent+1)

    if isinstance(d,np.ndarray) and d.dtype.names is not None:
        # Note: and short-circuits by default
        for n in d.dtype.names:
            # This means it's a struct, it's bit of a kludge test.
            print('\t' * indent + 'Field: ' + str(n))
            print_mat_nested(d[n], indent+1)


def loadmat(filename):
    """
    This function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects

    from: `StackOverflow <http://stackoverflow.com/questions/7008608/ ...
    scipy-io-loadmat-nested-structures-i-e-dictionaries>`_
    """
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    """
    Checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries.
    """
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries.
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict