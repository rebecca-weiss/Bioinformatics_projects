#!/usr/bin/env python3
# my_io.py
"""my_io module file"""

import os
from assignment5 import config

def get_fh(file=None, mode=None):
    """Open the file name passed in, and passes back a file object, i.e. the filehandle"""
    try:
        fobj = open(file, mode)
    except IOError as err:
        config.get_error_string_4_unable_to_open(file, mode)
        raise err
    except ValueError as err: #test something like my_io.get_fh("does_not_exist.txt", "rrr")
        config.get_error_string_4_ValueError()
        raise err
    except TypeError as err: #test something like my_io.get_fh([], "r")
        config.get_error_string_4_TypeError()
        raise err
    return fobj


def is_valid_gene_file_name(file):
    """Check to make sure the given file name exists, if it does it return True (else returns False)"""
    return os.path.exists(file)
