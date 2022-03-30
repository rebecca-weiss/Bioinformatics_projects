#!/usr/bin/env python3
# my_io.py
"""my_io module file"""

def get_fh(filename, action_rw):
    """Open the file name passed in, and passes back a a file object, i.e. the file handle"""
    try:
        fhandle = open(filename, action_rw)
    except IOError:
        raise IOError(f"{filename} does not exist in this directory")

    except ValueError:
        if action_rw != 'r' or action_rw != 'w':
            raise ValueError(f"Must indicate 'r' to read or 'w' to write; {action_rw} not accepted")

    return fhandle



