#!/usr/bin/env python3
# config.py
"""This module is used for configuration for my_io.py and get_gene_info.py scripts"""

# Error" doesn't conform to snake_case naming style
# pylint: disable=invalid-name


_UNIGENE_DIR = "/data/PROGRAMMING/assignment5"
_UNIGENE_FILE_ENDING = "unigene"


def get_unigene_directory():
    """returns file directory path"""
    return _UNIGENE_DIR


def get_unigene_extension():
    """returns file suffix"""
    return _UNIGENE_FILE_ENDING


def get_host_keywords():
    """Creates a dictionary for mapping common names to scientific names"""
    bos_tarus = "Bos_taurus"
    homo_sapiens = "Homo_sapiens"
    equus_caballus = "Equus_caballus"
    mus_musculus = "Mus_musculus"
    ovis_aries = "Ovis_aries"
    rattus_norvegicus = "Rattus_norvegicus"
    host_keywords = {
        "bos taurus": bos_tarus,
        "cow": bos_tarus,
        "cows": bos_tarus,
        "homo sapiens": homo_sapiens,
        "homo_sapiens": homo_sapiens,
        "human": homo_sapiens,
        "humans": homo_sapiens,
        "equus caballus": equus_caballus,
        "horse": equus_caballus,
        "horses": equus_caballus,
        "mus_musculus": mus_musculus,
        "mouse": mus_musculus,
        "mice": mus_musculus,
        "ovis_aries": ovis_aries,
        "sheep": ovis_aries,
        "sheeps": ovis_aries,
        "rattus_norvegicus": rattus_norvegicus,
        "rat": rattus_norvegicus,
        "rats": rattus_norvegicus
    }
    return host_keywords


def get_error_string_4_unable_to_open(file, mode):
    """ Print the invalid argument type message and exits the program """
    # error when used my_io.get_fh(file, "1234")
    print(f"Could not open the file: {file} for type '{mode}\n")


def get_error_string_4_ValueError():
    """ Print the invalid argument type message and exits the program """
    # error when used get_fh("file", "1234")
    print("Invalid argument Value for opening a file for reading/writing\n")


def get_error_string_4_TypeError():
    """ Print the invalid argument type message and exits the program """
    # error when used get_fh(file, "r", "w")
    print("Invalid argument Type passed in:\n")
