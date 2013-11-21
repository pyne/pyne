"""
Tools to generate, set and check the hashes of datasets in pyne. Specifically,
these are the high level tools to iterate over all lower level functions in
hasher

Author: crbates
"""

import ConfigParser
import os
import hashlib

import pyne.dbgen.hasher as hasher


def write_hash_config(nuc_data):
    config = ConfigParser.ConfigParser()
    config.add_section('Hashes')
    mhash = hashlib.md5()
    for item in dir(hasher):
        if "calc" in item:
            fn = getattr(hasher, item)
            hashv = fn(nuc_data)
            config.set('Hashes', item[5:-5], hashv)
            mhash.update(hashv)
    config.set('Hashes', 'global', mhash.hexdigest())
    local_filename = os.path.join(os.path.split(__file__)[0], 'nuc_hash.cfg')
    with open(local_filename, 'wb') as configfile:
        config.write(configfile)


def check_hashes(nuc_data):
    check_list = []
    config = ConfigParser.ConfigParser()
    local_filename = os.path.join(os.path.split(__file__)[0], 'nuc_hash.cfg')
    config.read(local_filename)
    mhash = hashlib.md5()
    for hashopt in config.options('Hashes'):
        for item in dir(hasher):
            if "calc_" + hashopt in item:
                fn = getattr(hasher, item)
                hashv = fn(nuc_data)
                val = (config.get('Hashes', hashopt) == hashv)
                check_list.append((hashopt, val))
                mhash.update(hashv)
    val = (config.get('Hashes', 'global') == mhash.hexdigest())
    check_list.append(("all", val))
    return check_list


def set_internal_hashes(nuc_data):
    for item in dir(hasher):
        if "set" in item:
            fn = getattr(hasher, item)
            fn(nuc_data)


def check_internal_hashes(nuc_data):
    for item in dir(hasher):
        if "check" in item:
            fn = getattr(hasher, item)
            hashv = fn(nuc_data)
            check_list.append((item[5:-5], hashv))