"""This is a command line interface for manipulating ORIGEN v2.2 TAPE9.INP files.
"""
from __future__ import print_function
import os
import sys
import argparse
from glob import glob
from pyne.utils import QA_warn

import numpy as np

from pyne import data
from pyne import ensdf
from pyne import nucname
from pyne import origen22
from pyne.xs.cache import XSCache
from pyne.xs.data_source import EAFDataSource, SimpleDataSource, NullDataSource
from pyne.dbgen.api import build_dir
from pyne.dbgen.decay import grab_ensdf_decay

QA_warn(__name__)


def parse_ensdf(files):
    """Parses a list of ensdf files for origen."""
    decays = []
    branches = []
    for f in files:
        decs, brs = ensdf.origen_data(f)
        decays.extend(decs)
        branches.extend(brs)
    decays = [x for x in decays if x[0] is not None]
    branches = [x for x in branches if x[0] is not None]
    return decays, branches


def _is_metastable_beta_decay_0(item, metastable_cutoff):
    return (
        (item[3] is not None)
        and (item[3] > metastable_cutoff)
        and (item[4] > 0)
        and (item[3] != np.inf)
        and (item[2] == 0)
    )


def _is_metastable_beta_decay_x(item, metastable_cutoff):
    return (
        (item[3] is not None)
        and (item[3] > metastable_cutoff)
        and (item[4] > 0)
        and (item[3] != np.inf)
        and (item[2] > 0)
    )


def _plus_eq_lib(lib, field, key, val):
    lib[field][key] = lib[field].get(key, 0.0) + val


def _plus_eq_decay_t9(t9, field, nuc, key, val):
    if nuc in origen22.ACTIVATION_PRODUCT_NUCS:
        _plus_eq_lib(t9[1], field, key, val)
    if nuc in origen22.ACTINIDE_AND_DAUGHTER_NUCS:
        _plus_eq_lib(t9[2], field, key, val)
    if nuc in origen22.FISSION_PRODUCT_NUCS:
        _plus_eq_lib(t9[3], field, key, val)

    if nuc not in origen22.NUCS:
        # guess its classification
        if nucname.anum in nucname.act:
            # nuc is actinide
            _plus_eq_lib(t9[2], field, key, val)
        else:
            # default to activation product
            _plus_eq_lib(t9[1], field, key, val)


def _eq_decay_t9(t9, field, nuc, key, val):
    if nuc in origen22.ACTIVATION_PRODUCT_NUCS:
        t9[1][field][key] = val
    if nuc in origen22.ACTINIDE_AND_DAUGHTER_NUCS:
        t9[2][field][key] = val
    if nuc in origen22.FISSION_PRODUCT_NUCS:
        t9[3][field][key] = val

    if nuc not in origen22.NUCS:
        # guess its classification
        if nucname.anum in nucname.act:
            # nuc is actinide
            t9[2][field][key] = val
        else:
            # default to activation product
            t9[1][field][key] = val


def _set_branch_item(t9, nuc, key, item):
    _eq_decay_t9(t9, "half_life", nuc, key, item[2] or 0.0)
    if "%SF" in item[3] and item[3]["%SF"] != "?":
        _eq_decay_t9(
            t9, "frac_spont_fiss", nuc, key, float(item[3]["%SF"] or 0.0) / 100.0
        )
    if "%EC" in item[3] and item[3]["%EC"] != "?":
        _plus_eq_decay_t9(
            t9,
            "frac_beta_plus_or_electron_capture",
            nuc,
            key,
            float(item[3]["%EC"] or 0.0) / 100.0,
        )
    if "%B+" in item[3] and item[3]["%B+"] != "?":
        _plus_eq_decay_t9(
            t9,
            "frac_beta_plus_or_electron_capture",
            nuc,
            key,
            float(item[3]["%B+"] or 0.0) / 100.0,
        )
    if "%EC+%B+" in item[3] and item[3]["%EC+%B+"] != "?":
        _plus_eq_decay_t9(
            t9,
            "frac_beta_plus_or_electron_capture",
            nuc,
            key,
            float(item[3]["%EC+%B+"] or 0.0) / 100.0,
        )
    if "%B-N" in item[3] and item[3]["%B-N"] != "?":
        _eq_decay_t9(t9, "frac_beta_n", nuc, key, float(item[3]["%B-N"] or 0.0) / 100.0)
    if "%A" in item[3] and item[3]["%A"] != "?":
        _eq_decay_t9(t9, "frac_alpha", nuc, key, float(item[3]["%A"] or 0.0) / 100.0)
    if "%IT" in item[3] and item[3]["%IT"] != "?":
        _eq_decay_t9(
            t9,
            "frac_isomeric_transition",
            nuc,
            key,
            float(item[3]["%IT"] or 0.0) / 100.0,
        )


def gendecay(decays, branches, metastable_cutoff=1.0):
    """This computes ORIGEN TAPE9 decay data based on ENSDF data.

    Parameters
    ----------
    decays : list
        decay list from parse_ensdf()
    branches : list
        branches list from parse_ensdf()
    metastable_cutoff : float, optional
        minimum lifetime of metastable state (in seconds) to be included.

    Returns
    -------
    t9 : dict
        a TAPE9 dictionary for the decay library
    """
    t9 = {
        1: {"_type": "decay", "title": "PyNE Decay Data for Activation Products"},
        2: {"_type": "decay", "title": "PyNE Decay Data for Actinides & Daughters"},
        3: {"_type": "decay", "title": "PyNE Decay Data for Fission Products"},
    }
    for nlb, lib in t9.items():
        for field in origen22.DECAY_FIELDS:
            lib[field] = {}

    longest = {}
    longest2 = {}
    for item in decays:
        nuc = nucname.id(item[0])
        key = nucname.zzaaam(nuc)
        if _is_metastable_beta_decay_0(item, metastable_cutoff):
            if "B-" in item[5]:
                _plus_eq_decay_t9(t9, "frac_beta_minus_x", nuc, key, item[6] / 100.0)
            if "B+" in item[5] or "EC" in item[5]:
                _plus_eq_decay_t9(
                    t9,
                    "frac_beta_plus_or_electron_capture_x",
                    nuc,
                    key,
                    item[6] / 100.0,
                )
        if _is_metastable_beta_decay_x(item, metastable_cutoff):
            key += 1
            longest2[key] = longest2.get(key, 0)
            if item[1] == longest2[key]:
                if "B-" in item[5]:
                    _plus_eq_decay_t9(
                        t9, "frac_beta_minus_x", nuc, key, item[6] / 100.0
                    )
                    # item[6]*item[8]/100.0)
                if "B+" in item[5] or "EC" in item[5]:
                    _plus_eq_decay_t9(
                        t9,
                        "frac_beta_plus_or_electron_capture_x",
                        nuc,
                        key,
                        item[6] / 100.0,
                    )
                    # key, item[6]*item[8]/100.0)
            elif item[1] > longest2[key]:
                longest2[key] = item[1]
                if "B-" in item[5]:
                    # _eq_decay_t9(t9, 'frac_beta_minus_x', nuc, key, item[6]*item[8]/100.0)
                    _eq_decay_t9(t9, "frac_beta_minus_x", nuc, key, item[6] / 100.0)
                if "B+" in item[5] or "EC" in item[5]:
                    _eq_decay_t9(
                        t9,
                        "frac_beta_plus_or_electron_capture_x",
                        nuc,
                        key,
                        item[6] / 100.0,
                    )
                    # key, item[6]*item[8]/100.0)
    for item in branches:
        nuc = nucname.id(item[0])
        key = nucname.zzaaam(nuc)
        if (item[1] == 0) and (item[2] > metastable_cutoff):
            _set_branch_item(t9, nuc, key, item)
        if (item[1] != 0) and (item[2] > metastable_cutoff):
            key += 1
            longest[key] = longest.get(key, 0)
            if item[2] <= longest[key]:
                continue
            _set_branch_item(t9, nuc, key, item)
    for nucs, hl in zip(
        [
            origen22.ACTIVATION_PRODUCT_NUCS,
            origen22.ACTINIDE_AND_DAUGHTER_NUCS,
            origen22.FISSION_PRODUCT_NUCS,
        ],
        [t9[i]["half_life"] for i in range(1, 4)],
    ):
        for nuc in nucs:
            key = nucname.zzaaam(nuc)
            if key not in hl:
                hl[key] = data.half_life(nuc)
    return t9


def main_gen(ns):
    """Generates an open TAPE9.INP file. by default this only uses completely open
    data.
    """
    files = glob(os.path.join(ns.build_dir, "ENSDF", "ensdf.*"))
    if len(files) == 0:
        grab_ensdf_decay(ns.build_dir)
        files = glob(os.path.join(ns.build_dir, "ENSDF", "ensdf.*"))
    print("parsing ENSDF decay data")
    decays, branches = parse_ensdf(files)
    print("creating ORIGEN decay libraries")
    t9decay = gendecay(decays, branches, metastable_cutoff=ns.metastable_cutoff)
    print("creating ORIGEN cross section libraries")
    xsc = XSCache(data_source_classes=[EAFDataSource, SimpleDataSource, NullDataSource])
    xsc.load()
    t9xsfpy = origen22.xslibs(xscache=xsc, verbose=True)
    t9 = origen22.merge_tape9([t9decay, t9xsfpy])
    origen22.write_tape9(t9, outfile=ns.filename)


_cmd_mains = {
    "gen": main_gen,
}


def main():
    parser = argparse.ArgumentParser(
        description="Manipulates ORIGEN v2.2 " "TAPE9.INP files."
    )
    subparsers = parser.add_subparsers(
        title="cmd",
        help="available sub-commands",
        description="the subcommands",
        dest="cmd",
    )

    gen = subparsers.add_parser(
        "gen", help="Creates a TAPE9 file based only on open data."
    )
    gen.add_argument("-o", dest="filename", default="TAPE9.INP", help="output filename")
    gen.add_argument(
        "-b",
        dest="build_dir",
        action="store",
        default=build_dir,
        help="path to the build directory.",
    )
    gen.add_argument(
        "--metastable-cutoff",
        dest="metastable_cutoff",
        default=1.0,
        type=float,
        help="minimum cutoff value to be considered " "metastable state, in sec.",
    )
    ns = parser.parse_args()
    if ns.cmd not in _cmd_mains:
        sys.exit("command {0!r} could not be found".format(ns.cmd))
    _cmd_mains[ns.cmd](ns)
    sys.exit()


if __name__ == "__main__":
    main()
