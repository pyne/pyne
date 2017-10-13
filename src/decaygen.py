#! /usr/bin/env python
"""This file generates a static C++ decayer function for use with PyNE.
It is suppossed to be fast.
"""
import os
import io
import sys
import pdb
import warnings
import traceback
from argparse import ArgumentParser, Namespace

import numpy as np
warnings.simplefilter('ignore', RuntimeWarning)
import tables as tb
import jinja2

from pyne.utils import QAWarning, toggle_warnings
warnings.simplefilter('ignore', QAWarning)
toggle_warnings()
from pyne import nuc_data
from pyne import nucname
from pyne.data import branch_ratio, half_life, decay_const, \
    decay_children, fpyield

ENV = jinja2.Environment(undefined=jinja2.StrictUndefined)

autogenwarn = """
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//              WARNING
// This file has been auto generated
// Do not modify directly. You have
// been warned. This is that warning
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
""".strip()

HEADER = ENV.from_string("""
#{{ dummy_ifdef }} PYNE_DECAY_IS_DUMMY
#ifndef PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM
#define PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM

{{ autogenwarn }}

// This file was generated with the following command:
// {{ args }}

#include <map>
//#include <cmath>

#ifndef PYNE_IS_AMALGAMATED
#include "data.h"
#include "nucname.h"
#endif

namespace pyne {
namespace decayers {

extern const int all_nucs[{{ nucs|length }}];

std::map<int, double> decay(std::map<int, double> comp, double t);

}  // namespace decayers
}  // namespace pyne

#endif  // PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM
#endif  // PYNE_DECAY_IS_DUMMY
""".strip())

SOURCE = ENV.from_string("""
#{{ dummy_ifdef }} PYNE_DECAY_IS_DUMMY
{{ autogenwarn }}

// This file was generated with the following command:
// {{ args }}

#ifndef PYNE_IS_AMALGAMATED
#include "decay.h"
#include "nucname.h"
#endif

namespace pyne {
namespace decayers {

{{ funcs }}

std::map<int, double> decay(std::map<int, double> comp, double t) {
  // setup
  using std::map;
  int nuc;
  int i = 0;
  double out [{{ nucs|length }}] = {};  // init to zero
  map<int, double> outcomp;

  // body
  map<int, double>::const_iterator it = comp.begin();
  for (; it != comp.end(); ++it) {
    switch (nucname::znum(it->first)) {
      {{ cases|indent(6) }}
      default:
        outcomp.insert(*it);
        break;
    }
  }

  // cleanup
  for (i = 0; i < {{ nucs|length }}; ++i)
    if (out[i] > 0.0)
      outcomp[nucname::state_id_to_id(all_nucs[i])] = out[i];
  return outcomp;
}

const int all_nucs [{{ nucs|length }}] = {
{{ nucs | join(", ") | wordwrap(width=78, break_long_words=False) | indent(2, True) }}
};

}  // namespace decayers
}  // namespace pyne

#endif  // PYNE_DECAY_IS_DUMMY
""".strip())


ELEM_FUNC = ENV.from_string("""
void decay_{{ elem|lower }}(double t, std::map<int, double>::const_iterator &it, std::map<int, double> &outcomp, double (&out)[{{ nucs|length }}]) {
  //using std::exp2;
  switch (nucname::id_to_state_id(it->first)) {
    {{ cases|indent(4) }}
    } default: {
      outcomp.insert(*it);
      break;
    }
  }
}
""".strip())


# Some strings that need not be redefined
BREAK = '  break;'
CHAIN_STMT = '  out[{0}] += {1};'
CHAIN_EXPR = '(it->second) * ({0})'
EXP_EXPR = 'exp2({a:.17e}*t)'
KEXP_EXPR = '{k:.17e}*' + EXP_EXPR
B_STMT = 'double b{b} = {exp};'
B_EXPR = 'b{b}'
KB_EXPR = '{k:.17e}*' + B_EXPR


def genfiles(nucs, short=1e-16, sf=False, dummy=False):
    ctx = Namespace(
        nucs=nucs,
        autogenwarn=autogenwarn,
        dummy_ifdef=('ifdef' if dummy else 'ifndef'),
        args=' '.join(sys.argv)
        )
    ctx.cases = gencases(nucs)
    ctx.funcs = genelemfuncs(nucs, short=short, sf=sf)
    hdr = HEADER.render(ctx.__dict__)
    src = SOURCE.render(ctx.__dict__)
    return hdr, src


def genchains(chains, sf=False):
    chain = chains[-1]
    children = decay_children(chain[-1])
    # filters spontaneous fission
    if not sf:
        children = {c for c in children if (0.0 == fpyield(chain[-1], c)) and (c not in chain) }
    if decay_const(chain[-1]) != 0:
        for child in children:
            if child not in chain:
                chains.append(chain + (child,))
                chains = genchains(chains, sf=sf)
    return chains


def k_from_hl_stable(hl, gamma):
    C = len(hl)
    outer = 1 / (hl[:C-1] - hl[:C-1, np.newaxis])
    # identity is ignored, set to unity
    mask = np.ones(C-1, dtype=bool)
    outer[mask, mask] = 1.0
    # end nuclide is stable so ignore
    # collapse by taking the product
    p = outer.prod(axis=0)
    k = -gamma * p * hl[:-1]**(C-2)
    k = np.append(k, gamma)
    return k


def k_from_hl_unstable(hl, gamma):
    outer = 1 / (hl - hl[:, np.newaxis])
    # identity is ignored, set to unity
    mask = np.ones(len(hl), dtype=bool)
    outer[mask, mask] = 1.0
    # collapse by taking the product
    p = outer.prod(axis=0)
    # get the other pieces
    C = len(hl)
    T_C = hl[-1]
    T_i_C = hl**(C - 2)
    # compute k
    k = (gamma * T_C) * T_i_C * p
    return k


def k_filter(k, short=1e-16):
    k_not_inf_or_nan = np.bitwise_and(~np.isinf(k), ~np.isnan(k))
    if k_not_inf_or_nan.sum() == 0:
        return k_not_inf_or_nan
    k_abs = np.abs(k)
    k_max = k_abs[k_not_inf_or_nan].max()
    k_filt = (k_abs / k_max) > short
    k_filt = np.bitwise_and(k_filt, k_not_inf_or_nan)
    return k_filt


def hl_filter(hl, short=1e-16):
    ends_stable = np.isinf(hl[-1])
    if ends_stable:
        hl_filt = (hl[:-1] / hl[:-1].sum()) > short
        hl_filt = np.append(hl_filt, True)
    else:
        hl_filt = (hl / hl.sum()) > short
    return hl_filt


def k_a_from_hl(chain, short=1e-16):
    hl = np.array([half_life(n, False) for n in chain])
    hl = hl[~np.isnan(hl)]
    a = -1.0 / hl
    gamma = np.prod([branch_ratio(p, c) for p, c in zip(chain[:-1], chain[1:])])
    if gamma == 0.0 or np.isnan(gamma):
        return None, None
    ends_stable = np.isinf(hl[-1])
    if ends_stable:
        k = k_from_hl_stable(hl, gamma)
    else:
        k = k_from_hl_unstable(hl, gamma)
    # filtering makes compiling faster by pre-ignoring negligible species
    # in this chain. They'll still be picked up in their own chains.
    mask = k_filter(k, short=short)
    if mask.sum() == 0:
        return None, None
    return k[mask], a[mask]


def kbexpr(k, b):
    if k == 1.0:
        return B_EXPR.format(b=b)
    else:
        return KB_EXPR.format(k=k, b=b)


def ensure_cse(a_i, b, cse):
    bkey = EXP_EXPR.format(a=a_i)
    if bkey not in cse:
        b += 1
        cse[bkey] = b
    return b

def b_from_a(cse, a_i):
    bkey = EXP_EXPR.format(a=a_i)
    return cse[bkey]


def chainexpr(chain, cse, b, bt, short=1e-16):
    child = chain[-1]
    if len(chain) == 1:
        a_i = -1.0 / half_life(child, False)
        b = ensure_cse(a_i, b, cse)
        terms = B_EXPR.format(b=b_from_a(cse, a_i))
    else:
        k, a = k_a_from_hl(chain, short=short)
        if k is None:
            return None, b, bt
        terms = []
        for k_i, a_i in zip(k, a):
            if k_i == 1.0 and a_i == 0.0:
                term = str(1.0 - bt)  # a slight optimization
                bt = 1
            elif a_i == 0.0:
                if not np.isnan(k_i):
                    if bt < 1:
                        if k_i + bt < 1:
                            term = '{0:.17e}'.format(k_i)  # another slight optimization
                            bt += k_i
                        else:
                            term = '{0:.17e}'.format(1.0 - bt)
                            bt = 1.0
                    else:
                        term = '0'
                else:
                    term = '0'
            else:
                b = ensure_cse(a_i, b, cse)
                term = kbexpr(k_i, b_from_a(cse, a_i))
            terms.append(term)
        terms = ' + '.join(terms)
    return CHAIN_EXPR.format(terms), b, bt


def gencase(nuc, idx, b, short=1e-16, sf=False):
    case = ['}} case {0}: {{'.format(nuc)]
    dc = decay_const(nuc, False)
    if dc == 0.0:
        # stable nuclide
        case.append(CHAIN_STMT.format(idx[nuc], 'it->second'))
    else:
        chains = genchains([(nuc,)], sf=sf)
        print('{} has {} chains'.format(nucname.name(nuc), len(set(chains))))
        cse = {}  # common sub-expression exponents to elimnate
        bt = 0
        for c in chains:
            if c[-1] not in idx:
                continue
            cexpr, b, bt = chainexpr(c, cse, b, bt, short=short)
            if cexpr is None:
                continue
            case.append(CHAIN_STMT.format(idx[c[-1]], cexpr))
        bstmts = ['  ' + B_STMT.format(exp=exp, b=bval) for exp, bval in \
                  sorted(cse.items(), key=lambda x: x[1])]
        case = case[:1] + bstmts + case[1:]
    case.append(BREAK)
    return case, b


def elems(nucs):
    return sorted(set(map(nucname.znum, nucs)))


def gencases(nucs):
    switches = []
    for i in elems(nucs):
        c = ['case {0}:'.format(i),
             '  decay_{0}(t, it, outcomp, out);'.format(nucname.name(i).lower()),
             '  break;']
        switches.append('\n'.join(c))
    return '\n'.join(switches)


def genelemfuncs(nucs, short=1e-16, sf=False):
    idx = dict(zip(nucs, range(len(nucs))))
    cases = {i: [-1, []] for i in elems(nucs)}
    for nuc in nucs:
        z = nucname.znum(nuc)
        case, cases[z][0] = gencase(nuc, idx, cases[z][0], short=short, sf=sf)
        cases[z][1] += case
    funcs = []
    for i, (b, kases) in cases.items():
        kases[0] = kases[0][2:]
        ctx = dict(nucs=nucs, elem=nucname.name(i), cases='\n'.join(kases))
        funcs.append(ELEM_FUNC.render(ctx))
    return "\n\n".join(funcs)


def load_default_nucs():
    with tb.open_file(nuc_data) as f:
        ll = f.root.decay.level_list
        stable = ll.read_where('(nuc_id%10000 == 0) & (nuc_id != 0)')
        metastable = ll.read_where('metastable > 0')
    nucs = set(int(nuc) for nuc in stable['nuc_id'])
    nucs |= set(int(nuc) for nuc in metastable['nuc_id'])
    nucs = sorted(nuc for nuc in nucs if not np.isnan(decay_const(nuc, False)))
    return nucs




def build_tarfile(ns):
    import tarfile
    with tarfile.open('decay.tar.gz', 'w:gz') as tar:
        tar.add(ns.hdr)
        tar.add(ns.src)


def write_if_diff(filename, contents):
    """Only writes the file if it is different. This prevents touching the file needlessly."""
    if not os.path.isfile(filename):
        existing = None
    else:
        with io.open(filename, 'r') as f:
            existing = f.read()
    if contents == existing:
        return
    with io.open(filename, 'w') as f:
        f.write(contents)


def build(hdr='decay.h', src='decay.cpp', nucs=None, short=1e-16, sf=False,
          dummy=False):
    nucs = load_default_nucs() if nucs is None else list(map(nucname.id, nucs))
    #nucs = nucs[:200]
    h, s = genfiles(nucs, short=short, sf=sf, dummy=dummy)
    write_if_diff(hdr, h)
    write_if_diff(src, s)


def main():
    parser = ArgumentParser('decay-gen')
    parser.add_argument('--hdr', default='decay.h', help='The header file name.')
    parser.add_argument('--src', default='decay.cpp', help='The source file name.')
    parser.add_argument('--nucs', nargs='+', default=None,
                        help='Nuclides to generate for.')
    parser.add_argument('--dummy', action='store_true', default=False,
                        dest='dummy', help='Makes dummy versions as '
                        'compile-time fallbacks.')
    parser.add_argument('--no-dummy', action='store_false', default=False,
                        dest='dummy', help='Makes regular files.')
    parser.add_argument('--filter-short', default=1e-16, type=float, dest='short',
                        help='Fraction of sum of all half-lives below which a '
                             'nuclide is filtered from a decay chain, default 1e-16.')
    parser.add_argument('--spontaneous-fission', default=False, action='store_true',
                        dest='sf', help='Includes spontaneous fission decay chains, '
                                        'default False.')
    parser.add_argument('--tar', action='store_true', default=False,
                        help='Builds decay.tar.gz')
    parser.add_argument('--cred', default='../rs.cred',
                        help='Path to credentials file.')
    parser.add_argument('--no-build', dest='build', default=True, action='store_false',
                       help='Does not build the source code.')
    ns = parser.parse_args()
    if ns.build:
        try:
            build(hdr=ns.hdr, src=ns.src, nucs=ns.nucs, short=ns.short, sf=ns.sf,
                  dummy=ns.dummy)
        except Exception:
            type, value, tb = sys.exc_info()
            traceback.print_exc()
            pdb.post_mortem(tb)

    if ns.tar:
        print("building decay.tar.gz ...")
        build_tarfile(ns)


if __name__ == '__main__':
    main()
