#! /usr/bin/env python
"""This file generates a static C++ decayer function for use with PyNE.
It is suppossed to be fast.
"""
import io
import warnings
from argparse import ArgumentParser, Namespace

import numpy as np
import tables as tb
import jinja2

from pyne.utils import QAWarning
warnings.simplefilter('ignore', QAWarning)
from pyne import nuc_data
from pyne import data
from pyne import rxname
from pyne import nucname
from pyne.data import branch_ratio, half_life, decay_const, decay_children
from pyne.material import Material

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
#ifndef PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM
#define PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM

{{ autogenwarn }}

#include <map>

#include "data.h"

namespace pyne {
namespace decayers {

extern const int all_nucs[{{ nucs|length }}];

std::map<int, double> decay(std::map<int, double> comp, double t);

}  // namespace decayers
}  // namespace pyne

#endif  // PYNE_GEUP5PGEJBFGNHGI36TRBB4WGM
""".strip())

SOURCE = ENV.from_string("""
{{ autogenwarn }}

#include "decay.h"

namespace pyne {
namespace decayers {

std::map<int, double> decay(std::map<int, double> comp, double t) {
  // setup
  using std::map;
  int nuc;
  double out [{{ nucs|length }}] = {};  // init to zero
  map<int, double> outcomp;
  
  // body
  map<int, double>::const_iterator it = comp.begin();
  for (; it != comp.end(); ++it) {
    switch (it->first) {
      {{ cases|indent(6) }}
      } default: {
        outcomp.insert(*it);
        break;
      }
    }
  }
  
  // cleanup
  for (int i = 0; i < {{ nucs|length }}; ++i)
    if (out[i] > 0.0)
      outcomp[all_nucs[i]] = out[i];
  return outcomp;
}

const int all_nucs [{{ nucs|length }}] = {
{{ nucs | join(", ") | wordwrap(width=78, break_long_words=False) | indent(2, True) }}
};

}  // namespace decayers
}  // namespace pyne
""".strip())


# In[44]:

# Some strings that need not be redefined
BREAK = '  break;'
CHAIN_STMT = '  out[{0}] += {1};'
CHAIN_EXPR = '(it->second) * ({0})'
EXP_EXPR = 'exp2({a}*t)'
KEXP_EXPR = '{k}*' + EXP_EXPR

def genfiles(nucs):
    ctx = Namespace(
        nucs=nucs,
        autogenwarn=autogenwarn,
        )
    ctx.cases = gencases(nucs)
    hdr = HEADER.render(ctx.__dict__)
    src = SOURCE.render(ctx.__dict__)
    return hdr, src


def genchains(chains):
    chain = chains[-1]
    children = decay_children(chain[-1])
    for child in children:
        chains.append(chain + (child,))
        chains = genchains(chains)
    return chains


def k_a(chain):
    a = [-1.0 / half_life(n) for n in chain]
    dc = np.array(list(map(decay_const, chain)))
    if np.isnan(dc).any():
        return None, None
    cij = dc[:, np.newaxis] / (dc[:, np.newaxis] - dc)
    mask = np.ones(len(chain), dtype=bool)
    cij[mask, mask] = 1.0
    ci = cij.prod(axis=0)
    k = (dc / dc[-1]) * ci
    return k, a


def kexpexpr(k, a):
    if k == 1.0:
        return EXP_EXPR.format(a=a)
    else:
        return KEXP_EXPR.format(k=k, a=a)


def chainexpr(chain):
    child = chain[-1]
    dc_child = decay_const(child)
    if len(chain) == 1:
        a = -1.0 / half_life(child)
        terms = EXP_EXPR.format(a=a)
    elif dc_child == 0.0:
        k, a = k_a(chain[:-1])
        if k is None:
            return None
        terms = ['1.0']
        terms += [kexpexpr(k_i, a_i) for k_i, a_i in zip(k, a)]
        terms = ' - '.join(terms)
    else:
        k, a = k_a(chain)
        if k is None:
            return None
        terms = [kexpexpr(k_i, a_i) for k_i, a_i in zip(k, a)]
        terms = ' + '.join(terms)
    return CHAIN_EXPR.format(terms)


def gencase(nuc, idx):
    case = ['}} case {0}: {{'.format(nuc)]
    dc = decay_const(nuc)
    if dc == 0.0:
        # stable nuclide
        case.append(CHAIN_STMT.format(idx[nuc], 'it->second'))
    else:
        chains = genchains([(nuc,)])
        for c in chains:
            if c[-1] not in idx:
                continue
            cexpr = chainexpr(c)
            if cexpr is None:
                continue
            case.append(CHAIN_STMT.format(idx[c[-1]], cexpr))
    case.append(BREAK)
    return case


def gencases(nucs):
    idx = dict(zip(nucs, range(len(nucs))))
    cases = []
    for nuc in nucs:
        cases += gencase(nuc, idx)
    cases[0] = cases[0][2:]
    return '\n'.join(cases)


def load_default_nucs():
    with tb.open_file(nuc_data) as f:
        ll = f.root.decay.level_list
        stable = ll.read_where('(nuc_id%10000 == 0) & (nuc_id != 0)')
        metastable = ll.read_where('metastable > 0')
    nucs = set(int(nuc) for nuc in stable['nuc_id']) 
    nucs |= set(nucname.state_id_to_id(int(nuc)) for nuc in metastable['nuc_id']) 
    nucs = sorted(nuc for nuc in nucs if not np.isnan(decay_const(nuc)))
    return nucs


def main():
    parser = ArgumentParser('decay-gen')
    parser.add_argument('--hdr', default='decay.h', help='The header file name.')
    parser.add_argument('--src', default='decay.cpp', help='The source file name.')
    parser.add_argument('--nucs', nargs='+', default=None, 
                        help='Nuclides to generate for.')
    ns = parser.parse_args()
    nucs = load_default_nucs() if ns.nucs is None else list(map(nucname.id, ns.nucs))
    hdr, src = genfiles(nucs)
    with io.open(ns.hdr, 'w') as f:
        f.write(hdr)
    with io.open(ns.src, 'w') as f:
        f.write(src)


if __name__ == '__main__':
    main()

