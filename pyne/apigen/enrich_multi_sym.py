"""Generate the API for symbolic multicomponent enrichment using sympy.
Note that there is a bug in sympy v0.7.2 that prevents the cse() function
from being used with infinities.  For a work around see [1].

1. https://groups.google.com/forum/#!msg/sympy/YL1R_hR6OKQ/axKrCsCSMQsJ
"""
from __future__ import print_function, division
import os
import logging
import multiprocessing
import time
from pyne.utils import QA_warn

from sympy import (
    Symbol,
    pprint,
    latex,
    diff,
    count_ops,
    simplify,
    cse,
    Eq,
    Q,
    log,
    logcombine,
    Abs,
    exp,
    sqrt,
    series,
    separate,
    powsimp,
    collect,
    expand,
    Abs,
)
from sympy.solvers import solve
from sympy.utilities.iterables import numbered_symbols

from utils import cse_to_c

QA_warn(__name__)

NPROCS = 10


def _aggstatus(stat, msg, aggstat):
    if not aggstat:
        print(msg)
    stat += msg + "\n"
    return stat


def cgen_ncomp(ncomp=3, nporder=2, aggstat=False, debug=False):
    """Generates a C function for ncomp (int) number of components.
    The jth key component is always in the first position and the kth
    key component is always in the second.  The number of enrichment
    stages (NP) is calculated via a taylor series approximation.  The
    order of this approximation may be set with nporder.  Only values
    of 1 or 2 are allowed. The aggstat argument determines whether the
    status messages should be aggreated and printed at the end or output
    as the function executes.
    """
    start_time = time.time()
    stat = _aggstatus("", "generating {0} component enrichment".format(ncomp), aggstat)
    r = range(0, ncomp)
    j = 0
    k = 1

    # setup-symbols
    alpha = Symbol("alpha", positive=True, real=True)
    LpF = Symbol("LpF", positive=True, real=True)
    PpF = Symbol("PpF", positive=True, real=True)
    TpF = Symbol("TpF", positive=True, real=True)
    SWUpF = Symbol("SWUpF", positive=True, real=True)
    SWUpP = Symbol("SWUpP", positive=True, real=True)
    NP = Symbol("NP", positive=True, real=True)  # Enrichment Stages
    NT = Symbol("NT", positive=True, real=True)  # De-enrichment Stages
    NP0 = Symbol("NP0", positive=True, real=True)  # Enrichment Stages Initial Guess
    NT0 = Symbol("NT0", positive=True, real=True)  # De-enrichment Stages Initial Guess
    NP1 = Symbol("NP1", positive=True, real=True)  # Enrichment Stages Computed Value
    NT1 = Symbol("NT1", positive=True, real=True)  # De-enrichment Stages Computed Value
    Mstar = Symbol("Mstar", positive=True, real=True)
    MW = [Symbol("MW[{0}]".format(i), positive=True, real=True) for i in r]
    beta = [alpha ** (Mstar - MWi) for MWi in MW]

    # np_closed helper terms
    NP_b = Symbol("NP_b", real=True)
    NP_2a = Symbol("NP_2a", real=True)
    NP_sqrt_base = Symbol("NP_sqrt_base", real=True)

    xF = [Symbol("xF[{0}]".format(i), positive=True, real=True) for i in r]
    xPi = [Symbol("xP[{0}]".format(i), positive=True, real=True) for i in r]
    xTi = [Symbol("xT[{0}]".format(i), positive=True, real=True) for i in r]
    xPj = Symbol("xPj", positive=True, real=True)
    xFj = xF[j]
    xTj = Symbol("xTj", positive=True, real=True)
    ppf = (xFj - xTj) / (xPj - xTj)
    tpf = (xFj - xPj) / (xTj - xPj)

    xP = [
        (
            ((xF[i] / ppf) * (beta[i] ** (NT + 1) - 1))
            / (beta[i] ** (NT + 1) - beta[i] ** (-NP))
        )
        for i in r
    ]
    xT = [
        (
            ((xF[i] / tpf) * (1 - beta[i] ** (-NP)))
            / (beta[i] ** (NT + 1) - beta[i] ** (-NP))
        )
        for i in r
    ]
    rfeed = xFj / xF[k]
    rprod = xPj / xP[k]
    rtail = xTj / xT[k]

    # setup constraint equations
    numer = [
        ppf * xP[i] * log(rprod) + tpf * xT[i] * log(rtail) - xF[i] * log(rfeed)
        for i in r
    ]
    denom = [log(beta[j]) * ((beta[i] - 1.0) / (beta[i] + 1.0)) for i in r]
    LoverF = sum([n / d for n, d in zip(numer, denom)])
    SWUoverF = -1.0 * sum(numer)
    SWUoverP = SWUoverF / ppf

    prod_constraint = (xPj / xFj) * ppf - (beta[j] ** (NT + 1) - 1) / (
        beta[j] ** (NT + 1) - beta[j] ** (-NP)
    )
    tail_constraint = (xTj / xFj) * (sum(xT)) - (1 - beta[j] ** (-NP)) / (
        beta[j] ** (NT + 1) - beta[j] ** (-NP)
    )
    # xp_constraint = 1.0 - sum(xP)
    # xf_constraint = 1.0 - sum(xF)
    # xt_constraint = 1.0 - sum(xT)

    # This is NT(NP,...) and is correct!
    # nt_closed = solve(prod_constraint, NT)[0]

    # However, this is NT(NP,...) rewritten (by hand) to minimize the number of NP
    # and M* instances in the expression.  Luckily this is only depends on the key
    # component and remains general no matter the number of components.
    nt_closed = (
        -MW[0] * log(alpha)
        + Mstar * log(alpha)
        + log(xTj)
        + log((-1.0 + xPj / xF[0]) / (xPj - xTj))
        - log(
            alpha ** (NP * (MW[0] - Mstar))
            * (xF[0] * xPj - xPj * xTj)
            / (-xF[0] * xPj + xF[0] * xTj)
            + 1
        )
    ) / ((MW[0] - Mstar) * log(alpha))

    # new expression for normalized flow rate
    # NOTE: not needed, solved below
    # loverf = LoverF.xreplace({NT: nt_closed})

    # Define the constraint equation with which to solve NP. This is chosen such to
    # minimize the number of ops in the derivatives (and thus np_closed).  Other,
    # more verbose possibilities are commented out.
    # np_constraint = (xP[j]/sum(xP) - xPj).xreplace({NT: nt_closed})
    # np_constraint = (xP[j]- sum(xP)*xPj).xreplace({NT: nt_closed})
    # np_constraint = (xT[j]/sum(xT) - xTj).xreplace({NT: nt_closed})
    np_constraint = (xT[j] - sum(xT) * xTj).xreplace({NT: nt_closed})

    # get closed form approximation of NP via symbolic derivatives
    stat = _aggstatus(stat, "  order-{0} NP approximation".format(nporder), aggstat)
    d0NP = np_constraint.xreplace({NP: NP0})
    d1NP = diff(np_constraint, NP, 1).xreplace({NP: NP0})
    if 1 == nporder:
        np_closed = NP0 - d1NP / d0NP
    elif 2 == nporder:
        d2NP = diff(np_constraint, NP, 2).xreplace({NP: NP0}) / 2.0
        # taylor series polynomial coefficients, grouped by order
        # f(x) = ax**2 + bx + c
        a = d2NP
        b = d1NP - 2 * NP0 * d2NP
        c = d0NP - NP0 * d1NP + NP0 * NP0 * d2NP
        # quadratic eq. (minus only)
        # np_closed = (-b - sqrt(b**2 - 4*a*c)) / (2*a)
        # However, we need to break up this expr as follows to prevent
        # a floating point arithmetic bug if b**2 - 4*a*c is very close
        # to zero but happens to be negative.  LAME!!!
        np_2a = 2 * a
        np_sqrt_base = b**2 - 4 * a * c
        np_closed = (-NP_b - sqrt(NP_sqrt_base)) / (NP_2a)
    else:
        raise ValueError("nporder must be 1 or 2")

    # generate cse for writing out
    msg = "  minimizing ops by eliminating common sub-expressions"
    stat = _aggstatus(stat, msg, aggstat)
    exprstages = [
        Eq(NP_b, b),
        Eq(NP_2a, np_2a),
        # fix for floating point sqrt() error
        Eq(NP_sqrt_base, np_sqrt_base),
        Eq(NP_sqrt_base, Abs(NP_sqrt_base)),
        Eq(NP1, np_closed),
        Eq(NT1, nt_closed).xreplace({NP: NP1}),
    ]
    cse_stages = cse(exprstages, numbered_symbols("n"))
    exprothers = (
        [
            Eq(LpF, LoverF),
            Eq(PpF, ppf),
            Eq(TpF, tpf),
            Eq(SWUpF, SWUoverF),
            Eq(SWUpP, SWUoverP),
        ]
        + [Eq(*z) for z in zip(xPi, xP)]
        + [Eq(*z) for z in zip(xTi, xT)]
    )
    exprothers = [e.xreplace({NP: NP1, NT: NT1}) for e in exprothers]
    cse_others = cse(exprothers, numbered_symbols("g"))
    exprops = count_ops(exprstages + exprothers)
    cse_ops = count_ops(cse_stages + cse_others)
    msg = "    reduced {0} ops to {1}".format(exprops, cse_ops)
    stat = _aggstatus(stat, msg, aggstat)

    # create function body
    ccode, repnames = cse_to_c(*cse_stages, indent=6, debug=debug)
    ccode_others, repnames_others = cse_to_c(*cse_others, indent=6, debug=debug)
    ccode += ccode_others
    repnames |= repnames_others

    msg = "  completed in {0:.3G} s".format(time.time() - start_time)
    stat = _aggstatus(stat, msg, aggstat)
    if aggstat:
        print(stat)
    return ccode, repnames, stat


_func_header1 = """pyne::enrichment::Cascade pyne::enrichment::solve_symbolic(pyne::enrichment::Cascade & orig_casc)
{
  pyne::enrichment::Cascade casc = orig_casc;
  int j = casc.j;
  int k = casc.k;
  double alpha = casc.alpha;
  double NP0 = casc.N;
  //double NT0 = casc.M;
  double Mstar = casc.Mstar;
  double xPj = casc.x_prod_j;
  //double xFj = casc.x_feed_j;
  double xTj = casc.x_tail_j;
  int ncomp = casc.mat_feed.comp.size();
  double LpF = -1.0, PpF = -1.0, TpF = -1.0, 
         SWUpF = -1.0, SWUpP = -1.0, 
         NP_b = -1.0, NP_sqrt_base = -1.0, NP_2a = -1.0, 
         NP1 = -1.0, NT1 = -1.0;
  double * MW = new double [ncomp];
  double * xP = new double [ncomp];
  double * xF = new double [ncomp];
  double * xT = new double [ncomp];
"""

_func_header2 = """ 
  int nuc;
  int i = 2;
  MW[0] = pyne::atomic_mass(j);
  MW[1] = pyne::atomic_mass(k);
  xF[0] = casc.mat_feed.comp[j];
  xF[1] = casc.mat_feed.comp[k];
  for(pyne::comp_iter ci = casc.mat_feed.comp.begin(); ci != casc.mat_feed.comp.end(); ci++)
  {
    nuc = (*ci).first;
    if (nuc == j || nuc == k)
        continue;
    MW[i] = pyne::atomic_mass(nuc);
    xF[i] = (*ci).second;
    i++;
  };

  switch (ncomp)
  {
"""

_func_footer = """ 
  };

  i = 2;
  casc.mat_prod.comp[j] = xP[0];
  casc.mat_prod.comp[k] = xP[1];
  casc.mat_tail.comp[j] = xT[0];
  casc.mat_tail.comp[k] = xT[1];
  for(pyne::comp_iter ci = casc.mat_feed.comp.begin(); ci != casc.mat_feed.comp.end(); ci++)
  {
    nuc = (*ci).first;
    if (nuc == j || nuc == k)
        continue;
    casc.mat_prod.comp[nuc] = xP[i];
    casc.mat_tail.comp[nuc] = xT[i];
    i++;
  };
  // must renormalize to eliminate numerical error
  casc.mat_prod.norm_comp();
  casc.mat_tail.norm_comp();
  casc.mat_prod.mass = PpF;
  casc.mat_tail.mass = TpF;

  casc.N = NP1;
  casc.M = NT1;
  casc.l_t_per_feed = LpF;
  casc.swu_per_feed = SWUpF;
  casc.swu_per_prod = SWUpP;

  delete [] MW;
  delete [] xP;
  delete [] xF;
  delete [] xT;

  return casc;
};
"""


def _mapable_cgen_ncomp(kwargs):
    return cgen_ncomp(**kwargs)


def cgen_func(max_ncomp=40, debug=False):
    """Generate C function to compute multicoponent enrichment cascades for
    a number of components between 3 and max_ncomp.
    """
    ncomps = range(3, max_ncomp + 1)
    if 1 == NPROCS:
        ncomp_kwargs = [{"ncomp": n, "debug": debug, "aggstat": False} for n in ncomps]
        cgened = map(_mapable_cgen_ncomp, ncomp_kwargs)
    elif 1 < NPROCS:
        ncomp_kwargs = [{"ncomp": n, "debug": debug, "aggstat": True} for n in ncomps]
        pool = multiprocessing.Pool(NPROCS)
        cgened = pool.map(_mapable_cgen_ncomp, ncomp_kwargs)
    else:
        raise ValueError("NPROCS must be greater than or equal to 1")
    cases = ""
    repnames = set()
    for ncomp, (ccode_ncomp, repnames_ncomp, statmsg) in zip(ncomps, cgened):
        cases += "    case {0}:\n".format(ncomp)
        cases += ccode_ncomp
        cases += "      break;\n"
        repnames |= repnames_ncomp
        logging.info(statmsg)
    repdeclare = "  double " + repnames.pop() + " = 0.0,\n"
    repdectemp = "         {0} = 0.0"
    repdeclare += ",\n".join([repdectemp.format(r) for r in repnames])
    repdeclare += ";\n"
    ccode = _func_header1
    ccode += repdeclare
    ccode += _func_header2
    ccode += cases
    ccode += _func_footer
    return ccode


_header_file_template = r""" 
/// \file {hfname}
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief A multicomponent enrichment cascade solver using
///     a symbolic solution to the mass flow rate equations.

/*********************************************************/
/***            Symbolic Enrichment Functions          ***/
/*** WARNING: This file is auto-generated.             ***/
/***                  DO NOT MODIFY!!!                 ***/
/*********************************************************/

#ifndef PYNE_OU4PO4TJDBDM5PY4VKAVL7JCSM
#define PYNE_OU4PO4TJDBDM5PY4VKAVL7JCSM

#include <math.h>

#ifndef PYNE_IS_AMALGAMATED
#include "enrichment_cascade.h"
#endif

namespace pyne {{
namespace enrichment {{

  /// A multicomponent enrichment cascade solver using     
  /// a symbolic solution to the mass flow rate equations.
  /// \param orig_casc The original state of the cascade.
  /// \return A cascade solved for new N, M, and total flow
  ///         rates.
  Cascade solve_symbolic(Cascade & orig_casc);

// end enrichment
}};
// end pyne
}};

#endif
"""

_source_file_header_template = """ 
/*********************************************************/
/***            Symbolic Enrichment Functions          ***/
/*** WARNING: This file is auto-generated.             ***/
/***                  DO NOT MODIFY!!!                 ***/
/*********************************************************/
#ifndef PYNE_IS_AMALGAMATED
#include "{hfname}"
#endif

"""


def cgen_header_file(hfname="temp"):
    """Generates a valid C/C++ header file for multicomponent enrichment cascades."""
    hcode = _header_file_template.format(hfname=os.path.split(hfname)[-1])
    return hcode


def cgen_source_file(hfname="temp", max_ncomp=40, debug=False):
    """Generates a valid C/C++ source file for multicomponent enrichment cascades."""
    ccode = _source_file_header_template.format(hfname=os.path.split(hfname)[-1])
    ccode += cgen_func(max_ncomp, debug=debug)
    return ccode


def cgen_file(
    filename="temp", header_filename=None, lang="C++", max_ncomp=40, debug=False
):
    """Generate C/C++ header and source file to compute multicoponent enrichment
    cascades for a number of components between 3 and max_ncomp. The filename
    argument should not end in extension ('.h', '.c', or '.cpp') as it will be
    appended automatically.
    """
    logfile = "sme{0}.log".format(max_ncomp)
    if os.path.exists(logfile):
        os.remove(logfile)
    logging.basicConfig(filename=logfile, level=logging.DEBUG)
    hfname = filename + ".h" if header_filename is None else header_filename
    sfname = filename + "." + {"C": "c", "C++": "cpp", "CPP": "cpp"}[lang.upper()]
    logging.info("header filename: " + hfname)
    logging.info("source filename: " + sfname)
    logging.info("language: " + lang)
    logging.info("maximum number of components: {0}".format(max_ncomp))
    logging.info("debug enabled: {0}".format(debug))
    hcode = cgen_header_file(hfname)
    ccode = cgen_source_file(hfname, max_ncomp, debug=debug)
    with open(hfname, "w") as f:
        f.write(hcode)
    with open(sfname, "w") as f:
        f.write(ccode)


if __name__ == "__main__":
    cgen_file(max_ncomp=3)
