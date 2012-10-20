"""Generate the API for symbolic multicomponent enrichment using sympy.
Note that there is a bug in sympy v0.7.2 that prevents the cse() function
from being used with infinities.  For a work around see [1].

1. https://groups.google.com/forum/#!msg/sympy/YL1R_hR6OKQ/axKrCsCSMQsJ
"""

from sympy import Symbol, pprint, latex, diff, count_ops, simplify, cse, Eq, Q, \
    log, logcombine, Abs, exp, sqrt, series, separate, powsimp, collect, expand
from sympy.solvers import solve
from sympy.utilities.iterables import numbered_symbols

def cgen_ncomp(ncomp=3, nporder=2):
    """Generates a C function for ncomp (int) number of components.
    The jth key component is always in the first position and the kth
    key component is always in the second.  The number of enrichment 
    stages (NP) is calculated via a taylor series approximation.  The
    order of this approximation may be set with nporder.  Only values
    of 1 or 2 are allowed.
    """
    r = range(0, ncomp)
    j = 0
    k = 1

    # setup-symbols
    alpha = Symbol('alpha', positive=True, real=True)
    LpF = Symbol('LpF', positive=True, real=True)
    NP = Symbol('NP', positive=True, real=True)    # Enrichment Stages
    NT = Symbol('NT', positive=True, real=True)    # De-enrichment Stages
    NP0 = Symbol('NP0', positive=True, real=True)  # Enrichment Stages Initial Guess
    NT0 = Symbol('NT0', positive=True, real=True)  # De-enrichment Stages Initial Guess
    Mstar = Symbol('Mstar', positive=True, real=True)
    MW = [Symbol('MW[{0}]'.format(i), positive=True, real=True) for i in r]
    beta = [alpha**(Mstar - MWi) for MWi in MW]

    xF = [Symbol('xF[{0}]'.format(i), positive=True, real=True) for i in r]
    xPi = [Symbol('xP[{0}]'.format(i), positive=True, real=True) for i in r]
    xTi = [Symbol('xT[{0}]'.format(i), positive=True, real=True) for i in r]
    xPj = Symbol('xPj', positive=True, real=True)
    xFj = xF[j]
    xTj = Symbol('xTj', positive=True, real=True)
    ppf = (xFj - xTj)/(xPj - xTj)
    tpf = (xFj - xPj)/(xTj - xFj)

    xP = [((ppf*xF[i]*(beta[i]**(NT+1) - 1))/(beta[i]**(NT+1) - beta[i]**(-NP))) \
                                                                            for i in r]
    xT = [((tpf*xF[i]*(1 - beta[i]**(-NP)))/(beta[i]**(NT+1) - beta[i]**(-NP))) \
                                                                            for i in r]
    rfeed = xFj / xF[k]
    rprod = xPj / xP[k]
    rtail = xTj / xT[k]

    # setup constraint equations
    numer = [ppf*xP[i]*log(rprod) + tpf*xT[i]*log(rtail) - xF[i]*log(rfeed) for i in r]
    denom = [log(beta[j]) * ((beta[i] - 1.0)/(beta[i] + 1.0)) for i in r]
    LoverF = sum([n/d for n, d in zip(numer, denom)])

    prod_constraint = (xPj/xFj)*ppf - (beta[j]**(NT+1) - 1)/\
                      (beta[j]**(NT+1) - beta[j]**(-NP))
    tail_constraint = (xTj/xFj)*(reduce(sumtwo, xT)) - (1 - beta[j]**(-NP))/\
                      (beta[j]**(NT+1) - beta[j]**(-NP))
    #xp_constraint = 1.0 - sum(xP)
    #xf_constraint = 1.0 - sum(xF)
    #xt_constraint = 1.0 - sum(xT)

    # This is NT(NP,...) and is correct!
    #nt_closed = solve(prod_constraint, NT)[0] 

    # However, this is NT(NP,...) rewritten (by hand) to minimize the number of NP 
    # and M* instances in the expression.  Luckily this is only depends on the key 
    # component and remains general no matter the number of components.
    nt_closed = (-MW[j]*log(alpha) + Mstar*log(alpha) + log(xTj) + log(xF[j] - xPj) - \
        log(alpha**((MW[j]-Mstar)*NP)*(xF[j]*xPj - xPj*xTj)/(xF[j]*xTj - xF[j]*xPj) + \
        1) - log(xF[j]*xTj - xF[j]*xPj))/((MW[j] - Mstar)*log(alpha))

    # new expression for normalized flow rate
    loverf = LoverF.xreplace({NT: nt_closed})

    # Define the constraint equation with which to solve NP. This is chosen such to 
    # minimize the number of ops in the derivatives (and thus np_closed).  Other, 
    # more verbose possibilities are commented out.
    #np_constraint = (xP[j]/reduce(sumtwo, xP) - xPj).xreplace({NT: nt_closed})
    #np_constraint = (xP[j]- reduce(sumtwo, xP)*xPj).xreplace({NT: nt_closed})
    #np_constraint = (xT[j]/reduce(sumtwo, xT) - xTj).xreplace({NT: nt_closed})
    np_constraint = (xT[j]- sum(xT)*xTj).xreplace({NT: nt_closed})

    # get closed form approximation of NP via symbolic derivatives
    d0NP = np_constraint.xreplace({NP: NP0})
    d1NP = diff(np_constraint, NP, 1).xreplace({NP: NP0})
    if 1 == nporder:
        np_closed = NP0 - d1NP / d0NP
    elif 2 == nporder:
        d2NP = diff(np_constraint, NP, 2).xreplace({NP: NP0})/2.0
        # taylor series polynomial coefficients, grouped by order
        # f(x) = ax**2 + bx + c
        a = d2NP
        b = d1NP - 2*NP0*d2NP
        c = d0NP - NP0*d1NP + NP0*NP0*d2NP
        # quadratic eq. (minus only)
        np_closed = (-b - sqrt(b**2 - 4*a*c)) / (2*a) 
    else:
        raise ValueError("nporder must be 1 or 2")

    # generate cse for writing out
    cse_stages = cse([Eq(NT, nt_closed), Eq(NP, np_closed)], numbered_symbols('n'))
    cse_others = cse([Eq(LpF, LoverF)] + [Eq(*z) for z in zip(xPi, xP)] + \
                     [Eq(*z) for z in zip(xTi, xT)], numbered_symbols('g'))

