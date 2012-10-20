"""Generate the API for symbolic multicomponent enrichment using sympy.
Note that there is a bug in sympy v0.7.2 that prevents the cse() function
from being used with infinities.  For a work around see [1].

1. https://groups.google.com/forum/#!msg/sympy/YL1R_hR6OKQ/axKrCsCSMQsJ
"""

from sympy import Symbol, pprint, latex, diff, count_ops, simplify, cse, Eq, Q, \
    log, logcombine, Abs, exp, sqrt, series, separate, powsimp, collect, expand
from sympy.solvers import solve

def cgen_ncomp(ncomp=3):
    """Generates a C function for ncomp (int) number of components.
    The jth key component is always in the first position and the kth
    key component is always in the second.
    """

r = range(0, 3)
j = 0
k = 1

alpha = Symbol('alpha', positive=True, real=True)
LpF = Symbol('LpF', positive=True, real=True)
NP = Symbol('NP', positive=True, real=True)  # Enrichment Stages
NT = Symbol('NT', positive=True, real=True)  # De-enrichment Stages
NP0 = Symbol('NP0', positive=True, real=True)  # Enrichment Stages Initial Guess
NT0 = Symbol('NT0', positive=True, real=True)  # De-enrichment Stages Initial Guess
Mstar = Symbol('Mstar', positive=True, real=True)
MW = [Symbol('MW[{0}]'.format(i), positive=True, real=True) for i in r]
beta = [alpha**(Mstar - MWi) for MWi in MW]

xF = [Symbol('xF[{0}]'.format(i), positive=True, real=True) for i in r]
xPj = Symbol('xPj', positive=True, real=True)
xFj = xF[j]
xTj = Symbol('xTj', positive=True, real=True)
ppf = (xFj - xTj)/(xPj - xTj)
tpf = (xFj - xPj)/(xTj - xFj)

xP = [((ppf*xF[i]*(beta[i]**(NT+1) - 1))/(beta[i]**(NT+1) - beta[i]**(-NP))) for i in r]
xT = [((tpf*xF[i]*(1 - beta[i]**(-NP)))/(beta[i]**(NT+1) - beta[i]**(-NP))) for i in r]

rfeed = xFj / xF[k]
rprod = xPj / xP[k]
rtail = xTj / xT[k]

numer = [ppf*xP[i]*log(rprod) + tpf*xT[i]*log(rtail) - xF[i]*log(rfeed) for i in r]
denom = [log(beta[j]) * ((beta[i] - 1.0)/(beta[i] + 1.0)) for i in r]
LoverF = reduce(sumtwo, [n/d for n, d in zip(numer, denom)])

prod_constraint = (xPj/xFj)*ppf - (beta[j]**(NT+1) - 1)/(beta[j]**(NT+1) - beta[j]**(-NP))
tail_constraint = (xTj/xFj)*(reduce(sumtwo, xT)) - (1 - beta[j]**(-NP))/(beta[j]**(NT+1) - beta[j]**(-NP))
xp_constraint = 1.0 - reduce(sumtwo, xP)
xf_constraint = 1.0 - reduce(sumtwo, xF)
xt_constraint = 1.0 - reduce(sumtwo, xT)

duc = {alpha:1.05, Mstar: 236.5, NP0:30.0, NT0:10.0, NP: 27.1835088212, NT: 13.3875092512, xPj: 0.05, xTj:0.0025}
duc.update({mwi:x for mwi, x in zip(MW, [235.043931368, 238.050789466, 234.040953616])})
duc.update({xf:x for xf, x in zip(xF, [0.0072, 0.992745, 0.000055])})

# <codecell>

d = dict(duc)
del d[Mstar]
evalLoverF = LoverF.xreplace(d)
x = np.linspace(235.5, 237.5, 21)
y = np.array([evalLoverF.evalf(subs={Mstar:xi}) for xi in x], dtype=float)

# <codecell>

fig = plt.figure(figsize=(7,7))
plt.plot(x, y, figure=fig)

# <codecell>

count_ops(LoverF)

# <codecell>

dLdMstar = diff(LoverF, Mstar)
count_ops(dLdMstar)

# <codecell>

d = dict(duc)
del d[Mstar]
evaldLdMstar = dLdMstar.xreplace(d)
x = np.linspace(235.5, 237.5, 21)
y = np.array([evaldLdMstar.evalf(subs={Mstar:xi}) for xi in x], dtype=float)

# <codecell>

fig = plt.figure(figsize=(7,7))
plt.plot(x, y, figure=fig)

# <codecell>

# This is NT(NP,...) and is correct!
nt_closed = solve(prod_constraint, NT)[0] 
nt_closed

# <codecell>

# This is nt_closed rewritten (by hand) to minimize the number of NP and M* instances
# in the expression.  Luckily this is only depends on the key component and remains 
# general no matter the number of components
nt_closed = (-MW[0]*log(alpha) + Mstar*log(alpha) + log(xTj) + log(xF[0] - xPj) - \
            log(alpha**((MW[0]-Mstar)*NP)*(xF[0]*xPj - xPj*xTj)/(xF[0]*xTj - xF[0]*xPj) + 1) - \
            log(xF[0]*xTj - xF[0]*xPj))/((MW[0] - Mstar)*log(alpha))
nt_closed

# <codecell>

d = dict(duc)
del d[NP]
nt_of_np = nt_closed.xreplace(d)
x = np.linspace(1.0, 50.0, 21)
y = np.array([nt_of_np.evalf(subs={NP:xi}) for xi in x], dtype=float)

# <codecell>

fig = plt.figure(figsize=(7,7))
plt.plot(x, y, figure=fig)

# <codecell>

loverf = LoverF.xreplace({NT: nt_closed})

# <codecell>

d = dict(duc)
del d[NP]
del d[Mstar]
loverfval = loverf.xreplace(d)
x = np.linspace(235.5, 237.5, 51), np.linspace(1.0, 50.0, 51)[::-1]
y = np.array([[loverfval.evalf(subs={NP:x2i, Mstar:x1i}) for x1i in x[0]] for x2i in x[1]], dtype=float)

# <codecell>

fig = plt.figure(figsize=(7,7))
plt.imshow(y, figure=fig, extent=(x[0][0], x[0][-1], x[1][-1], x[1][0]), aspect='auto', cmap='jet')
cb = plt.colorbar()

# <codecell>

# Define the constraint equation with which to solve NP
# This is chosen such to minimize the number of ops in the 
# derivatives (and thus np_closed).
#np_constraint = (xP[j]/reduce(sumtwo, xP) - xPj).xreplace({NT: nt_closed})
#np_constraint = (xP[j]- reduce(sumtwo, xP)*xPj).xreplace({NT: nt_closed})
#np_constraint = (xT[j]/reduce(sumtwo, xT) - xTj).xreplace({NT: nt_closed})
np_constraint = (xT[j]- reduce(sumtwo, xT)*xTj).xreplace({NT: nt_closed})

# <codecell>

d = dict(duc)
del d[NP]
del d[Mstar]
npcval = np_constraint.xreplace(d)
x = np.linspace(235.5, 237.5, 51), np.linspace(1.0, 50.0, 51)[::-1]
y = np.array([[npcval.evalf(subs={NP:x2i, Mstar:x1i}) for x1i in x[0]] for x2i in x[1]], dtype=float)

# <codecell>

fig = plt.figure(figsize=(7,7))
plt.imshow(y, figure=fig, extent=(x[0][0], x[0][-1], x[1][-1], x[1][0]), aspect='auto', cmap='jet',)
cb = plt.colorbar()
plt.contour(x[0], x[1], y, [-0.0001, 0.0, 0.0001], figure=fig, colors='k', linestyles='-')

# <codecell>

# Get second order, closed form approximation of NP
# symbolic derivatives
d1NP = diff(np_constraint, NP, 1).xreplace({NP: NP0})
d2NP = diff(np_constraint, NP, 2).xreplace({NP: NP0})/2.0
# taylor series polynomial coefficients, grouped by order
# f(x) = ax**2 + bx + c
a = d2NP
b = d1NP - 2*NP0*d2NP
c = np_constraint.xreplace({NP: NP0}) - NP0*d1NP + NP0*NP0*d2NP
# quadratic eq. (minus only)
np_closed2 = (-b - sqrt(b**2 - 4*a*c)) / (2*a) 

# <codecell>

print "np_constraint", count_ops(np_constraint)
print "d1NP", count_ops(d1NP)
print "d2NP", count_ops(d2NP)
print "np_closed2", count_ops(np_closed2)

# <codecell>

d = dict(duc)
del d[Mstar]
npc2val = np_closed2.xreplace(d)
x = np.linspace(235.5, 237.5, 21)
y = np.array([npc2val.evalf(subs={Mstar:xi}) for xi in x], dtype=complex)
y.real[y.imag != 0.0] = 0.0
y = y.real

# <codecell>

fig = plt.figure(figsize=(7,7))
plt.plot(x, y, figure=fig)

# <codecell>

# Get first order, closed form approximation of NP
np_closed1 = NP0 - d1NP / np_constraint.xreplace({NP: NP0})
print "np_closed1", count_ops(np_closed1)

# <codecell>

d = dict(duc)
del d[Mstar]
npc1val = np_closed1.xreplace(d)
x = np.linspace(235.5, 237.5, 51)
y = np.array([npc1val.evalf(subs={Mstar:xi}) for xi in x], dtype=complex)
y.real[y.imag != 0.0] = 0.0
y = y.real

# <codecell>

fig = plt.figure(figsize=(7,7))
plt.plot(x, y, figure=fig)

# <codecell>

# Generate common sub-expresions to write out
#cse_np_closed = cse(np_closed1)
cse_np_closed = cse(np_closed2)
cse_loverf = cse(loverf)

# <codecell>

print "cse_np_closed", count_ops(cse_np_closed)
print "cse_loverf", count_ops(cse_loverf)

# <codecell>

cse_system = cse([nt_closed, Eq(NP, np_closed2), loverf])

# <codecell>

print "cse_system", count_ops(cse_system)

# <codecell>



