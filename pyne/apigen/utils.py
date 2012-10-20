"""Utility functions for pyne.apigen"""

import re

from sympy.utilities.codegen import codegen

def cse_to_c(replacements, reduced_exprs, indent=2):
    """Converts the return value sympy.cse() to a single C code snippet.
    """
    ccode = ""
    ws = " " * indent

    rtn_pattern = re.compile('\s*?return (.*)')
    equ_pattern = re.compile('\s*?(\w+)\s*?=\s*?(\w.*)')
    repline_template = '{ind}double {name} = {expr}\n'

    for repsym, repexpr in replacements:
        repname = str(repsym)
        gencode = codegen((repname, repexpr), "C", repname, header=False, empty=False)
        genexpr = gencode[0][1]
        if repexpr.is_Equality:
            genexpr = equ_pattern.search(genexpr).group(2)
        else:
            genexpr = rtn_pattern.search(genexpr).group(1)
        genline = repline_template.format(ind=ws, name=repname, expr=genexpr)
        ccode += genline

    for redexpr in reduced_exprs:
        gencode = codegen((repname, repexpr), "C", repname, header=False, empty=False)
        genexpr = gencode[0][1]
        if repexpr.is_Equality:
            genexpr = equ_pattern.search(genexpr).group(0)
        else:
            genexpr = rtn_pattern.search(genexpr).group(0)
        ccode += genexpr + "\n"
    return ccode
