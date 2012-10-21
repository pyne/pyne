"""Utility functions for pyne.apigen"""

import re

from sympy.utilities.codegen import codegen

def cse_to_c(replacements, reduced_exprs, indent=2, debug=False):
    """Converts the return value sympy.cse() to a single C code snippet.
    """
    ccode = ""
    ws = " " * indent

    rtn_pattern = re.compile('\s*?return (.*)')
    equ_pattern = re.compile('\s*?([\w\[\]]+)\s*?=\s*?(\S.*)')
    repline_template = '{ind}double {name} = {expr}\n'
    redline_template = '{ind}{name} = {expr}\n'
    redrtnline_template = '{ind}return {expr}\n'

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
        if debug:
            ccode += ws + 'std::cout << "{0} = " << {0} << "\\n";\n'.format(repname)

    for redexpr in reduced_exprs:
        gencode = codegen(("redname", redexpr), "C", "temp", header=False, empty=False)
        genexpr = gencode[0][1]
        if redexpr.is_Equality:
            m = equ_pattern.search(genexpr)
            genname = m.group(1)
            genexpr = m.group(2)
            genline = redline_template.format(ind=ws, name=genname, expr=genexpr)
            if debug:
                ccode += ws+'std::cout << "{0} = " << {0} << "\\n";\n'.format(genname)
        else:
            genexpr = rtn_pattern.search(genexpr).group(1)
            genline = redrtnline_template.format(ind=ws, expr=genexpr)
        ccode += genline
    return ccode
