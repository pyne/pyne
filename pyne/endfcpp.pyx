from cython.operator cimport dereference

cimport pyne.cpp_endf

def read_endf(fname):
    pyne.cpp_endf.read_endf(fname)
    
def print_mt_list():
    print pyne.cpp_endf.library.mt_451.mt_list