from cython.operator cimport dereference

cimport pyne.cpp_endf

def read_endf(fname):
    pyne.cpp_endf.read_endf(fname)
    
def get_library_contents():
    print pyne.cpp_endf.get_library_contents()

def load_mt(mat, mf, mt):
    pyne.cpp_endf.load_dataset_to_API(mat, mf, mt)
    
def get_mt_list():
    return pyne.cpp_endf.library.mt_451.mt_list