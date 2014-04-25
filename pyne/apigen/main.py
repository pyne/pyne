from __future__ import print_function
import argparse
from warnings import warn
from pyne.utils import VnVWarning

warn(__name__ + " is not yet V&V compliant.", VnVWarning)

def main():
    parser = argparse.ArgumentParser("Generates PyNE API")
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--enrich', action='store_true', default=False)
    ns = parser.parse_args()

    if ns.enrich:
        print("generating symbolic enrichment cascades")
        import enrich_multi_sym
        enrich_multi_sym.cgen_file(filename="cpp/enrichment_symbolic15", 
                                   header_filename="cpp/enrichment_symbolic.h",
                                   max_ncomp=15, 
                                   debug=ns.debug)

if __name__ == '__main__':
    main()
