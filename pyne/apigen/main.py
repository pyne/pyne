import argparse

import enrich_multi_sym

def main():
    debug = False
    enrich_multi_sym.cgen_file(filename="cpp/enrichment_symbolic05", 
                               header_filename="cpp/enrichment_symbolic.h",
                               max_ncomp=5, 
                               debug=debug)

if __name__ == '__main__':
    main()
