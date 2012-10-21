import argparse

import enrich_multi_sym

def main():
    debug = False
    enrich_multi_sym.cgen_file(filename="cpp/enrichment_symbolic", max_ncomp=40, 
                               debug=debug)

if __name__ == '__main__':
    main()
