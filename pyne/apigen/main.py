import argparse

def main():
    parser = argparse.ArgumentParser("Generates PyNE API")
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--enrich', action='store_true', default=False)
    parser.add_argument('--stlwrap', action='store_true', default=False)
    ns = parser.parse_args()

    if ns.enrich:
        print "generating symbolic enrichment cascades"
        import enrich_multi_sym
        enrich_multi_sym.cgen_file(filename="cpp/enrichment_symbolic15", 
                                   header_filename="cpp/enrichment_symbolic.h",
                                   max_ncomp=15, 
                                   debug=ns.debug)
    if ns.stlwrap:
        print "generating C++ standard library wrappers & converters"
        import stlwrap
        template = [('py2c_set', 'int'),
                    ('py2c_set', 'str'),
                    ('py2c_map', 'int', 'int'),
                    ('py2c_map', 'int', 'double'),
                    ('py2c_map', 'str', 'int'),
                    ('py2c_map', 'int', 'str'),
                    ('py2c_map', 'str', 'double'),
                    ('set', 'int'),
                    ('set', 'str'),
                    ('map', 'str', 'str'),
                    ('map', 'str', 'int'),
                    ('map', 'int', 'str'),
                    ('map', 'str', 'uint'),
                    ('map', 'uint', 'str'),
                    ('map', 'uint', 'uint'),
                    ('map', 'str', 'double'),
                    ('map', 'int', 'int'),
                    ('map', 'int', 'double'),
                    ('map', 'uint', 'double'),
                    ('map', 'int', 'complex'),
                    ('map', 'int', 'vector[double]'),
                    ('map', 'str', 'vector[double]'),
                    ]
        stlwrap.genfiles(template, 
                         fname='pyne/stlcontainers', 
                         testname="pyne/tests/test_stlcontainers")

if __name__ == '__main__':
    main()
