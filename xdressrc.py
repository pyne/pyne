from xdress.utils import apiname
from xdress.doxygen import default_doxygen_config

package = 'pyne'
packagedir = 'pyne'

plugins = ('xdress.stlwrap', 'xdress.autoall', 'xdress.autodescribe',
           'xdress.descfilter', 'xdress.doxygen', 'xdress.cythongen')

doxygen_config = {'PROJECT_NAME': 'PYNE',
                  'EXTRACT_ALL': False,  # Note usage of python False
                  'GENERATE_DOCBOOK': False,
                  'GENERATE_LATEX': True  # Could be 'YES' or True
                  }

doxyfile_name = './build/tally_doxyfile'

includes = ['/opt/local/include']

extra_types = 'extra_types'  # non-default value

stlcontainers = [
    ('set', 'str'),
    ('set', 'int32'),
    ('pair', 'double', 'double'),
    ('map', 'str', 'str'),
    ('map', 'str', 'int32'),
    ('map', 'int32', 'str'),
    ('map', 'str', 'uint32'),
    ('map', 'uint32', 'str'),
    ('map', 'str', 'float64'),
    ('map', 'uint32', 'uint32'),
    ('map', 'int32', 'int32'),
    ('map', 'int32', 'float64'),
    ('map', 'int32', 'complex'),
    ('map', 'uint32', 'float64'),
    ('map', 'str', ('vector', 'float64')),
    ('map', 'int32', ('vector', 'float64')),
    ('vector', ('vector', 'int')),
    ('vector', ('vector', 'double')),
    ('vector', ('vector', ('vector', 'double'))),
    ('vector', ('pair', 'double','double')),
    ]

# stlcontainers_module = 'stlcontainers'

classes = [apiname('Tally',  'src/tally.*', incfiles='tally.h'),
           apiname('mt_base', 'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('mt451',  'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('mt452_mf1', 'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('mt455_mf1', 'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('mt456_mf1', 'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('mt458_mf1', 'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('mt460_mf1', 'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('mtfpy_mf8', 'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('mt457_mf8', 'src/endf_mt.*', tarbase='endf2',
                   incfiles=['endf_mt.h']),
           apiname('endf_id', 'src/endf.[hc]*', tarbase='endf2',
                   incfiles=['endf.h']),
           apiname('library', 'src/endf.[hc]*', tarbase='endf2',
                   incfiles=['endf.h'])]

functions = []

skipmethods = {'library': ['gen_content_list', ]}
skipattrs = {'library': ['contents', ]}
