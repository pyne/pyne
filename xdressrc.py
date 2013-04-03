package = 'pyne'
packagedir = 'pyne'
sourcedir = 'cpp'

extra_types = 'extra_types'  # non-default value

stlcontainers = [
    ('set', 'str'),
    ('set', 'int32'),
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
    ]

#stlcontainers_module = 'stlcontainers'

classes = []

functions = []
