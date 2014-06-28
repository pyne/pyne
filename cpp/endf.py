from xdress.utils import Arg

mod = {
    'docstring': "Python wrapper for endf library.",
    'library': {
        'methods': {
            (('get', 'mt_451'), ('comp', 'endf_id')): {
                'return': 'mt_451',
                'defaults': ((Arg.NONE, None),)},
            (('get', 'mt_451'), ('mat','int'), ('mf', 'int'), ('mt', 'int')): {
                'return': 'mt_451',
                'defaults': ((Arg.NONE, None),(Arg.NONE, None),(Arg.NONE, None))},
            }
        },
    }
