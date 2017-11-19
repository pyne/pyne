$PROJECT = 'pyne'
$ACTIVITIES = [#'version_bump', 'changelog',
               'nose',
               #'sphinx', 'tag', 'conda_forge', 'ghpages', 'ghrelease'
               ]

$VERSION_BUMP_PATTERNS = [
    ('setup.py', "VERSION\s*=.*", "VERSION = '$VERSION'"),
    ('setup_sub.py', "VERSION\s*=.*", "VERSION = '$VERSION'"),
    ('pyne/__init__.py', "__version__\s*=.*", "__version__ = '$VERSION'"),
    ('src/utils.cpp', "std::string pyne::VERSION\s*=.*",
     'std::string pyne::VERSION = "$VERSION"'),
]

$CHANGELOG_FILENAME = 'CHANGELOG.rst'
$CHANGELOG_LATEST = 'docs/previous/$VERSION-release-notes.rst'
$CHANGELOG_TEMPLATE = 'TEMPLATE.rst'

$DOCKER_CONDA_DEPS = ['cmake', 'pkg-config', 'setuptools', 'gcc',
    'libgcc', 'libgfortran', 'python', 'blas', 'openblas', 'boost-cpp',
    'hdf5', 'bzip2', 'xz', 'moab', 'cython', 'numpy', 'pytables',
    'jinja2', 'nose', 'sphinx', 'numpydoc', 'cloud_sptheme',
    'sphinxcontrib-bibtex',
    ]
$DOCKER_INSTALL_COMMAND = './setup.py install && nuc_data_make'
$DOCKER_GIT_NAME = 'pyne'
$DOCKER_GIT_EMAIL = 'pyne@pyne.io'

$NOSE_ARGS = ['-w', 'tests']