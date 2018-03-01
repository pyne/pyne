$PROJECT = $GITHUB_ORG = $GITHUB_REPO = 'pyne'
$ACTIVITIES = ['version_bump', 'changelog', 'nose', 'sphinx',
               'tag', 'push_tag', 'conda_forge', 'ghpages', 'ghrelease'
               ]

$VERSION_BUMP_PATTERNS = [
    ('setup.py', "VERSION\s*=.*", "VERSION = '$VERSION'"),
    ('setup_sub.py', "VERSION\s*=.*", "VERSION = '$VERSION'"),
    ('pyne/__init__.py', "__version__\s*=.*", "__version__ = '$VERSION'"),
    ('src/utils.cpp', "std::string pyne::VERSION\s*=.*",
     'std::string pyne::VERSION = "$VERSION";'),
]

$CHANGELOG_FILENAME = 'CHANGELOG.rst'
$CHANGELOG_LATEST = 'docs/previous/$VERSION-release-notes.rst'
$CHANGELOG_TEMPLATE = 'TEMPLATE.rst'

$DOCKER_APT_DEPS = ['libc6', 'libc6-i386', 'libc6-dev', 'libc-dev', 'gcc']
$DOCKER_CONDA_DEPS = ['cmake', 'pkg-config', 'setuptools', 'gcc',
    'libgcc', 'libgfortran', 'make', 'scipy',
    'python', 'blas', 'openblas', 'boost-cpp',
    'hdf5', 'bzip2', 'xz', 'moab=4', 'cython', 'numpy', 'pytables',
    'jinja2', 'nose', 'sphinx', 'numpydoc', 'cloud_sptheme',
    'sphinxcontrib-bibtex', 'prettytable','nbconvert'
    ]
$DOCKER_INSTALL_COMMAND = ('rm -rf build && ./setup.py install --clean && '
                           'cd $HOME && nuc_data_make')
$DOCKER_GIT_NAME = 'pyne'
$DOCKER_GIT_EMAIL = 'pyne@pyne.io'

$NOSE_ARGS = ['-w', 'tests']

$PUSH_TAG_REMOTE = 'git@github.com:pyne/pyne.git'
$PUSH_TAG_TARGET = 'master'
$GHPAGES_REPO = 'git@github.com:pyne/pyne.github.com.git'