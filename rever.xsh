$PROJECT = 'pyne'
$ACTIVITIES = ['version_bump', 'changelog']

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

