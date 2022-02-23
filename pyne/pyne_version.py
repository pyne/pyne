import os

# Set PYNE version
PYNE_MAJOR_VERSION = 0
PYNE_MINOR_VERSION = 7
PYNE_PATCH_VERSION = 6

PYNE_VERSION = '{0}.{1}.{2}'.format(PYNE_MAJOR_VERSION, PYNE_MINOR_VERSION, PYNE_PATCH_VERSION)

def write_cpp_header(path="."):

    with open(os.path.join(path , 'pyne_version.h') , 'w') as pvh:
        header_text = '''#ifndef PYNE_VERSION_HEADER
#define PYNE_VERSION_HEADER

#define PYNE_VERSION {0}
#define PYNE_VERSION_STRING "{0}"

#endif // PYNE_VERSION_HEADER
'''.format(PYNE_VERSION)
        pvh.write(header_text)


if __name__ == "__main__":
    write_cpp_header()

