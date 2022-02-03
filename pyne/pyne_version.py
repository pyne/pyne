# Set PYNE version
PYNE_MAJOR_VERSION = 0
PYNE_MINOR_VERSION = 7
PYNE_PATCH_VERSION = 6

PYNE_VERSION = f'{PYNE_MAJOR_VERSION}.{PYNE_MINOR_VERSION}.{PYNE_PATCH_VERSION}'

def write_cpp_header():

    with open('pyne_version.h','w') as pvh:
        header_text = f'''#ifndef PYNE_VERSION_HEADER
#define PYNE_VERSION_HEADER

#define PYNE_VERSION {PYNE_VERSION}
#define PYNE_VERSION_STRING "{PYNE_VERSION}"

#endif // PYNE_VERSION_HEADER
'''
        pvh.write(header_text)


if __name__ == "__main__":
    write_cpp_header()

