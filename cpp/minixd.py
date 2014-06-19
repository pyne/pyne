def get_after(before, line):
    parts = line.split()
    getnext = False
    after = None
    for item in parts:
        if getnext == True:
            after = item
            getnext = False
        if item == before:
            getnext = True
    return after


supported_types = {'int' : 'int', 'double' : 'double',
                   'std::vector<int>' : 'vector[int]',
                   'std::vector<double>' : 'vector[double]',
                   'std::vector<std::vector<int> >' : 'vector[vector[int]]',
                   'std::vector<std::vector<double> >' : 'vector[vector[double]]',
                   'std::vector<std::vector<std::vector<double> > >' : 'vector[vector[vector[double]]]',
                   }

def process_header(fname):
    with open(fname, 'r') as f:
        namespace = None
        fillstruct = False
        structs = []
        objects = dict()
        openbr = 0
        brnumber = 0
        for line in f.readlines():
            openbr += line.count("{")
            openbr -= line.count("}")
            if "namespace" in line:
                namespace = get_after('namespace', line)
            if "typedef struct" in line:
                structname = get_after('struct', line)
                if ':' in line:
                    parent = get_after(':', line)
                else:
                    parent = None
                struct = dict()
                struct['name'] = structname
                struct['parent'] = parent
                struct['items'] = dict()
                struct['docs'] = dict()
                struct['unmatched'] = dict()
                brnumber = openbr
                fillstruct = True
            if brnumber <= openbr and fillstruct == True and ";" in line:
                parts = line.split()
                name = None
                docs = ''
                if len(line.split('///<')) == 2:
                    docs = line.split('///<')[1]
                dtype = ''
                for part in parts:
                    if ';' in part:
                        name = part.rstrip(';')
                        break
                    dtype += part + " "
                dtype = dtype.rstrip()
                if name is not None:
                    if dtype in supported_types:
                        struct['items'][name] = supported_types[dtype]
                        struct['docs'][name] = docs.rstrip('\n')
                    else:
                        struct['unmatched'][name] = dtype
                        struct['docs'][name] = docs.rstrip('\n')
            if fillstruct == True and brnumber > openbr:
                fillstruct = False
                structs.append(struct)
            if ";" in line:
                parts = line.split()
                name = None
                docs = ''
                if len(line.split('///<')) == 2:
                    docs = line.split('///<')[1]
                dtype = ''
                for part in parts:
                    if ';' in part:
                        name = part.rstrip(';')
                        break
                    dtype += part + " "
                dtype = dtype.rstrip()
                if name is not None:
                    if dtype in supported_types:
                        objects['types'] = supported_types[dtype]
                        objects['docs'] = docs.rstrip('\n')
                    else:
                        struct['unmatched'][name] = dtype
                        struct['docs'][name] = docs.rstrip('\n')
        return namespace, structs


def build_pxd(header, altname= None):
    namespace, structs = process_header(header)
    if altname is not None:
        pxdname = altname
    else:
        pxdname = header.rstrip('.h')
    with open('cpp_' + pxdname + ".pxd", 'w') as f:
        f.write("from libcpp.vector cimport vector\n")
        f.write("\n")
        f.write('cdef extern from "' + header + '" namespace "' + namespace
                + '":\n')
        for struct in structs:
            if struct['parent'] is not None:
                f.write('    cdef cppclass ' + struct['name'] + '(' + struct['parent'] + '):\n')
            else:
                f.write('    cdef cppclass ' + struct['name'] + ':\n')
            for item in struct['items']:
                f.write(8 * ' ' + struct['items'][item] + ' ' + item + " #" + struct['docs'][item] + '\n')
            for item in struct['unmatched']:
                f.write(8 * ' ' + '# ' + struct['unmatched'][item] + " #" + struct['docs'][item] + '\n')
            f.write('\n')

def build_pyx(header, pxd, pyxname):
    namespace, structs = process_header(header)
    with open(pyxname + ".pyx", 'w') as f:
        f.write("cimport " + pxd + "\n")
        f.write("\n")
        for struct in structs:
            f.write('cdef class _' + struct['name'] + ':\n')
            f.write('    cdef ' + pxd + '.' + struct['name'] + ' *thisptr\n')
            # don't allocate or de-allocate we'll do this by hand later
            f.write('    def __cinit__(self):\n')
            f.write('        pass\n')
            f.write('    def __dealloc__(self):\n')
            f.write('        pass\n')
            getstr = 8 * ' ' + 'def __get__(self): return self.thisptr.'
            if struct['parent'] is not None:
                for pstruct in structs:
                    if struct['parent'] == pstruct['name']:
                        for item in pstruct['items']:
                            f.write('    property ' + item + ':\n')
                            f.write(8 * ' ' + '"""' + pstruct['docs'][item] + '"""' + '\n')
                            f.write(getstr + item + '\n')
                        for item in pstruct['unmatched']:
                            f.write('    #property ' + item + ':\n')
            for item in struct['items']:
                f.write('    property ' + item + ':\n')
                f.write(8 * ' ' + '"""' + struct['docs'][item] + '"""' + '\n')
                f.write(getstr + item + '\n')
            for item in struct['unmatched']:
                f.write('    #property ' + item + ':\n')
            f.write("\n")
