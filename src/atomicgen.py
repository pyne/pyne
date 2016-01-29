#!/usr/bin/python
import re

# print the header
def print_header_file(filename):
    header_file = ""
    header_file += "/// \/file atomic_nuclear_data.h\n"
    header_file += "/// \/author Andrew Davis (andrew.davis@wisc.edu)\n"
    header_file += "///\n"
    header_file += "/// \/brief Implements all the fundamental atomic & nuclear data data\n"
    header_file += "#include <map>\n"
    header_file += "\n"
    header_file += "namespace pyne\n"
    header_file += "{\n"
    header_file += "  /// main function to be called when you wish to load the nuclide data \n"
    header_file += "  /// into memory \n"
    header_file += "  void _load_atomic_mass_map_memory();\n"
    header_file += "  /// function to create mapping from nuclides in id form\n" 
    header_file += "  /// to their atomic masses\n"
    header_file += "  \n"
    header_file += "  void _insert_atomic_mass_map();\n"
    header_file += "  \n"
    header_file += "  /// function to create mapping from nuclides in id form \n"
    header_file += "  /// to their natural abundances\n"
    header_file += "  void _insert_abund_map();\n"
    header_file += "  \n"
    header_file += "  /// Mapping from nuclides in id form to their natural abundances\n"
    header_file += "  extern std::map<int,double> natural_abund_map;\n"
    header_file += "  \n"
    header_file += "  /// Mapping from nuclides in id form to their atomic masses.\n"
    header_file += "  extern std::map<int,double> atomic_mass_map;\n"
    header_file += "  \n"
    header_file += "  /// Mapping from nuclides in id form to the associated error in \n"
    header_file += "  /// abdundance \n"
    header_file += "  extern std::map<int,double> atomic_mass_error_map;\n"
    header_file += "} // namespace pyne\n"
    
    f = open(filename,'w')
    f.write(header_file)
    f.close()

# print the masses map
def print_atomic_mass_errors():
    atomic_mass_error = ""

    file = open('../pyne/dbgen/mass.mas12','r')
    amdc_regex = re.compile('[ \d-]*? (\d{1,3})[ ]{1,4}(\d{1,3}) [A-Z][a-z]? .*? (\\d{1,3}) ([ #.\d]{5,12}) ([ #.\d]+)[ ]*?$')
    for line in file:
        m = amdc_regex.search(line)
        if m is None:
            continue
            
        nuc = (10000000 * int(m.group(1))) + (10000 * int(m.group(2)))
        error = 1E-6 * float(m.group(5).strip().replace('#', ''))
        atomic_mass_error = "  atomic_mass_error["+str(nuc)+"] = "+str(error)+";\n"
    file.close()
    return atomic_mass_error

# print the masses map
def print_atomic_mass():
    atomic_mass = ""
    file = open('../pyne/dbgen/mass.mas12','r')
    amdc_regex = re.compile('[ \d-]*? (\d{1,3})[ ]{1,4}(\d{1,3}) [A-Z][a-z]? .*? (\\d{1,3}) ([ #.\d]{5,12}) ([ #.\d]+)[ ]*?$')
    for line in file:
        m = amdc_regex.search(line)
        if m is None:
            continue
            
        nuc = (10000000 * int(m.group(1))) + (10000 * int(m.group(2)))
        mass = float(m.group(3)) + 1E-6 * float(m.group(4).strip().replace('#',''))
        atomic_mass += "  atomic_mass_map["+str(nuc)+"] = "+str(mass)+";\n"
    file.close()
    return atomic_mass

# print the abundances map
def print_abundances():
    isotopic_abundances = ""
    file = open('../pyne/dbgen/abundances.txt','r')

    for line in file:
        # tokenise the line
        splitted = line.split(' ')
        tokens = []
        for item in splitted:
            if item != '':
                tokens.append(item)
                # ignore the comment line
        if tokens[0] != "#":
            tokens[3] = tokens[3].replace("\n","")
            isotopic_abundances += "  natural_abund_map["+str((int(tokens[0])*10000000)+(int(tokens[2])*10000))+"] = "+str(tokens[3])+";\n"
    return isotopic_abundances

def print_cpp_file(filename):
    #empty string
    cpp_file = ""
    cpp_file += "// Implements basic nuclear data functions.\n"
    cpp_file += "#ifndef PYNE_IS_AMALGAMATED\n"
    cpp_file += "#include \"atomic_data.h\"\n"
    cpp_file += "#endif\n"
    cpp_file += "  \n"
    cpp_file += "void pyne::_load_atomic_mass_map_memory() { \n"
    cpp_file += "  // header version of atomic weight table data \n"
    cpp_file += "  //see if the data table is already loaded\n"
    cpp_file += "  if(!atomic_mass_map.empty()) {\n"
    cpp_file += "    return;\n"
    cpp_file += "  } else { \n"
    cpp_file += "    _insert_atomic_mass_map();\n"
    cpp_file += "  }\n"
    cpp_file += "  //see if the data table is already loaded\n"
    cpp_file += "  if(!natural_abund_map.empty()) {\n"
    cpp_file += "    return;\n"
    cpp_file += "  } else { \n"
    cpp_file += "    _insert_abund_map();\n"
    cpp_file += "  }\n"
    cpp_file += "}\n"
    cpp_file += "\n"
    cpp_file += "void pyne::_insert_atomic_mass_map() { \n"
    cpp_file += print_atomic_mass()
    cpp_file += "}\n"
    cpp_file += "\n"
    cpp_file += "void pyne::_insert_abund_map() { \n"
    cpp_file += print_abundances()
    cpp_file += "}\n"

    f = open(filename,'w')
    f.write(cpp_file)
    f.close()

if __name__ == '__main__':
    print_cpp_file("atomic_data.cpp")
    print_header_file("atomic_data.h")
