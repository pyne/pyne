#ifndef PYNE_IS_AMALGAMATED
#include "endf.h"
#endif

//Core functions for parsing endf lines and datasets
pyne::endf::endf_library_struct pyne::endf::library;

pyne::endf::endf_library_struct::~endf_library_struct() {
    //delete structs that have been allocated
}

pyne::endf::control_struct pyne::endf::read_cont(std::ifstream &infile) {
    control_struct cts;
    std::string line;
    getline(infile, line);
    char * cstr = new char [11];
    std::strcpy (cstr, (line.substr(0,11)).c_str());
    cts.c1 = pyne::endftod(cstr);
    std::strcpy (cstr, (line.substr(11,11)).c_str());
    cts.c2 = pyne::endftod(cstr);
    cts.l1 = stoi(line.substr(22,11));
    cts.l2 = stoi(line.substr(33,11));
    cts.n1 = stoi(line.substr(44,11));
    cts.n2 = stoi(line.substr(55,11));
    cts.mat = stoi(line.substr(66, 4));
    cts.mf = stoi(line.substr(70, 2));
    cts.mt = stoi(line.substr(72, 3));
    return cts;
}

pyne::endf::list_struct pyne::endf::read_list(std::ifstream &infile) {
    list_struct lst;
    std::string line;
    getline(infile, line);
    char * cstr = new char [11];
    std::strcpy (cstr, (line.substr(0,11)).c_str());
    lst.c1 = pyne::endftod(cstr);
    std::strcpy (cstr, (line.substr(11,11)).c_str());
    lst.c2 = pyne::endftod(cstr);
    lst.l1 = stoi(line.substr(22,11));
    lst.l2 = stoi(line.substr(33,11));
    lst.npl = stoi(line.substr(44,11));
    lst.n2 = stoi(line.substr(55,11));
    lst.data = std::vector<double>(lst.npl,0.0);
    int npl = lst.npl;
    int n = 0;
    while (npl > 0) {      
        getline(infile, line);
        for (int i = 0; i < 6; ++i){
            std::strcpy (cstr, (line.substr(i*11,11)).c_str());
            lst.data[n]=pyne::endftod(cstr);
            ++n;
            --npl;
            if (npl == 0) 
                break;
        }
    }
    return lst;
}


void pyne::endf::read_451(std::ifstream &infile) {
    //library.mt_451 = new mt_451_struct();
    library.mt_451.mt_list.clear();

    control_struct cs = read_cont(infile);
    library.mt_451.nuc_id = cs.c1;
    library.mt_451.awr = cs.c2;
    library.mt_451.mat = cs.mat;
    library.mt_451.mf = cs.mf;
    library.mt_451.mt = cs.mt;
    
    library.mt_451.lrp = cs.l1;
    library.mt_451.lfi = cs.l2;
    library.mt_451.nlib = cs.n1;
    library.mt_451.nmod = cs.n2;
    
    
    cs = read_cont(infile);
    library.mt_451.elis = cs.c1;
    library.mt_451.sta = cs.c2;
    library.mt_451.lis = cs.l1;
    library.mt_451.liso = cs.l2;
    //library.mt_451.elis = cs.n1;
    library.mt_451.nfor = cs.n2;
    
    
    cs = read_cont(infile);
    library.mt_451.awi = cs.c1;
    library.mt_451.emax = cs.c2;
    library.mt_451.lrel = cs.l1;
    //library.mt_451.liso = cs.l2;
    library.mt_451.nsub = cs.n1;
    library.mt_451.nver = cs.n2;
    
    
    cs = read_cont(infile);
    library.mt_451.temp = cs.c1;
    //library.mt_451.emax = cs3.c2;
    library.mt_451.ldrv = cs.l1;
    //library.mt_451.liso = cs.l2;
    library.mt_451.nwd = cs.n1;
    library.mt_451.nxc = cs.n2;
    std::string line;
    for (int i = 0; i < library.mt_451.nwd; ++i) {
        getline(infile, line);
    }
    
    for (int i = 0; i < library.mt_451.nxc; ++i) {
        cs = read_cont(infile);
        int tmp_arr[] = {cs.l1,cs.l2,cs.n1,cs.n2};
        std::vector<int> tmp (tmp_arr, tmp_arr + sizeof(tmp_arr) / sizeof(int) );
        library.mt_451.mt_list.push_back(tmp);
    }
    getline(infile, line); // get send line
}

void pyne::endf::read_454(std::ifstream &infile) {
    //library.mt_454 = new mt_fpy_struct();
    library.mt_454.yields.clear();
    
    control_struct cs = read_cont(infile);
    library.mt_454.nuc_id = cs.c1;
    library.mt_454.awr = cs.c2;
    library.mt_454.mat = cs.mat;
    library.mt_454.mf = cs.mf;
    library.mt_454.mt = cs.mt;
 
    library.mt_454.le = cs.l1;
    for (int i = 0; i < library.mt_454.le; ++i) {
        list_struct lst = read_list(infile);
        library.mt_454.e.push_back(lst.c1);
        library.mt_454.i.push_back(lst.l1);
        library.mt_454.yields.push_back(lst.data);
    }
}

void pyne::endf::read_459(std::ifstream &infile) {
    //library.mt_459 = new mt_fpy_struct();

    library.mt_459.yields.clear();
    control_struct cs = read_cont(infile);
    library.mt_459.nuc_id = cs.c1;
    library.mt_459.awr = cs.c2;
    library.mt_459.mat = cs.mat;
    library.mt_459.mf = cs.mf;
    library.mt_459.mt = cs.mt;
    
    library.mt_459.le = cs.l1;
    for (int i = 0; i < library.mt_459.le; ++i) {
        list_struct lst = read_list(infile);
        library.mt_459.e.push_back(lst.c1);
        library.mt_459.i.push_back(lst.l1);
        library.mt_459.yields.push_back(lst.data);
    }
}

void pyne::endf::read_endf(std::string filenm) {
    std::ifstream infile;
    std::string line;
    infile.open(filenm.c_str());
    getline(infile,line);//revision and other stuff
    
    int place;
    int mt;
    int nmat = 1;// 451 should be first if it's not don't run 
    std::string temp;
    while(nmat >0) {
         place = infile.tellg();
         getline(infile,line);
         infile.seekg(place);
         mt = stoi(temp.assign(line, 72, 3));
         
         //This should check to see what this is the end of but for now
         //we'll assume it's a modern single mat per file situation.
         if (mt == 0)
            continue;
         else
            --nmat;
         if (mt == 451) {
            //~library();
            read_451(infile);
            nmat = library.mt_451.nxc - 1;
            //TODO: determine section positions so we can check on parsers
         }else if (mt == 454) {
            read_454(infile);
         }else if (mt == 459) {
            read_459(infile);
         }
         
         
    }
}