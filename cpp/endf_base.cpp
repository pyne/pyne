#ifndef PYNE_IS_AMALGAMATED
  #include "endf_base.h"
  #include "utils.h"
#endif

pyne::endf::control_struct pyne::endf::read_cont(std::ifstream &infile) {
    control_struct cts;
    std::string line;
    getline(infile, line);
    cts.c1 = pyne::endftod(line.substr(0, 11).c_str());
    cts.c2 = pyne::endftod(line.substr(11, 11).c_str());
    cts.l1 = atoi(line.substr(22, 11).c_str());
    cts.l2 = atoi(line.substr(33, 11).c_str());
    cts.n1 = atoi(line.substr(44, 11).c_str());
    cts.n2 = atoi(line.substr(55, 11).c_str());
    cts.mat = atoi(line.substr(66, 4).c_str());
    cts.mf = atoi(line.substr(70, 2).c_str());
    cts.mt = atoi(line.substr(72, 3).c_str());
    return cts;
}

pyne::endf::list_struct pyne::endf::read_list(std::ifstream &infile) {
    list_struct lst;
    std::string line;
    getline(infile, line);
    lst.c1 = pyne::endftod(line.substr(0, 11).c_str());
    lst.c2 = pyne::endftod(line.substr(11, 11).c_str());
    lst.l1 = atoi(line.substr(22, 11).c_str());
    lst.l2 = atoi(line.substr(33, 11).c_str());
    lst.npl = atoi(line.substr(44, 11).c_str());
    lst.n2 = atoi(line.substr(55, 11).c_str());
    lst.mat = atoi(line.substr(66, 4).c_str());
    lst.mf = atoi(line.substr(70, 2).c_str());
    lst.mt = atoi(line.substr(72, 3).c_str());
    lst.data = std::vector<double>(lst.npl, 0.0);
    int npl = lst.npl;
    int n = 0;
    while (npl > 0) {
        getline(infile, line);
        for (int i = 0; i < 6; ++i){
            lst.data[n]=pyne::endftod((line.substr(i*11, 11)).c_str());
            ++n;
            --npl;
            if (npl == 0)
                break;
        }
    }
    return lst;
}

pyne::endf::tab1_struct pyne::endf::read_tab1(std::ifstream &infile){
  tab1_struct tab1;
  std::string line;
  getline(infile, line);
  tab1.c1 = pyne::endftod(line.substr(0, 11).c_str());
  tab1.c2 = pyne::endftod(line.substr(11, 11).c_str());
  tab1.l1 = atoi(line.substr(22, 11).c_str());
  tab1.l2 = atoi(line.substr(33, 11).c_str());
  tab1.nr = atoi(line.substr(44, 11).c_str());
  tab1.np = atoi(line.substr(55, 11).c_str());
  tab1.mat = atoi(line.substr(66, 4).c_str());
  tab1.mf = atoi(line.substr(70, 2).c_str());
  tab1.mt = atoi(line.substr(72, 3).c_str());
  tab1.nbt = std::vector<int> (tab1.nr, 0); ///<
  tab1.intn = std::vector<int> (tab1.nr, 0);
  tab1.x = std::vector<double> (tab1.np, 0.0);
  tab1.y = std::vector<double> (tab1.np, 0.0);
  int nr = tab1.nr;
  int n = 0;
  while (nr > 0) {
      getline(infile, line);
      for (int i = 0; i < 3; ++i){
          tab1.nbt[n]=atoi(line.substr(i*22, 11).c_str());
          tab1.intn[n]=atoi(line.substr(i*22+11, 11).c_str());
          ++n;
          --nr;
          if (nr == 0)
              break;
      }
  }
  int np = tab1.np;
  n = 0;
  while (np > 0) {
    getline(infile, line);
    for (int i = 0; i < 3; ++i){
        tab1.x[n]=pyne::endftod((line.substr(i*22, 11)).c_str());
        tab1.y[n]=pyne::endftod((line.substr(i*22+11, 11)).c_str());
        ++n;
        --np;
        if (np == 0)
            break;
    }
  }
  return tab1;
}

pyne::endf::tab2_struct pyne::endf::read_tab2(std::ifstream &infile){
  tab2_struct tab2;
  std::string line;
  getline(infile, line);
  tab2.c1 = pyne::endftod(line.substr(0, 11).c_str());
  tab2.c2 = pyne::endftod(line.substr(11, 11).c_str());
  tab2.l1 = atoi(line.substr(22, 11).c_str());
  tab2.l2 = atoi(line.substr(33, 11).c_str());
  tab2.nr = atoi(line.substr(44, 11).c_str());
  tab2.nz = atoi(line.substr(55, 11).c_str());
  tab2.mat = atoi(line.substr(66, 4).c_str());
  tab2.mf = atoi(line.substr(70, 2).c_str());
  tab2.mt = atoi(line.substr(72, 3).c_str());
  tab2.nbt = std::vector<int> (tab2.nr, 0); ///<
  tab2.intn = std::vector<int> (tab2.nr, 0);
  int nr = tab2.nr;
  int n = 0;
  while (nr > 0) {
      getline(infile, line);
      for (int i = 0; i < 3; ++i){
          tab2.nbt[n]=atoi(line.substr(i*22, 11).c_str());
          tab2.intn[n]=atoi(line.substr(i*22+11, 11).c_str());
          ++n;
          --nr;
          if (nr == 0)
              break;
      }
  }
  return tab2;
}
