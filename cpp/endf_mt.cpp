#ifndef PYNE_IS_AMALGAMATED
  #include "endf_base.h"
  #include "endf_mt.h"
#endif

pyne::endf::mt_451 pyne::endf::read_451(std::ifstream &infile) {
    mt_451 mt_451;

    control cs = read_cont(infile);
    mt_451.nuc_id = cs.c1;
    mt_451.awr = cs.c2;
    mt_451.mat = cs.mat;
    mt_451.mf = cs.mf;
    mt_451.mt = cs.mt;

    mt_451.lrp = cs.l1;
    mt_451.lfi = cs.l2;
    mt_451.nlib = cs.n1;
    mt_451.nmod = cs.n2;


    cs = read_cont(infile);
    mt_451.elis = cs.c1;
    mt_451.sta = cs.c2;
    mt_451.lis = cs.l1;
    mt_451.liso = cs.l2;
    //mt_451.elis = cs.n1;
    mt_451.nfor = cs.n2;


    cs = read_cont(infile);
    mt_451.awi = cs.c1;
    mt_451.emax = cs.c2;
    mt_451.lrel = cs.l1;
    //mt_451.liso = cs.l2;
    mt_451.nsub = cs.n1;
    mt_451.nver = cs.n2;


    cs = read_cont(infile);
    mt_451.temp = cs.c1;
    //mt_451.emax = cs3.c2;
    mt_451.ldrv = cs.l1;
    //mt_451.liso = cs.l2;
    mt_451.nwd = cs.n1;
    mt_451.nxc = cs.n2;
    std::string line;
    for (int i = 0; i < mt_451.nwd; ++i) {
        getline(infile, line);
    }

    for (int i = 0; i < mt_451.nxc; ++i) {
        cs = read_cont(infile);
        int tmp_arr[] = {cs.l1, cs.l2, cs.n1, cs.n2};
        std::vector<int> tmp (tmp_arr, tmp_arr + sizeof(tmp_arr) / sizeof(int) );
        mt_451.mt_list.push_back(tmp);
    }
    getline(infile, line); // get send line
    return mt_451;
}

pyne::endf::mt_fpy_8 pyne::endf::read_fpy_8(std::ifstream &infile) {
    mt_fpy_8 mt_fpy;

    control cs = read_cont(infile);
    mt_fpy.nuc_id = cs.c1;
    mt_fpy.awr = cs.c2;
    mt_fpy.mat = cs.mat;
    mt_fpy.mf = cs.mf;
    mt_fpy.mt = cs.mt;

    mt_fpy.le = cs.l1;
    std::vector<std::vector<double> > tmp;
    for (int i = 0; i < mt_fpy.le; ++i) {
        list lst = read_list(infile);
        mt_fpy.e.push_back(lst.c1);
        mt_fpy.i.push_back(lst.l1);

        tmp.clear();
        for (int j = 0; j < lst.data.size()/4; ++j) {
          double dt[] = {lst.data[j*4], lst.data[j*4+1],
                lst.data[j*4+2], lst.data[j*4+3]};
          tmp.push_back(std::vector<double> (dt, dt + sizeof(dt)/sizeof(double)));
        }
        mt_fpy.yields.push_back(tmp);
    }
    std::string line;
    getline(infile, line);
    return mt_fpy;
}

pyne::endf::mt_452_1 pyne::endf::read_452_1(std::ifstream &infile) {
  mt_452_1 mt_452;

  control cs = read_cont(infile);
  mt_452.nuc_id = cs.c1;
  mt_452.awr = cs.c2;
  mt_452.mat = cs.mat;
  mt_452.mf = cs.mf;
  mt_452.mt = cs.mt;

  mt_452.lnu = cs.l2;

  if (mt_452.lnu == 1) {
    list lst = read_list(infile);
    mt_452.poly = lst.data;
  }else if (mt_452.lnu == 2) {
    tab1 tab = read_tab1(infile);
    mt_452.nbt = tab.nbt;
    mt_452.intn = tab.intn;
    mt_452.eint = tab.x;
    mt_452.nu_e = tab.y;
  }
  std::string line;
  getline(infile, line);
  return mt_452;
}

pyne::endf::mt_455_1 pyne::endf::read_455_1(std::ifstream &infile) {
  mt_455_1 mt_455;

  control cs = read_cont(infile);
  mt_455.nuc_id = cs.c1;
  mt_455.awr = cs.c2;
  mt_455.mat = cs.mat;
  mt_455.mf = cs.mf;
  mt_455.mt = cs.mt;

  mt_455.ldg = cs.l1;
  mt_455.lnu = cs.l2;

  if ((mt_455.lnu == 2) && (mt_455.ldg == 0)) {
    list lst = read_list(infile);
    tab1 tab = read_tab1(infile);
    mt_455.lambdas = lst.data;
    mt_455.nbt = tab.nbt;
    mt_455.intn = tab.intn;
    mt_455.eint = tab.x;
    mt_455.nu_d = tab.y;

  }
  if ((mt_455.lnu == 2) && (mt_455.ldg == 1)) {
    tab2 tab = read_tab2(infile);
    int nv;
    for (int i = 0; i < tab.l2; ++i) {
      list lst = read_list(infile);
      nv = lst.npl/2;
      std::vector<double> lambdas(nv, 0.0);
      std::vector<double> alphas(nv, 0.0);
      for (int j = 0; j < nv; ++j) {
        lambdas[j] = lst.data[j*2];
        alphas[j] = lst.data[j*2 + 1];
      }
      mt_455.lambda_arr.push_back(lambdas);
      mt_455.alpha_arr.push_back(alphas);
    }
    tab1 tab1 = read_tab1(infile);
    mt_455.nbt = tab1.nbt;
    mt_455.intn = tab1.intn;
    mt_455.eint = tab1.x;
    mt_455.nu_d = tab1.y;
  }
  if ((mt_455.lnu == 1) && (mt_455.ldg == 0)) {
    list lst = read_list(infile);
    mt_455.lambda_arr.push_back(lst.data);
    tab1 tab = read_tab1(infile);
    mt_455.nbt = tab.nbt;
    mt_455.intn = tab.intn;
    mt_455.nu_d = tab.y;
  }
  if ((mt_455.lnu == 1) && (mt_455.ldg == 1)) {
    tab2 tab = read_tab2(infile);
    mt_455.ne = tab.nbt;
    mt_455.einti = tab.intn;
    int nv;
    for (int i = 0; i < tab.l2; ++i) {
      list lst = read_list(infile);
      nv = lst.npl/2;
      std::vector<double> lambdas(nv, 0.0);
      std::vector<double> alphas(nv, 0.0);
      for (int j = 0; j < nv; ++j) {
        lambdas[j] = lst.data[j*2];
        alphas[j] = lst.data[j*2 + 1];
      }
      mt_455.lambda_arr.push_back(lambdas);
      mt_455.alpha_arr.push_back(alphas);
    }
    tab1 tab1 = read_tab1(infile);
    mt_455.nbt = tab1.nbt;
    mt_455.intn = tab1.intn;
    mt_455.eint = tab1.x;
    mt_455.nu_d = tab1.y;
  }
  std::string line;
  getline(infile, line);
  return mt_455;
}

pyne::endf::mt_456_1 pyne::endf::read_456_1(std::ifstream &infile) {
  mt_456_1 mt_456;
  control cs = read_cont(infile);
  mt_456.nuc_id = cs.c1;
  mt_456.awr = cs.c2;
  mt_456.mat = cs.mat;
  mt_456.mf = cs.mf;
  mt_456.mt = cs.mt;

  mt_456.lnu = cs.l2;

  if (mt_456.lnu == 1) {
    list lst = read_list(infile);
    mt_456.nu = lst.data;
  }else if (mt_456.lnu == 2) {
    tab1 tab = read_tab1(infile);
    mt_456.nbt = tab.nbt;
    mt_456.intn = tab.intn;
    mt_456.eint = tab.x;
    mt_456.nu_e = tab.y;
  }

  std::string line;
  getline(infile, line);
  return mt_456;
}

pyne::endf::mt_458_1 pyne::endf::read_458_1(std::ifstream &infile) {
  mt_458_1 mt_458;
  control cs = read_cont(infile);
  mt_458.nuc_id = cs.c1;
  mt_458.awr = cs.c2;
  mt_458.mat = cs.mat;
  mt_458.mf = cs.mf;
  mt_458.mt = cs.mt;

  list lst = read_list(infile);

  mt_458.efr = std::vector<double> (lst.l2 + 1, 0);
  mt_458.defr = std::vector<double> (lst.l2 + 1, 0);
  mt_458.enp = std::vector<double> (lst.l2 + 1, 0);
  mt_458.denp = std::vector<double> (lst.l2 + 1, 0);
  mt_458.end = std::vector<double> (lst.l2 + 1, 0);
  mt_458.dend = std::vector<double> (lst.l2 + 1, 0);
  mt_458.egp = std::vector<double> (lst.l2 + 1, 0);
  mt_458.degp = std::vector<double> (lst.l2 + 1, 0);
  mt_458.egd = std::vector<double> (lst.l2 + 1, 0);
  mt_458.degd = std::vector<double> (lst.l2 + 1, 0);
  mt_458.eb = std::vector<double> (lst.l2 + 1, 0);
  mt_458.deb = std::vector<double> (lst.l2 + 1, 0);
  mt_458.enu = std::vector<double> (lst.l2 + 1, 0);
  mt_458.denu = std::vector<double> (lst.l2 + 1, 0);
  mt_458.er = std::vector<double> (lst.l2 + 1, 0);
  mt_458.der = std::vector<double> (lst.l2 + 1, 0);
  mt_458.et = std::vector<double> (lst.l2 + 1, 0);
  mt_458.det = std::vector<double> (lst.l2 + 1, 0);

  for (int i = 0; i < lst.l2 + 1; ++i) {
    mt_458.efr[i] = lst.data[i*18];
    mt_458.defr[i] = lst.data[i*18+1];
    mt_458.enp[i] = lst.data[i*18+2];
    mt_458.denp[i] = lst.data[i*18+3];
    mt_458.end[i] = lst.data[i*18+4];
    mt_458.dend[i] = lst.data[i*18+5];
    mt_458.egp[i] = lst.data[i*18+6];
    mt_458.degp[i] = lst.data[i*18+7];
    mt_458.egd[i] = lst.data[i*18+8];
    mt_458.degd[i] = lst.data[i*18+9];
    mt_458.eb[i] = lst.data[i*18+10];
    mt_458.deb[i] = lst.data[i*18+11];
    mt_458.enu[i] = lst.data[i*18+12];
    mt_458.denu[i] = lst.data[i*18+13];
    mt_458.er[i] = lst.data[i*18+14];
    mt_458.der[i] = lst.data[i*18+15];
    mt_458.et[i] = lst.data[i*18+16];
    mt_458.det[i] = lst.data[i*18+17];
  }

  std::string line;
  getline(infile, line);
  return mt_458;
}

pyne::endf::mt_460_1 pyne::endf::read_460_1(std::ifstream &infile) {
  mt_460_1 mt_460;
  control cs = read_cont(infile);

  mt_460.nuc_id = cs.c1;
  mt_460.awr = cs.c2;
  mt_460.mat = cs.mat;
  mt_460.mf = cs.mf;
  mt_460.mt = cs.mt;

  mt_460.lo = cs.l1;
  mt_460.ng = cs.n1;
  mt_460.elist = std::vector<double>(mt_460.ng, 0.0);

  if (mt_460.lo == 1) {
    for (int i = 0; i < mt_460.ng; ++i) {
      tab1 tab = read_tab1(infile);
      mt_460.elist[i] = tab.c1;
      mt_460.nbt.push_back(tab.nbt);
      mt_460.intn.push_back(tab.intn);
      mt_460.tint.push_back(tab.x);
      mt_460.t.push_back(tab.y);
    }
  }else if(mt_460.lo == 2) {
    list lst = read_list(infile);
    mt_460.lambdas = lst.data;
  }

  std::string line;
  getline(infile, line);
  return mt_460;
}
