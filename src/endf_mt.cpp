#ifndef PYNE_IS_AMALGAMATED
  #include "endf_base.h"
  #include "endf_mt.h"
#endif

pyne::mt451 pyne::read_mt451(std::ifstream &infile) {
    mt451 mt451;

    control cs = read_cont(infile);
    mt451.nuc_id = cs.c1;
    mt451.awr = cs.c2;
    mt451.mat = cs.mat;
    mt451.mf = cs.mf;
    mt451.mt = cs.mt;

    mt451.lrp = cs.l1;
    mt451.lfi = cs.l2;
    mt451.nlib = cs.n1;
    mt451.nmod = cs.n2;


    cs = read_cont(infile);
    mt451.elis = cs.c1;
    mt451.sta = cs.c2;
    mt451.lis = cs.l1;
    mt451.liso = cs.l2;
    //mt451.elis = cs.n1;
    mt451.nfor = cs.n2;


    cs = read_cont(infile);
    mt451.awi = cs.c1;
    mt451.emax = cs.c2;
    mt451.lrel = cs.l1;
    //mt451.liso = cs.l2;
    mt451.nsub = cs.n1;
    mt451.nver = cs.n2;


    cs = read_cont(infile);
    mt451.temp = cs.c1;
    //mt451.emax = cs3.c2;
    mt451.ldrv = cs.l1;
    //mt451.liso = cs.l2;
    mt451.nwd = cs.n1;
    mt451.nxc = cs.n2;
    std::string line;
    for (int i = 0; i < mt451.nwd; ++i) {
        getline(infile, line);
    }

    for (int i = 0; i < mt451.nxc; ++i) {
        cs = read_cont(infile);
        int tmp_arr[] = {cs.l1, cs.l2, cs.n1, cs.n2};
        std::vector<int> tmp (tmp_arr, tmp_arr + sizeof(tmp_arr) / sizeof(int) );
        mt451.mt_list.push_back(tmp);
    }
    getline(infile, line); // get send line
    return mt451;
}

pyne::mtfpy_mf8 pyne::read_mtfpy_mf8(std::ifstream &infile) {
    mtfpy_mf8 mt_fpy;

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

pyne::mt452_mf1 pyne::read_mt452_mf1(std::ifstream &infile) {
  mt452_mf1 mt452;

  control cs = read_cont(infile);
  mt452.nuc_id = cs.c1;
  mt452.awr = cs.c2;
  mt452.mat = cs.mat;
  mt452.mf = cs.mf;
  mt452.mt = cs.mt;

  mt452.lnu = cs.l2;

  if (mt452.lnu == 1) {
    list lst = read_list(infile);
    mt452.poly = lst.data;
  }else if (mt452.lnu == 2) {
    tab1 tab = read_tab1(infile);
    mt452.nbt = tab.nbt;
    mt452.intn = tab.intn;
    mt452.eint = tab.x;
    mt452.nu_e = tab.y;
  }
  std::string line;
  getline(infile, line);
  return mt452;
}

pyne::mt455_mf1 pyne::read_mt455_mf1(std::ifstream &infile) {
  mt455_mf1 mt455;

  control cs = read_cont(infile);
  mt455.nuc_id = cs.c1;
  mt455.awr = cs.c2;
  mt455.mat = cs.mat;
  mt455.mf = cs.mf;
  mt455.mt = cs.mt;

  mt455.ldg = cs.l1;
  mt455.lnu = cs.l2;

  if ((mt455.lnu == 2) && (mt455.ldg == 0)) {
    list lst = read_list(infile);
    tab1 tab = read_tab1(infile);
    mt455.lambdas = lst.data;
    mt455.nbt = tab.nbt;
    mt455.intn = tab.intn;
    mt455.eint = tab.x;
    mt455.nu_d = tab.y;

  }
  if ((mt455.lnu == 2) && (mt455.ldg == 1)) {
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
      mt455.lambda_arr.push_back(lambdas);
      mt455.alpha_arr.push_back(alphas);
    }
    tab1 tab1 = read_tab1(infile);
    mt455.nbt = tab1.nbt;
    mt455.intn = tab1.intn;
    mt455.eint = tab1.x;
    mt455.nu_d = tab1.y;
  }
  if ((mt455.lnu == 1) && (mt455.ldg == 0)) {
    list lst = read_list(infile);
    mt455.lambda_arr.push_back(lst.data);
    tab1 tab = read_tab1(infile);
    mt455.nbt = tab.nbt;
    mt455.intn = tab.intn;
    mt455.nu_d = tab.y;
  }
  if ((mt455.lnu == 1) && (mt455.ldg == 1)) {
    tab2 tab = read_tab2(infile);
    mt455.ne = tab.nbt;
    mt455.einti = tab.intn;
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
      mt455.lambda_arr.push_back(lambdas);
      mt455.alpha_arr.push_back(alphas);
    }
    tab1 tab1 = read_tab1(infile);
    mt455.nbt = tab1.nbt;
    mt455.intn = tab1.intn;
    mt455.eint = tab1.x;
    mt455.nu_d = tab1.y;
  }
  std::string line;
  getline(infile, line);
  return mt455;
}

pyne::mt456_mf1 pyne::read_mt456_mf1(std::ifstream &infile) {
  mt456_mf1 mt456;
  control cs = read_cont(infile);
  mt456.nuc_id = cs.c1;
  mt456.awr = cs.c2;
  mt456.mat = cs.mat;
  mt456.mf = cs.mf;
  mt456.mt = cs.mt;

  mt456.lnu = cs.l2;

  if (mt456.lnu == 1) {
    list lst = read_list(infile);
    mt456.nu = lst.data;
  }else if (mt456.lnu == 2) {
    tab1 tab = read_tab1(infile);
    mt456.nbt = tab.nbt;
    mt456.intn = tab.intn;
    mt456.eint = tab.x;
    mt456.nu_e = tab.y;
  }

  std::string line;
  getline(infile, line);
  return mt456;
}

pyne::mt458_mf1 pyne::read_mt458_mf1(std::ifstream &infile) {
  mt458_mf1 mt458;
  control cs = read_cont(infile);
  mt458.nuc_id = cs.c1;
  mt458.awr = cs.c2;
  mt458.mat = cs.mat;
  mt458.mf = cs.mf;
  mt458.mt = cs.mt;

  list lst = read_list(infile);

  mt458.efr = std::vector<double> (lst.l2 + 1, 0);
  mt458.defr = std::vector<double> (lst.l2 + 1, 0);
  mt458.enp = std::vector<double> (lst.l2 + 1, 0);
  mt458.denp = std::vector<double> (lst.l2 + 1, 0);
  mt458.end = std::vector<double> (lst.l2 + 1, 0);
  mt458.dend = std::vector<double> (lst.l2 + 1, 0);
  mt458.egp = std::vector<double> (lst.l2 + 1, 0);
  mt458.degp = std::vector<double> (lst.l2 + 1, 0);
  mt458.egd = std::vector<double> (lst.l2 + 1, 0);
  mt458.degd = std::vector<double> (lst.l2 + 1, 0);
  mt458.eb = std::vector<double> (lst.l2 + 1, 0);
  mt458.deb = std::vector<double> (lst.l2 + 1, 0);
  mt458.enu = std::vector<double> (lst.l2 + 1, 0);
  mt458.denu = std::vector<double> (lst.l2 + 1, 0);
  mt458.er = std::vector<double> (lst.l2 + 1, 0);
  mt458.der = std::vector<double> (lst.l2 + 1, 0);
  mt458.et = std::vector<double> (lst.l2 + 1, 0);
  mt458.det = std::vector<double> (lst.l2 + 1, 0);

  for (int i = 0; i < lst.l2 + 1; ++i) {
    mt458.efr[i] = lst.data[i*18];
    mt458.defr[i] = lst.data[i*18+1];
    mt458.enp[i] = lst.data[i*18+2];
    mt458.denp[i] = lst.data[i*18+3];
    mt458.end[i] = lst.data[i*18+4];
    mt458.dend[i] = lst.data[i*18+5];
    mt458.egp[i] = lst.data[i*18+6];
    mt458.degp[i] = lst.data[i*18+7];
    mt458.egd[i] = lst.data[i*18+8];
    mt458.degd[i] = lst.data[i*18+9];
    mt458.eb[i] = lst.data[i*18+10];
    mt458.deb[i] = lst.data[i*18+11];
    mt458.enu[i] = lst.data[i*18+12];
    mt458.denu[i] = lst.data[i*18+13];
    mt458.er[i] = lst.data[i*18+14];
    mt458.der[i] = lst.data[i*18+15];
    mt458.et[i] = lst.data[i*18+16];
    mt458.det[i] = lst.data[i*18+17];
  }

  std::string line;
  getline(infile, line);
  return mt458;
}

pyne::mt460_mf1 pyne::read_mt460_mf1(std::ifstream &infile) {
  mt460_mf1 mt460;
  control cs = read_cont(infile);

  mt460.nuc_id = cs.c1;
  mt460.awr = cs.c2;
  mt460.mat = cs.mat;
  mt460.mf = cs.mf;
  mt460.mt = cs.mt;

  mt460.lo = cs.l1;
  mt460.ng = cs.n1;
  mt460.elist = std::vector<double>(mt460.ng, 0.0);

  if (mt460.lo == 1) {
    for (int i = 0; i < mt460.ng; ++i) {
      tab1 tab = read_tab1(infile);
      mt460.elist[i] = tab.c1;
      mt460.nbt.push_back(tab.nbt);
      mt460.intn.push_back(tab.intn);
      mt460.tint.push_back(tab.x);
      mt460.t.push_back(tab.y);
    }
  }else if(mt460.lo == 2) {
    list lst = read_list(infile);
    mt460.lambdas = lst.data;
  }

  std::string line;
  getline(infile, line);
  return mt460;
}

pyne::mt457_mf8 pyne::read_mt457_mf8(std::ifstream &infile) {
  mt457_mf8 mt457;
  control cs = read_cont(infile);

  mt457.nuc_id = cs.c1;
  mt457.awr = cs.c2;
  mt457.mat = cs.mat;
  mt457.mf = cs.mf;
  mt457.mt = cs.mt;

  mt457.lis = cs.l1;
  mt457.liso = cs.l2;
  mt457.nst = cs.n1;
  mt457.nsp = cs.n2;

  if (mt457.nst == 1) {
    //boring stable nucleus
  } else {
    list lst = read_list(infile); //decay heat data
    mt457.erel = std::vector<std::pair<double,double> >(lst.npl/2);
    for (int i = 0; i < lst.npl/2; ++i) {
      mt457.erel[i] = std::make_pair(lst.data[i * 2], lst.data[i * 2 + 1]);
    }
    lst = read_list(infile); //Q values and branching ratios
    mt457.styp = std::vector<double>(mt457.nsp);
    mt457.lcon = std::vector<int>(mt457.nsp);
    mt457.ner = std::vector<int>(mt457.nsp);
    mt457.fd = std::vector<std::pair<double,double> >(mt457.nsp);
    mt457.eav = std::vector<std::pair<double,double> >(mt457.nsp);
    mt457.fc = std::vector<std::pair<double,double> >(mt457.nsp);
    for (int i = 0; i < mt457.nsp; ++i) {
      list spec = read_list(infile); //Basic spectrum information
      mt457.styp[i] = spec.c2;
      mt457.lcon[i] = spec.l1;
      mt457.ner[i] = spec.l2;
      mt457.fd[i] = std::make_pair(spec.data[0], spec.data[1]);
      mt457.eav[i] = std::make_pair(spec.data[2], spec.data[3]);
      mt457.fc[i] = std::make_pair(spec.data[4], spec.data[5]);

      if (spec.l1 != 1) { //LCON != 1
        for (int j = 0; i < spec.n2; ++j) {
          list rad = read_list(infile); //info for each discrete radiation
          mt457.er.push_back(std::make_pair(rad.c1, rad.c2));
          mt457.rtyp.push_back(rad.data[0]);
          mt457.type.push_back(rad.data[1]);
          mt457.ri.push_back(std::make_pair(rad.data[2], rad.data[3]));
          mt457.ris.push_back(std::make_pair(rad.data[4], rad.data[5]));
          if (rad.npl == 12) {
            mt457.ricc.push_back(std::make_pair(rad.data[6], rad.data[7]));
            mt457.rick.push_back(std::make_pair(rad.data[8], rad.data[9]));
            mt457.ricl.push_back(std::make_pair(rad.data[10], rad.data[11]));
          }
        }
      }
      if (spec.l1 != 0) { //LCON != 0
        tab1 tab = read_tab1(infile);//continuum radiation
        if (tab.l2 != 0) { //LCON != 0 and LCOV != 0
          list cov = read_list(infile); //covariance information
        }
      }
    }
  }
  std::string line;
  getline(infile, line);
  return mt457;
}
