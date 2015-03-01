#ifndef PYNE_IS_AMALGAMATED
  #include "endf_base.h"
  #include "endf_mt.h"
#endif

pyne::endf::mt451 pyne::endf::read_mt451(std::ifstream &infile) {
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

pyne::endf::mtfpy_mf8 pyne::endf::read_mtfpy_mf8(std::ifstream &infile) {
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

pyne::endf::mt452_mf1 pyne::endf::read_mt452_mf1(std::ifstream &infile) {
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

pyne::endf::mt455_mf1 pyne::endf::read_mt455_mf1(std::ifstream &infile) {
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

pyne::endf::mt456_mf1 pyne::endf::read_mt456_mf1(std::ifstream &infile) {
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

pyne::endf::mt458_mf1 pyne::endf::read_mt458_mf1(std::ifstream &infile) {
  mt458_mf1 mt458;
  control cs = read_cont(infile);
  mt458.nuc_id = cs.c1;
  mt458.awr = cs.c2;
  mt458.mat = cs.mat;
  mt458.mf = cs.mf;
  mt458.mt = cs.mt;

  list lst = read_list(infile);

  mt458.efr = std::vector<std::pair<double,double> > (lst.l2 + 1);
  mt458.pen = std::vector<std::pair<double,double> > (lst.l2 + 1);
  mt458.den = std::vector<std::pair<double,double> > (lst.l2 + 1);
  mt458.egp = std::vector<std::pair<double,double> > (lst.l2 + 1);
  mt458.egd = std::vector<std::pair<double,double> > (lst.l2 + 1);
  mt458.eb = std::vector<std::pair<double,double> > (lst.l2 + 1);
  mt458.enu = std::vector<std::pair<double,double> > (lst.l2 + 1);
  mt458.er = std::vector<std::pair<double,double> > (lst.l2 + 1);
  mt458.et = std::vector<std::pair<double,double> > (lst.l2 + 1);

  for (int i = 0; i < lst.l2 + 1; ++i) {
    mt458.efr[i] = std::make_pair(lst.data[i*18], lst.data[i*18+1]);
    mt458.pen[i] = std::make_pair(lst.data[i*18+2], lst.data[i*18+3]);
    mt458.den[i] = std::make_pair(lst.data[i*18+4], lst.data[i*18+5]);
    mt458.egp[i] = std::make_pair(lst.data[i*18+6], lst.data[i*18+7]);
    mt458.egd[i] = std::make_pair(lst.data[i*18+8], lst.data[i*18+9]);
    mt458.eb[i] = std::make_pair(lst.data[i*18+10], lst.data[i*18+11]);
    mt458.enu[i] = std::make_pair(lst.data[i*18+12], lst.data[i*18+13]);
    mt458.er[i] = std::make_pair(lst.data[i*18+14], lst.data[i*18+15]);
    mt458.et[i] = std::make_pair(lst.data[i*18+16], lst.data[i*18+17]);
  }

  std::string line;
  getline(infile, line);
  return mt458;
}

pyne::endf::mt460_mf1 pyne::endf::read_mt460_mf1(std::ifstream &infile) {
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

pyne::endf::mt457_mf8 pyne::endf::read_mt457_mf8(std::ifstream &infile) {
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
    mt457.nd = std::vector<int>(mt457.nsp);
    mt457.fd = std::vector<std::pair<double,double> >(mt457.nsp);
    mt457.eav = std::vector<std::pair<double,double> >(mt457.nsp);
    mt457.fc = std::vector<std::pair<double,double> >(mt457.nsp);
    for (int i = 0; i < mt457.nsp; ++i) {
      list spec = read_list(infile); //Basic spectrum information
      mt457.styp[i] = spec.c2;
      mt457.lcon[i] = spec.l1;
      mt457.ner[i] = spec.l2;
      mt457.nd[i] = spec.n2;
      mt457.fd[i] = std::make_pair(spec.data[0], spec.data[1]);
      mt457.eav[i] = std::make_pair(spec.data[2], spec.data[3]);
      mt457.fc[i] = std::make_pair(spec.data[4], spec.data[5]);
      if (spec.l1 != 1) { //LCON != 1
        for (int j = 0; j < spec.n2; ++j) {
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


pyne::endf::mf3 pyne::endf::read_mf3(std::ifstream &infile) {
  mf3 mf3_ob;
  control cs = read_cont(infile);
  mf3_ob.nuc_id = cs.c1;
  mf3_ob.awr = cs.c2;
  mf3_ob.mat = cs.mat;
  mf3_ob.mf = cs.mf;
  mf3_ob.mt = cs.mt;

  tab1 tab = read_tab1(infile);
  mf3_ob.nbt = tab.nbt;
  mf3_ob.intn = tab.intn;
  mf3_ob.energy = tab.x;
  mf3_ob.sigma = tab.y;

  std::string line;
  getline(infile, line);
  return mf3_ob;
}


pyne::endf::mf2 pyne::endf::read_mf2(std::ifstream &infile, int lrp) {
  mf2 mf2_ob;
  control cs = read_cont(infile);
  mf2_ob.nuc_id = cs.c1;
  mf2_ob.awr = cs.c2;
  mf2_ob.mat = cs.mat;
  mf2_ob.mf = cs.mf;
  mf2_ob.mt = cs.mt;

  mf2_ob.nis = cs.n1;
  tab1 tab;
  list lst;
  for(int i = 0; i < mf2_ob.nis; ++i) {
    cs = read_cont(infile);//isotope info
    mf2_ob.zai[i] = cs.c1;
    mf2_ob.abn[i] = cs.c2;
    mf2_ob.lfw[i] = cs.l2;
    mf2_ob.ner[i] = cs.n1;

    for(int j = 0; j < mf2_ob.ner; ++j) {
      cs = read_cont(infile);// range info
      mf2_ob.el[i][j] = cs.c1;
      mf2_ob.eh[i][j] = cs.c2;
      mf2_ob.lru[i][j] = cs.l1;
      mf2_ob.lrf[i][j] = cs.l2;
      mf2_ob.nro[i][j] = cs.n1;
      mf2_ob.naps[i][j] = cs.n2;
      if ((mf2_ob.lru[i][j] != 0) && (mf2_ob.nro[i][j] != 0)) {
        if (mf2_ob.lrf[i][j] < 5) {
          if ((mf2_ob.lru[i][j] == 1) && (mf2_ob.nro[i][j] != 0)) 
            tab = read_tab1(infile);
            mf2_ob.ap_nbt[i][j] = tab.nbt;
            mf2_ob.ap_intn[i][j] = tab.intn;
            mf2_ob.ap_e[i][j] = tab.x;
            mf2_ob.ap_val[i][j] = tab.y;
          if !((mf2_ob.lfw[i] == 1) && (mf2_ob.lrf[i][j] == 1)) {
            cs = read_cont(infile);// range info
            mf2_ob.spi[i][j] = cs.c1;
            mf2_ob.ap[i][j] = cs.c2;
            mf2_ob.lad[i][j] = cs.l1;
            mf2_ob.nls[i][j] = cs.n1;
            mf2_ob.nlsc[i][j] = cs.n2;
          }
          if (mf2_ob.lru[i][j] == 1) {
            if (mf2_ob.lrf[i][j] < 3) {
              for (int k = 0; k < mf2_ob.nls[i][j]; ++k) {
                lst = read_list(infile); //slbw and mlbw           
                mf2_ob.awri[i][j][k] = lst.c1;
                mf2_ob.qx[i][j][k] = lst.c2;
                mf2_ob.l[i][j][k] = lst.l1;
                mf2_ob.lrx[i][j][k] = lst.l2;
                for (int m = 0; m < lst.n2; ++m) {
                  mf2_ob.er[i][j][k][m] = lst.data[m*6];
                  mf2_ob.aj[i][j][k][m] = lst.data[m*6 + 1];
                  mf2_ob.gt[i][j][k][m] = lst.data[m*6 + 2];
                  mf2_ob.gn[i][j][k][m] = lst.data[m*6 + 3];
                  mf2_ob.gg[i][j][k][m] = lst.data[m*6 + 4];
                  mf2_ob.gf[i][j][k][m] = lst.data[m*6 + 5];
                }
              }
            } else if (mf2_ob.lrf[i][j] == 3) {
              for (int k = 0; k < mf2_ob.nls[i][j]; ++k) {
                lst = read_list(infile); // Reiche-Moore          
                mf2_ob.awri[i][j][k] = lst.c1;
                mf2_ob.apl[i][j][k] = lst.c2;
                mf2_ob.l[i][j][k] = lst.l1;
                for (int m = 0; m < lst.n2; ++m) {
                  mf2_ob.er[i][j][k][m] = lst.data[m*6];
                  mf2_ob.aj[i][j][k][m] = lst.data[m*6 + 1];
                  mf2_ob.gn[i][j][k][m] = lst.data[m*6 + 2];
                  mf2_ob.gg[i][j][k][m] = lst.data[m*6 + 3];
                  mf2_ob.gfa[i][j][k][m] = lst.data[m*6 + 4];
                  mf2_ob.gfb[i][j][k][m] = lst.data[m*6 + 5];
                }
              }
            } else {
              for (int k = 0; k < mf2_ob.nls[i][j]; ++k) {
                lst = read_list(infile); // Adler-Adler          
                mf2_ob.awri[i][j][k] = lst.c1;
                mf2_ob.li[i][j][k] = lst.l1;
                mf2_ob.abt[i][j][k] = std::vector<double>(&lst.data[0],&lst.data[6]);
                mf2_ob.abf[i][j][k] = std::vector<double>(&lst.data[6],&lst.data[12]);
                mf2_ob.abc[i][j][k] = std::vector<double>(&lst.data[12],&lst.data[18]);

                cs = read_cont(infile);
                mf2_ob.l[i][j][k] = cs.l1;
                mf2_ob.njs[i][j][k] = cs.n1;
                for (int m = 0; m < mf2_ob.njs[i][j]; ++m) {
                  lst = read_list(infile);
                  mf2_ob.aj[i][j][k][m] = lst.c1;
                  for (int n = 0; n < lst.n2; ++n) {
                    mf2_ob.det[i][j][k][m][n] = lst.data[n*12 + 0];
                    mf2_ob.dwt[i][j][k][m][n] = lst.data[n*12 + 1];
                    mf2_ob.grt[i][j][k][m][n] = lst.data[n*12 + 2];
                    mf2_ob.git[i][j][k][m][n] = lst.data[n*12 + 3];
                    mf2_ob.def[i][j][k][m][n] = lst.data[n*12 + 4];
                    mf2_ob.dwf[i][j][k][m][n] = lst.data[n*12 + 5];
                    mf2_ob.grf[i][j][k][m][n] = lst.data[n*12 + 6];
                    mf2_ob.gif[i][j][k][m][n] = lst.data[n*12 + 7];
                    mf2_ob.dec[i][j][k][m][n] = lst.data[n*12 + 8];
                    mf2_ob.dwc[i][j][k][m][n] = lst.data[n*12 + 9];
                    mf2_ob.grc[i][j][k][m][n] = lst.data[n*12 + 10];
                    mf2_ob.gic[i][j][k][m][n] = lst.data[n*12 + 11];
                  }

                }
              }
            }
          } else if ((mf2_ob.lfw[i] == 0) && (mf2_ob.lrf[i][j] == 1)) {
            for (int k = 0; k < mf2_ob.nls[i][j]; ++k) {
              lst = read_list(infile); // 
              mf2_ob.awri[i][j][k] = lst.c1;
              mf2_ob.l[i][j][k] = lst.l1;
              for (int m = 0; m < lst.n2; ++m) {
                mf2_ob.d[i][j][k][m] = lst.data[m*6];
                mf2_ob.aj[i][j][k][m] = lst.data[m*6 + 1];
                mf2_ob.amun[i][j][k][m] = lst.data[m*6 + 2];
                mf2_ob.gno[i][j][k][m] = lst.data[m*6 + 3];
                mf2_ob.gg[i][j][k][m] = lst.data[m*6 + 4];
              }
            }
          } else if ((mf2_ob.lfw[i] == 1) && (mf2_ob.lrf[i][j] == 1)) {
            lst = read_list(infile); // 
            mf2_ob.spi[i][j] = lst.c1;
            mf2_ob.ap[i][j] = lst.c2;
            mf2_ob.lssf[i][j] = lst.l1;
            mf2_ob.ne[i][j] = lst.n1;
            mf2_ob.nls[i][j] = lst.n2;
            mf2_ob.es[i][j] = lst.data;

            for (int k = 0; k < mf2_ob.nls[i][j]; ++k) {
              cs = read_cont(infile); 
              mf2_ob.awri[i][j][k] = cs.c1;
              mf2_ob.l[i][j][k] = cs.l1;
              mf2_ob.njs[i][j][k] = cs.n1;
              for (int m = 0; m < mf2_ob.njs[i][j][k]; ++m) {
                lst = read_list(infile);
                mf2_ob.muf[i][j][k][m] = lst.l2;
                mf2_ob.d[i][j][k][m] = lst.data[0];
                mf2_ob.aj[i][j][k][m] = lst.data[1];
                mf2_ob.amun[i][j][k][m] = lst.data[2];
                mf2_ob.gno[i][j][k][m] = lst.data[3];
                mf2_ob.gg[i][j][k][m] = lst.data[4];
                mf2_ob.gf[i][j][k][m] = std::vector<double>(&lst.data[6], &lst.data[lst.npl]);
              }
          } else if (mf2_ob.lrf[i][j] == 2) {
            for (int k = 0; k < mf2_ob.nls[i][j]; ++k) {
              cs = read_cont(infile);
              mf2_ob.awri[i][j][k] = cs.c1;
              mf2_ob.l[i][j][k] = cs.l1;
              mf2_ob.njs[i][j][k] = cs.n1;
              for (int m = 0; m < mf2_ob.njs[i][j][k]; ++m) {
                lst = read_list(infile);
                mf2_ob.muf[i][j][k][m] = lst.l2;
                mf2_ob.aj[i][j][k][m] = lst.c1;
                mf2_ob.amux[i][j][k][m] = lst.data[2];
                mf2_ob.amun[i][j][k][m] = lst.data[3];
                mf2_ob.amug[i][j][k][m] = lst.data[4];
                mf2_ob.amuf[i][j][k][m] = lst.data[5];
                for (int n = 6; n < lst.ne; ++n) {
                  mf2_ob.es[i][j][k][l][m][n] = lst.data(n*6 + 6);
                  mf2_ob.d[i][j][k][m][n] = lst.data(n*6 + 7);
                  mf2_ob.gx[i][j][k][m][n] = lst.data(n*6 + 8);
                  mf2_ob.gno[i][j][k][m][n] = lst.data(n*6 + 9);
                  mf2_ob.gg[i][j][k][m][n] = lst.data(n*6 + 10);
                  mf2_ob.gf[i][j][k][m][n] = lst.data(n*6 + 11);
                }
              }
            }
          }
        } else if (mf2_ob.lrf[i][j] == 7) {
          cs = read_cont(infile);

          mf2_ob.ifg[i][j] = cs.l1;
          mf2_ob.krm[i][j] = cs.l2;
          mf2_ob.njs[i][j] = cs.n1;
          mf2_ob.krl[i][j] = cs.n2;
          lst = read_list(infile);
          for (int k = 0; k < mf2_ob.njs[i][j]; ++k) {
            lst = read_list(infile);
            lst = read_list(infile);
            lst = read_list(infile);
            mf2_ob.lbk[i][j][k] = lst.npl;
            if (mf2_ob.lbk[i][j][k] == 1) {
              tab = read_tab1(infile);
              tab = read_tab1(infile);
            } else if (mf2_ob.lbk[i][j][k] == 2) {
              lst = read_list(infile);
            } else if (mf2_ob.lbk[i][j][k] == 3) {
              lst = read_list(infile);
            }
            lst = read_list(infile);
            mf2_ob.lps[i][j][k] = lst.npl;
            if (mf2_ob.lps[i][j][k] == 1) {
              tab = read_tab1(infile);
              tab = read_tab1(infile);
            }
          }
        }
      } else {
        cs = read_cont(infile);// range info
        mf2_ob.spi[i][j] = cs.c1;
        mf2_ob.ap[i][j] = cs.c2;
        mf2_ob.nls[i][j] = cs.n1;

      }
    }
  }
  std::string line;
  getline(infile, line);
  return mf2_ob;
}


