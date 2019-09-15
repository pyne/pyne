// Enrichment
#ifndef PYNE_IS_AMALGAMATED
#include "enrichment.h"
#endif

namespace pyne_enr = pyne::enrichment;

pyne_enr::Cascade pyne_enr::_fill_default_uranium_cascade() {
  // Default cascade for uranium-based enrichment
  Cascade duc;

  duc.alpha = 1.05;
  duc.Mstar = 236.5;

  duc.j = 922350000;
  duc.k = 922380000;

  duc.N = 30.0;
  duc.M = 10.0;

  duc.x_feed_j = 0.0072;
  duc.x_prod_j = 0.05;
  duc.x_tail_j = 0.0025;

  pyne::comp_map cm;
  cm[922340000] = 0.000055;
  cm[922350000] = 0.00720;
  cm[922380000] = 0.992745;
  duc.mat_feed = pyne::Material(cm, 1.0, 1.0);

  return duc;
}
pyne_enr::Cascade pyne_enr::default_uranium_cascade(pyne_enr::_fill_default_uranium_cascade());



double pyne_enr::feed_per_prod(double x_feed, double x_prod, double x_tail) {
  return 1 / prod_per_feed(x_feed, x_prod, x_tail);
}

double pyne_enr::feed_per_tail(double x_feed, double x_prod, double x_tail) {
  return 1 / tail_per_feed(x_feed, x_prod, x_tail);
}

double pyne_enr::prod_per_tail(double x_feed, double x_prod, double x_tail) {
  return 1 / tail_per_prod(x_feed, x_prod, x_tail);
}

double pyne_enr::prod_per_feed(double x_feed, double x_prod, double x_tail) {
  return (x_feed - x_tail) / (x_prod - x_tail);
}

double pyne_enr::tail_per_feed(double x_feed, double x_prod, double x_tail) {
  return (x_feed - x_prod) / (x_tail - x_prod);
}

double pyne_enr::tail_per_prod(double x_feed, double x_prod, double x_tail) {
  return (x_feed - x_prod) / (x_tail - x_feed);
}

double pyne_enr::value_func(double x) {
  return (2 * x - 1) * log(x / (1 - x));
}

double pyne_enr::swu_per_feed(double x_feed, double x_prod, double x_tail) {
  return
      prod_per_feed(x_feed, x_prod, x_tail) * value_func(x_prod) + \
      tail_per_feed(x_feed, x_prod, x_tail) * value_func(x_tail) - \
      value_func(x_feed);
}

double pyne_enr::swu_per_prod(double x_feed, double x_prod, double x_tail) {
  return
      value_func(x_prod) + \
      tail_per_prod(x_feed, x_prod, x_tail) * value_func(x_tail) - \
      feed_per_prod(x_feed, x_prod, x_tail) * value_func(x_feed);
}

double pyne_enr::swu_per_tail(double x_feed, double x_prod, double x_tail) {
  return
      prod_per_tail(x_feed, x_prod, x_tail) * value_func(x_prod) + \
      value_func(x_tail) - \
      feed_per_tail(x_feed, x_prod, x_tail) * value_func(x_feed);
}


double pyne_enr::alphastar_i(double alpha, double Mstar, double M_i) {
  // M_i is the mass of the ith nuclide
  return pow(alpha, (Mstar - M_i));
}


void pyne_enr::_recompute_nm(pyne_enr::Cascade & casc, double tolerance) {

  double ppf = prod_per_feed(casc.mat_feed.comp[casc.j], casc.x_prod_j, casc.x_tail_j);
  double tpf = tail_per_feed(casc.mat_feed.comp[casc.j], casc.x_prod_j, casc.x_tail_j);
  double astar_j = alphastar_i(casc.alpha, casc.Mstar, pyne::atomic_mass(casc.j));

  // Save original state of N & M
  double N = casc.N;
  double M = casc.M;
  double origN = casc.N;
  double origM = casc.M;

  double lhs_prod = ppf * casc.x_prod_j / casc.mat_feed.comp[casc.j];
  double rhs_prod = (pow(astar_j, M+1.0) - 1.0) / (pow(astar_j, M+1.0) - pow(astar_j, -N));
  double lhs_tail = tpf * casc.x_tail_j / casc.mat_feed.comp[casc.j];
  double rhs_tail = (1.0 - pow(astar_j, -N)) / (pow(astar_j, M+1.0) - pow(astar_j, -N));

  double n = 1.0;
  double delta_prod = lhs_prod - rhs_prod;
  double delta_tail = lhs_tail - rhs_tail;
  while (tolerance < fabs(delta_prod) && tolerance < fabs(delta_tail)) {
    delta_prod = lhs_prod - rhs_prod;
    delta_tail = lhs_tail - rhs_tail;

    if (tolerance < fabs(delta_prod)) {
      N = N - (delta_prod * N);
      rhs_prod = (pow(astar_j, M+1.0) - 1.0) / (pow(astar_j, M+1.0) - pow(astar_j, -N));
    }

    if (tolerance < fabs(delta_tail)) {
      M = M - (delta_tail * M);
      rhs_tail = (1.0 - pow(astar_j, -N)) / (pow(astar_j, M+1.0) - pow(astar_j, -N));
    }

    if (N < tolerance) {
      N = origN + n;
      M = origM + n;
      n = n + 1.0;
    }

    if (M < tolerance) {
      N = origN + n;
      M = origM + n;
      n = n + 1.0;
    }
  }

  casc.N = N;
  casc.M = M;
  return;
}



void pyne_enr::_recompute_prod_tail_mats(pyne_enr::Cascade & casc) {
  //This function takes a given initial guess number of enriching and stripping stages
  //for a given composition of fuel with a given jth key component, knowing the values
  //that are desired in both Product and Tails streams.  Having this it solves for what
  //the actual N and M stage numbers are and also what the product and waste streams
  //compositions are.  It returns precisely these.
  pyne::comp_map comp_prod;
  pyne::comp_map comp_tail;

  int nuc;
  double astar_i, numer_prod, numer_tail, denom_prod, denom_tail;

  double N = casc.N;
  double M = casc.M;

  double ppf = prod_per_feed(casc.mat_feed.comp[casc.j], casc.x_prod_j, casc.x_tail_j);
  double tpf = tail_per_feed(casc.mat_feed.comp[casc.j], casc.x_prod_j, casc.x_tail_j);

  for (pyne::comp_iter i = casc.mat_feed.comp.begin(); i != casc.mat_feed.comp.end(); i++) {
    nuc = (i->first);
    astar_i = alphastar_i(casc.alpha, casc.Mstar, pyne::atomic_mass(nuc));

    // calc prod comp
    numer_prod = casc.mat_feed.comp[nuc] * (pow(astar_i, M+1.0) - 1.0);
    denom_prod = (pow(astar_i, M+1.0) - pow(astar_i, -N)) / ppf;
    comp_prod[nuc] = numer_prod / denom_prod;

    // calc tail comp
    numer_tail = casc.mat_feed.comp[nuc] * (1.0 - pow(astar_i, -N));
    denom_tail = (pow(astar_i, M+1.0) - pow(astar_i, -N)) / tpf;
    comp_tail[nuc] = numer_tail / denom_tail;
  }

  casc.mat_prod = pyne::Material(comp_prod);
  casc.mat_tail = pyne::Material(comp_tail);
  return;
}



pyne_enr::Cascade pyne_enr::_norm_comp_secant(pyne_enr::Cascade & casc, \
                                              double tolerance, int max_iter) {
  // This function actually solves the whole system of equations.  It uses _recompute_prod_tail_mats
  // to find the roots for the enriching and stripping stage numbers.  It then
  // checks to see if the product and waste streams meet their target enrichments
  // for the jth component like they should.  If they don't then it trys other values
  // of N and M varied by the Secant ethod.  Rinse and repeat as needed.
  int j = casc.j;
  pyne_enr::Cascade prev_casc = casc;
  pyne_enr::Cascade curr_casc = casc;

  // Is the history of N and M that has been input
  unsigned int h;
  int niter = 0;
  int max_hist = max_iter / 10;
  std::vector<double> historyN;
  std::vector<double> historyM;

  // Initialize prev point
  prev_casc.N += 1.0;
  prev_casc.M += 1.0;
  _recompute_nm(prev_casc, tolerance);
  _recompute_prod_tail_mats(prev_casc);
  historyN.push_back(prev_casc.N);
  historyM.push_back(prev_casc.M);


  // Initialize current point
  _recompute_nm(curr_casc, tolerance);
  _recompute_prod_tail_mats(curr_casc);
  historyN.push_back(curr_casc.N);
  historyM.push_back(curr_casc.M);

  // My guess is that what we are checkin here is that the isotopic compositions
  // make sense with abs(1.0 - masscurr_P) rather than calculatign the
  // relative product to watse mass streams.
  double prev_N = prev_casc.N;
  double prev_M = prev_casc.M;
  double curr_N = curr_casc.N;
  double curr_M = curr_casc.M;
  double temp_prev_N = 0.0;
  double temp_prev_M = 0.0;
  double temp_curr_N = 0.0;
  double temp_curr_M = 0.0;

  double delta_x_prod_j = casc.x_prod_j - curr_casc.mat_prod.comp[j];
  double delta_x_tail_j = casc.x_tail_j - curr_casc.mat_tail.comp[j];

  while ((tolerance < fabs(delta_x_prod_j) / curr_casc.mat_prod.comp[j]  || \
          tolerance < fabs(delta_x_tail_j) / curr_casc.mat_tail.comp[j]) && \
          niter < max_iter) {
    delta_x_prod_j = casc.x_prod_j - curr_casc.mat_prod.comp[j];
    delta_x_tail_j = casc.x_tail_j - curr_casc.mat_tail.comp[j];

    if (tolerance <= fabs(delta_x_prod_j)/curr_casc.mat_prod.comp[j]) {
      // Make a new guess for N
      temp_curr_N = curr_N;
      temp_prev_N = prev_N;
      curr_N = curr_N + delta_x_prod_j*\
              ((curr_N - prev_N)/(curr_casc.mat_prod.comp[j] - prev_casc.mat_prod.comp[j]));
      prev_N = temp_curr_N;

      // If the new value of N is less than zero, reset.
      if (curr_N < 0.0)
        curr_N = (temp_curr_N + temp_prev_N)/2.0;
    }

    if (tolerance <= fabs(delta_x_tail_j)/curr_casc.mat_tail.comp[j]) {
      // Make a new guess for M
      temp_curr_M = curr_M;
      temp_prev_M = prev_M;
      curr_M = curr_M + delta_x_tail_j*\
               ((curr_M - prev_M)/(curr_casc.mat_tail.comp[j] - prev_casc.mat_tail.comp[j]));
      prev_M = temp_curr_M;

      // If the new value of M is less than zero, reset.
      if (curr_M < 0.0)
        curr_M = (temp_curr_M + temp_prev_M)/2.0;
    }

    // Check for infinite loops
    for (h = 0; h < historyN.size(); h++) {
      if (historyN[h] == curr_N && historyM[h] == curr_M) {
        curr_N = curr_N + delta_x_prod_j * \
              ((curr_N - prev_N)/(curr_casc.mat_prod.comp[j] - prev_casc.mat_prod.comp[j]));
        curr_M = curr_M + delta_x_tail_j * \
               ((curr_M - prev_M)/(curr_casc.mat_tail.comp[j] - prev_casc.mat_tail.comp[j]));;
        break;
      }
    }

    if (max_hist <= historyN.size()) {
      historyN.erase(historyN.begin());
      historyM.erase(historyM.begin());
    }
    historyN.push_back(curr_N);
    historyM.push_back(curr_M);

    niter += 1;

    // Calculate new isotopics for valid (N, M)
    prev_casc = curr_casc;
    curr_casc.N = curr_N;
    curr_casc.M = curr_M;
    _recompute_nm(curr_casc, tolerance);
    _recompute_prod_tail_mats(curr_casc);
  }
  return curr_casc;
}



double pyne_enr::_deltaU_i_OverG(pyne_enr::Cascade & casc, int i) {
  // Solves for a stage separative power relevant to the ith component
  // per unit of flow G.  This is taken from Equation 31 divided by G
  // from the paper "Wood, Houston G., Borisevich, V. D. and Sulaberidze, G. A.,
  // 'On a Criterion Efficiency for Multi-Isotope Mixtures Separation',
  // Separation Science and Technology, 34:3, 343 - 357"
  // To link to this article: DOI: 10.1081/SS-100100654
  // URL: http://dx.doi.org/10.1081/SS-100100654

  double astar_i = alphastar_i(casc.alpha, casc.Mstar, pyne::atomic_mass(i));
  return log(pow(casc.alpha, (casc.Mstar - pyne::atomic_mass(casc.j)) )) * \
                             ((astar_i - 1.0)/(astar_i + 1.0));
}


pyne_enr::Cascade pyne_enr::solve_numeric(pyne_enr::Cascade & orig_casc, \
                                          double tolerance, int max_iter) {
  // This function finds the total flow rate (L) over the feed flow rate (F)
  pyne_enr::Cascade casc = orig_casc;

  casc = _norm_comp_secant(casc, tolerance, max_iter);

  int nuc;
  int j = casc.j;
  int k = casc.k;
  double ppf = prod_per_feed(casc.mat_feed.comp[j], casc.x_prod_j, casc.x_tail_j);
  double tpf = tail_per_feed(casc.mat_feed.comp[j], casc.x_prod_j, casc.x_tail_j);

  // Matched Flow Ratios
  double rfeed = casc.mat_feed.comp[j] / casc.mat_feed.comp[k];
  double rprod = casc.mat_prod.comp[j] / casc.mat_prod.comp[k];
  double rtail = casc.mat_tail.comp[j] / casc.mat_tail.comp[k];

  double ltotpf = 0.0;
  double swupf = 0.0;
  double temp_numer = 0.0;

  for (pyne::comp_iter i = casc.mat_feed.comp.begin(); i != casc.mat_feed.comp.end(); i++) {
    nuc = (i->first);
    temp_numer = (ppf*casc.mat_prod.comp[nuc]*log(rprod) + \
                  tpf*casc.mat_tail.comp[nuc]*log(rtail) - \
                      casc.mat_feed.comp[nuc]*log(rfeed));
    ltotpf = ltotpf + (temp_numer / _deltaU_i_OverG(casc, nuc));
    swupf = swupf + temp_numer;
  }

  // Assign flow rates
  casc.l_t_per_feed = ltotpf;

  // The -1 term is put in the SWU calculation because otherwise swupf
  // represents the SWU that would be undone if you were to deenrich the
  // whole process.  Thus the SWU to enrich is -1x this number.  This is
  // a by-product of the value function used as a constraint.
  casc.swu_per_feed = -1 * swupf;       // This is the SWU for 1 kg of Feed material.
  casc.swu_per_prod = -1 * swupf / ppf;	// This is the SWU for 1 kg of Product material.

  // Assign isotopic streams the proper masses.
  casc.mat_prod.mass = casc.mat_feed.mass * ppf;
  casc.mat_tail.mass = casc.mat_feed.mass * tpf;
  return casc;
}

pyne_enr::Cascade pyne_enr::multicomponent(pyne_enr::Cascade & orig_casc, \
                                    char * solver, double tolerance, int max_iter) {
  std::string strsolver(solver);
  return multicomponent(orig_casc, strsolver, tolerance, max_iter);
}

pyne_enr::Cascade pyne_enr::multicomponent(pyne_enr::Cascade & orig_casc, \
                                    std::string solver, double tolerance, int max_iter) {
  // The multicomponent() function finds a value of Mstar by minimzing the seperative power.
  // Note that Mstar0 represents an intial guess at what Mstar might be.
  // This is the final function that actually solves for an optimized M* that makes the cascade!
  pyne_enr::Cascade temp_casc;
  pyne_enr::Cascade prev_casc = orig_casc;
  pyne_enr::Cascade curr_casc = orig_casc;

  // define the solver to use
  int solver_code;
  if (solver == "symbolic")
    solver_code = 0;
  else if (solver == "numeric")
    solver_code = 1;
  else
    throw "solver not known: " + solver;

  // validate Mstar or pick new value
  if ((orig_casc.Mstar < pyne::atomic_mass(orig_casc.j) &&  \
       orig_casc.Mstar < pyne::atomic_mass(orig_casc.k)) || \
      (orig_casc.Mstar > pyne::atomic_mass(orig_casc.j) &&  \
       orig_casc.Mstar > pyne::atomic_mass(orig_casc.k))) {
    double ms = (pyne::atomic_mass(orig_casc.j) + pyne::atomic_mass(orig_casc.k)) / 2.0;
    prev_casc.Mstar = ms;
    curr_casc.Mstar = ms;
  }

  // xpn is the exponential index
  double ooe = -log10(tolerance);
  double xpn = 1.0;

  // Initialize previous point
  switch (solver_code) {
    case 0:
      prev_casc = solve_symbolic(prev_casc);
      break;
    case 1:
      prev_casc = solve_numeric(prev_casc, tolerance, max_iter);
      break;
  }

  // Initialize curr_ent point
  curr_casc.Mstar = (pyne::atomic_mass(curr_casc.j) + curr_casc.Mstar) / 2.0;
  switch (solver_code) {
    case 0:
      curr_casc = solve_symbolic(curr_casc);
      break;
    case 1:
      curr_casc = solve_numeric(curr_casc, tolerance, max_iter);
      break;
  }

  double m = pyne::slope(curr_casc.Mstar, curr_casc.l_t_per_feed, \
                         prev_casc.Mstar, prev_casc.l_t_per_feed);
  double m_sign = m / fabs(m);

  double temp_m;
  double temp_m_sign;

  while (tolerance < fabs(curr_casc.l_t_per_feed - prev_casc.l_t_per_feed) / curr_casc.l_t_per_feed) {
    // Check that parameters are still well-formed
    if (isnan(curr_casc.Mstar) || isnan(curr_casc.l_t_per_feed) || \
        isnan(prev_casc.Mstar) || isnan(prev_casc.l_t_per_feed))
      throw EnrichmentIterationNaN();

    prev_casc = curr_casc;

    curr_casc.Mstar = curr_casc.Mstar - (m_sign * pow(10.0, -xpn));
    switch (solver_code) {
      case 0:
        curr_casc = solve_symbolic(curr_casc);
        break;
      case 1:
        curr_casc = solve_numeric(curr_casc, tolerance, max_iter);
        break;
    }

    if (prev_casc.l_t_per_feed < curr_casc.l_t_per_feed) {
      temp_casc = curr_casc;
      temp_casc.Mstar = temp_casc.Mstar - (m_sign * pow(10.0, -xpn));
      switch (solver_code) {
        case 0:
          temp_casc = solve_symbolic(temp_casc);
          break;
        case 1:
          temp_casc = solve_numeric(temp_casc, tolerance, max_iter);
          break;
      }

      temp_m = pyne::slope(curr_casc.Mstar, curr_casc.l_t_per_feed, \
                           temp_casc.Mstar, temp_casc.l_t_per_feed);
      if (temp_m == 0.0) {
        prev_casc = curr_casc;
        curr_casc = temp_casc;
        break;
      }

      temp_m_sign = temp_m / fabs(temp_m);
      if (m_sign != temp_m_sign) {
        xpn = xpn + 1;

        temp_casc = prev_casc;
        temp_casc.Mstar = temp_casc.Mstar + (m_sign * pow(10.0, -xpn));
        switch (solver_code) {
          case 0:
            temp_casc = solve_symbolic(temp_casc);
            break;
          case 1:
            temp_casc = solve_numeric(temp_casc, tolerance, max_iter);
            break;
        }
        temp_m = pyne::slope(prev_casc.Mstar, prev_casc.l_t_per_feed, \
                             temp_casc.Mstar, temp_casc.l_t_per_feed);

        if (temp_m == 0.0) {
          prev_casc = curr_casc;
          curr_casc = temp_casc;
          break;
        }

        m_sign = temp_m / fabs(temp_m);
        m = temp_m;
        prev_casc = curr_casc;
        curr_casc = temp_casc;
      }
    }
  }

  return curr_casc;
}
