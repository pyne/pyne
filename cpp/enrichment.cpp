// Enrichment 

#include "enrichment.h"

namespace pyne_enr = pyne::enrichment;

/*********************************/
/*** Enrichment Helper Classes ***/
/*********************************/
pyne_enr::Cascade::Cascade()
{
  alpha = 0.0;
  Mstar = 0.0;

  j = 0;
  k = 0;

  N = 0.0;
  M = 0.0;

  x_feed_j = 0.0;
  x_prod_j = 0.0;
  x_tail_j = 0.0;

  mat_feed = pyne::Material();
  mat_prod = pyne::Material();
  mat_tail = pyne::Material();
};


pyne_enr::Cascade::~Cascade()
{
};


void pyne_enr::Cascade::_reset_xjs()
{
  // resets the key enriment member variables
  x_feed_j = mat_feed.comp[j];
  x_prod_j = mat_prod.comp[j];
  x_tail_j = mat_tail.comp[j];
};

pyne_enr::Cascade pyne_enr::_fill_default_uranium_cascade()
{
  // Default cascade for uranium-based enrichment
  Cascade duc; 

  duc.alpha = 1.05;
  duc.Mstar = 236.5;

  duc.j = 922350;
  duc.k = 922380;

  duc.N = 30.0;
  duc.M = 10.0;

  duc.x_feed_j = 0.00711;
  duc.x_prod_j = 0.05;
  duc.x_tail_j = 0.0025;

  return duc;
};
pyne_enr::Cascade pyne_enr::default_uranium_cascade(pyne_enr::_fill_default_uranium_cascade());



double pyne_enr::prod_per_feed(double x_feed, double x_prod, double x_tail)
{
  // Product over Feed Enrichment Ratio
  return ((x_feed - x_tail)/(x_prod - x_tail));
}

double pyne_enr::tail_per_feed(double x_feed, double x_prod, double x_tail)
{
  // Tails over Feed Enrichment Ratio
  return ((x_feed - x_prod)/(x_tail - x_prod));
}

double pyne_enr::tail_per_prod(double x_feed, double x_prod, double x_tail)
{
  // Tails over Feed Enrichment Ratio
  return ((x_feed - x_prod)/(x_tail - x_feed));
}


double pyne_enr::alphastar_i(double alpha, double Mstar, double M_i)
{
  // M_i is the mass of the ith nuclide
  return pow(alpha, (Mstar - M_i));
}

double pyne_enr::Ei(double astar_i, double N)
{
  return ((astar_i - 1.0) / (1.0 - pow(astar_i, -N)));
};


double pyne_enr::Si(double astar_i, double M)
{
  return ((astar_i - 1.0)/(pow(astar_i, M+1) - 1.0));
};

void pyne_enr::_recompute_nm(pyne_enr::Cascade & casc, double tolerance)
{

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
  while (tolerance < fabs(delta_prod) && tolerance < fabs(delta_tail))
  {
    delta_prod = lhs_prod - rhs_prod;
    delta_tail = lhs_tail - rhs_tail;

    if (tolerance < fabs(delta_prod))
    {
      N = N - (delta_prod * N);
      rhs_prod = (pow(astar_j, M+1.0) - 1.0) / (pow(astar_j, M+1.0) - pow(astar_j, -N));
    };

    if (tolerance < fabs(delta_tail))
    {
      M = M - (delta_tail * M);
      rhs_tail = (1.0 - pow(astar_j, -N)) / (pow(astar_j, M+1.0) - pow(astar_j, -N));
    };

    if (N < tolerance)
    {
      N = origN + n;
      M = origM + n;
      n = n + 1.0;
    };

    if (M < tolerance)
    {
      N = origN + n;
      M = origM + n;
      n = n + 1.0;
    };
  };

  casc.N = N;
  casc.M = M;
  return; 
};



void pyne_enr::_recompute_prod_tail_mats(pyne_enr::Cascade & casc)
{
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
  double x_feed_j = casc.mat_feed.comp[casc.j];
  double x_prod_j = casc.x_prod_j;
  double x_tail_j = casc.x_tail_j;

  for (pyne::comp_iter i = casc.mat_feed.comp.begin(); i != casc.mat_feed.comp.end(); i++)
  {
    nuc = (i->first);
    astar_i = alphastar_i(casc.alpha, casc.Mstar, pyne::atomic_mass(nuc));

    // calc prod comp
    numer_prod = casc.mat_feed.comp[nuc] * (pow(astar_i, M+1.0) - 1.0);
    denom_prod = (pow(astar_i, M+1.0) - pow(astar_i, -N)) / \
                  prod_per_feed(x_feed_j, x_prod_j, x_tail_j);
    comp_prod[nuc] = numer_prod / denom_prod;

    // calc tail comp
    numer_tail = casc.mat_feed.comp[nuc] * (1.0 - pow(astar_i, -N));
	  denom_tail = (pow(astar_i, M+1.0) - pow(astar_i, -N)) / \
                  tail_per_feed(x_feed_j, x_prod_j, x_tail_j);
    comp_tail[nuc] = numer_tail / denom_tail;
  };

  casc.mat_prod = pyne::Material(comp_prod);
  casc.mat_prod.mass *= casc.mat_feed.mass;
  casc.x_prod_j = casc.mat_prod.comp[casc.j];

  casc.mat_tail = pyne::Material(comp_tail);
  casc.mat_tail.mass *= casc.mat_feed.mass;
  casc.x_tail_j = casc.mat_tail.comp[casc.j];

  return;
};



pyne_enr::Cascade pyne_enr::_norm_comp_secant(pyne_enr::Cascade & casc, double tolerance)
{
  // This function actually solves the whole system of equations.  It uses _recompute_prod_tail_mats 
  // to find the roots for the enriching and stripping stage numbers.  It then 
  // checks to see if the product and waste streams meet their target enrichments
  // for the jth component like they should.  If they don't then it trys other values 
  // of N and M varied by Newton's Method.  Rinse and repeat as needed.
  pyne_enr::Cascade prev_casc = casc;
  pyne_enr::Cascade curr_casc = casc;

  // Is the history of N and M that has been input
  uint h;
  std::vector<double> historyN;
  std::vector<double> historyM;

  // Start iteration Counter
  int counter = 0;

  // Initialize prev point
  prev_casc.N += 1.0;
  prev_casc.M += 1.0;
  _recompute_nm(prev_casc, tolerance);
  _recompute_prod_tail_mats(prev_casc);
  historyN.push_back(prev_casc.N);
  historyN.push_back(prev_casc.M);


  // Initialize current point
  _recompute_nm(curr_casc, tolerance);
  _recompute_prod_tail_mats(curr_casc);
  historyN.push_back(curr_casc.N);
  historyN.push_back(curr_casc.M);

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

  double delta_x_prod_j = casc.x_prod_j - curr_casc.x_prod_j;
  double delta_x_tail_j = casc.x_tail_j - curr_casc.x_tail_j;

  while (tolerance < fabs(delta_x_prod_j) || tolerance < fabs(delta_x_tail_j))
  {
    delta_x_prod_j = casc.x_prod_j - curr_casc.x_prod_j;
    delta_x_tail_j = casc.x_tail_j - curr_casc.x_tail_j;

    if (tolerance <= fabs(delta_x_prod_j))
    {
      // Make a new guess for N
      temp_curr_N = curr_N;
      temp_prev_N = prev_N;
      curr_N = curr_N + delta_x_prod_j*\
              ((curr_N - prev_N)/(curr_casc.x_prod_j - prev_casc.x_prod_j));
      prev_N = temp_curr_N;

      // If the new value of N is less than zero, reset.
      if (curr_N < 0.0)
        curr_N = (temp_curr_N + temp_prev_N)/2.0;
    };

    if (tolerance <= fabs(delta_x_tail_j))
    {
      // Make a new guess for M
      temp_curr_M = curr_M;
      temp_prev_M = prev_M;
      curr_M = curr_M + delta_x_tail_j*\
               ((curr_M - prev_M)/(curr_casc.x_tail_j - prev_casc.x_tail_j));
      prev_M = temp_curr_M;

      // If the new value of M is less than zero, reset.
      if (curr_M < 0.0)
        curr_M = (temp_curr_M + temp_prev_M)/2.0;
    };

    // Check for infinite loops
    for (h = 0; h < historyN.size(); h++)
      if (historyN[h] == curr_N && historyM[h] == curr_M)
        throw EnrichmentInfiniteLoopError();

    if (150 <= historyN.size())
    {
      historyN.erase(historyN.begin());
      historyM.erase(historyM.begin());
    };
    historyN.push_back(curr_N);
    historyM.push_back(curr_M);

    if (10000 < counter)
      throw EnrichmentIterationLimit();
    else
      counter += 1;

    // Calculate new isotopics for valid (N, M)        
    prev_casc = curr_casc;
    curr_casc.N = curr_N;
    curr_casc.M = curr_M;
    _recompute_nm(curr_casc, tolerance);
    _recompute_prod_tail_mats(curr_casc);
  };

  return curr_casc;
};


// I have serious doubts that this works...
pyne_enr::Cascade pyne_enr::_norm_comp_other(pyne_enr::Cascade & orig_casc, double tolerance)
{
  pyne_enr::Cascade casc = orig_casc;

  // Is the history of N and M that has been input
  uint h;
  std::vector<double> historyN;
  std::vector<double> historyM;

  // Initial point
  double N = casc.N;
  double M = casc.M;
  _recompute_nm(casc, tolerance);
  _recompute_prod_tail_mats(casc);
  double mass_prod = casc.mat_prod.mass;
  double mass_tail = casc.mat_tail.mass;

  while (tolerance < fabs(1.0 - mass_prod) && tolerance < fabs(1.0 - mass_tail) )
  {
    if (tolerance <= fabs(1.0 - mass_prod))
      N = N - (1.0 - mass_prod)/(1.0 + mass_prod);

    if (tolerance <= fabs(1.0 - mass_tail))
      M = M + (1.0 - mass_tail)/(1.0 + mass_tail);

    // Note this infinite loop checker does not raise an exception
    // Thus the exception cannot be caught by a try statement and then another
    // root finding Comp2Unity method tried.
    // This simply leaves N and M in whatever their curr_ent state is at the time.
    // This is fine for M* = 235.1 for U, since M* won't be optimized here as an outlier 
    // and since it is looping it is probably signalling around some actual value.

    // Check for infinite loops
    for (h = 0; h < historyN.size(); h++)
      if (historyN[h] == N && historyM[h] == M)
        //throw EnrichmentInfiniteLoopError(); //Possible future use.
        return casc;

    if (150 <= historyN.size())
    {
      historyN.erase(historyN.begin());
      historyM.erase(historyM.begin());
    };
    historyN.push_back(N);
    historyM.push_back(M);

    // Calculate new masses
    casc.N = N;
    casc.M = M;
    _recompute_nm(casc, tolerance);
    _recompute_prod_tail_mats(casc);
    mass_prod = casc.mat_prod.mass;
    mass_tail = casc.mat_tail.mass;
  };

  return casc;
};

/*

double pyne_enr::deltaU_i_OverG(int i)
{
  // Solves for a stage separative power relevant to the ith component
  // per unit of flow G.  This is taken from Equation 31 divided by G 
  // from the paper "Wood, Houston G., Borisevich, V. D. and Sulaberidze, G. A.,
  // 'On a Criterion Efficiency for Multi-Isotope Mixtures Separation', 
  // Separation Science and Technology, 34:3, 343 - 357"
  // To link to this article: DOI: 10.1081/SS-100100654
  // URL: http://dx.doi.org/10.1081/SS-100100654

  double alphastar_i = alphastar_i(pyne::atomic_mass(i));
  return log(pow( alpha_0, (Mstar - pyne::atomic_mass(j)) )) * ((alphastar_i - 1.0)/(alphastar_i + 1.0));
};


void pyne_enr::ltot_per_feed()
{
  // This function finds the total flow rate (L) over the feed flow rate (F)
  bool comp_converged = false;

  try
  {
    // Try secant method first
    _norm_comp_secant();
    comp_converged = true;
  }
  catch (...)
  {
    try
    {
      // Then try other cr8zy method
      _norm_comp_other();
    	comp_converged = true;
    }
		catch (...)
    {
      // No other methods to try!
      comp_converged = false;
    };
  };

  if (comp_converged)
  {
    double PoF = prod_per_feed(mat_feed.comp[j], x_prod_j, x_tail_j);
    double WoF = tail_per_feed(mat_feed.comp[j], x_prod_j, x_tail_j);

    // Matched Flow Ratios
    double RF = mat_feed.comp[j] / mat_feed.comp[k];
    double RP = mat_prod.comp[j] / mat_prod.comp[k];
    double RW = mat_tail.comp[j] / mat_tail.comp[k];

    double LtotalOverF = 0.0;
    double SWUoverF = 0.0;
    double tempNumerator = 0.0; 

    for (pyne::comp_iter i = mat_feed.comp.begin(); i != mat_feed.comp.end(); i++)
    {
      tempNumerator = (PoF*mat_prod.comp[i->first]*log(RP) + WoF*mat_tail.comp[i->first]*log(RW) - mat_feed.comp[i->first]*log(RF));
      LtotalOverF = LtotalOverF + (tempNumerator / deltaU_i_OverG(i->first));
      SWUoverF = SWUoverF + tempNumerator;
    };

    if (0 < bright::verbosity)
      std::cout << "    L/F = " << LtotalOverF << "\n";        

    // Assign flow rates
    TotalPerFeed = LtotalOverF;

    // The -1 term is put in the SWU calculation because otherwise SWUoverF   
    // represents the SWU that would be undone if you were to deenrich the 
    // whole process.  Thus the SWU to enrich is -1x this number.  This is 
    // a by-product of the value function used as a constraint.
    swu_per_feed    = -1 * SWUoverF;          //This is the SWU for 1 kg of Feed material.
    swu_per_prod = -1 * SWUoverF / PoF;	//This is the SWU for 1 kg of Product material.

    // Assign Isotopic streams the proper masses.
    mat_prod.mass  = mat_feed.mass * PoF;
    mat_tail.mass = mat_feed.mass * WoF;
  };

  return;
};


void pyne_enr::multicomponent(double Mstar_0, double tolerance)
{
  // The multicomponent() function finds a value of Mstar by minimzing the seperative power.  
  // Note that Mstar0 represents an intial guess at what Mstar might be.
  // This is the final function that actually solves for an optimized M* that makes the cascade!

  // History table that has Mstar, LoF, and slope between this point and the last_ one 
  // hist = []


  // xpn is the exponential index index
  double xpn = 1.0;

  // Initialize previous point
  double prev_Mstar = Mstar_0;
  double prev_ltot_per_feed = ltot_per_feed();

  // Initialize curr_ent point
  double curr__Mstar = Mstar_0 + 0.1;
  double curr__ltot_per_feed = ltot_per_feed();

  double m = pyne::slope(curr__Mstar, curr__ltot_per_feed, prev_Mstar, prev_ltot_per_feed);
  double m_sign = m / fabs(m);

  double tempMstar;
  double templtot_per_feed;
  double tempm;
  double tempm_sign;

  if (0.0 < m_sign)
  {
    tempMstar  = prev_Mstar;
    templtot_per_feed = prev_ltot_per_feed;
    prev_Mstar  = curr__Mstar;
    prev_ltot_per_feed = curr__ltot_per_feed;
    curr__Mstar  = tempMstar;
    curr__ltot_per_feed = templtot_per_feed;
  };

  // print points
  if (0 < bright::verbosity)
  {
    std::cout << "last_ Point: M* = " << prev_Mstar << "\tL/F = " << prev_ltot_per_feed << "\n";
    std::cout << "curr_ Point: M* = " << curr__Mstar << "\tL/F = " << curr__ltot_per_feed << "\n";
  };

  // Start iterations.    
  while (xpn < ooe)
  {
    // Check that parameters are still well-formed
    if ( isnan(curr__Mstar) || isnan(curr__ltot_per_feed) || isnan(prev_Mstar) || isnan(prev_ltot_per_feed) )
      throw EnrichmentIterationNaN();

    prev_Mstar  = curr__Mstar;
    prev_ltot_per_feed = curr__ltot_per_feed;

    curr__Mstar = curr__Mstar - (m_sign * pow(10.0, -xpn));
    Mstar = curr__Mstar;
    ltot_per_feed();
    curr__ltot_per_feed = TotalPerFeed;

    if (prev_ltot_per_feed < curr__ltot_per_feed)
    {
      tempMstar = curr__Mstar - (m_sign * pow(10.0, -xpn));
      Mstar = tempMstar;
      ltot_per_feed();
      templtot_per_feed = TotalPerFeed; 

      tempm = pyne::slope(curr__Mstar, curr__ltot_per_feed, tempMstar, templtot_per_feed);
      if (tempm == 0.0)
      {
        prev_Mstar  = curr__Mstar;
        prev_ltot_per_feed = curr__ltot_per_feed;
        curr__Mstar  = tempMstar;
        curr__ltot_per_feed = templtot_per_feed;

        // print Point
        if (0 < bright::verbosity)
          std::cout << "Next Point: M* = " << curr__Mstar << "\tL/F = " << curr__ltot_per_feed << "\n";
        break;
      };

      tempm_sign = tempm / fabs(tempm);
      if (m_sign != tempm_sign)
      {
        xpn = xpn + 1;

        tempMstar = prev_Mstar + (m_sign * pow(10.0, -xpn));
        Mstar = tempMstar;
        ltot_per_feed();
        templtot_per_feed = TotalPerFeed;
        tempm = pyne::slope(prev_Mstar, prev_ltot_per_feed, tempMstar, templtot_per_feed);

        if (tempm == 0.0)
        {
          curr__Mstar  = tempMstar;
          curr__ltot_per_feed = templtot_per_feed;

          // print Point
          if (0 < bright::verbosity)
            std::cout << "Next Point: M* = " << curr__Mstar << "\tL/F = " << curr__ltot_per_feed << "\n";
          break;
        };

        m_sign = tempm / fabs(tempm);
      };
    };
  };

  Mstar        = curr__Mstar;
  TotalPerFeed = curr__ltot_per_feed;
  return;
};

*/
