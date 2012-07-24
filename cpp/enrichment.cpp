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


double pyne_enr::alphastar_i(double alpha_0, double Mstar, double M_i)
{
  // M_i is the mass of the ith nuclide
  return pow(alpha_0, (Mstar - M_i));
}

double pyne_enr::Ei(double astar_i, double N)
{
  return ((astar_i - 1.0) / (1.0 - pow(astar_i, -N)));
};


double pyne_enr::Si(double astar_i, double M)
{
  return ((astar_i - 1.0)/(pow(astar_i, M+1) - 1.0));
};

/*

void pyne_enr::FindNM()
{
  // This give the order-of-exactness to which N and M are solved for.
  double ooe = 7.0;
  double tolerance = pow(10.0, -ooe);

  double PoF = prod_per_feed(mat_feed.comp[j], xP_j, xW_j);
  double WoF = tail_per_feed(mat_feed.comp[j], xP_j, xW_j);
  double alphastar_j = alphastar_i(pyne::atomic_mass(j));

  // Save original state of N & M
  double origN = N;
  double origM = M;

  double lhsP = PoF * xP_j / mat_feed.comp[j];
  double rhsP = (pow(alphastar_j, M+1.0) - 1.0) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
  double lhsW = WoF * xW_j / mat_feed.comp[j];
  double rhsW = (1.0 - pow(alphastar_j, -N)) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));

  double n = 1.0;
  while (tolerance < fabs(lhsP - rhsP) && tolerance < fabs(lhsW - rhsW))
  {
    if (tolerance < fabs(lhsP - rhsP))
    {
      N = N - ((lhsP - rhsP) * N);
      rhsP = (pow(alphastar_j, M+1.0) - 1.0) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
    };

    if (tolerance < fabs(lhsW - rhsW))
    {
      M = M - ((lhsW - rhsW) * M);
      rhsW = (1.0 - pow(alphastar_j, -N)) / (pow(alphastar_j, M+1.0) - pow(alphastar_j, -N));
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
  return; 
};
  

double pyne_enr::xP_i(int i)
{
  double alphastar_i = alphastar_i(pyne::atomic_mass(i));
  double numerator = mat_feed.comp[i]*(pow(alphastar_i, M+1.0) - 1.0);
  double denominator = (pow(alphastar_i, M+1.0) - pow(alphastar_i, -N)) / prod_per_feed(mat_feed.comp[j], xP_j, xW_j);
  return numerator / denominator;
};


double pyne_enr::xW_i(int i)
{
  double alphastar_i = alphastar_i(pyne::atomic_mass(i));
  double numerator = mat_feed.comp[i] * (1.0 - pow(alphastar_i, -N));
	double denominator = (pow(alphastar_i, M+1.0) - pow(alphastar_i, -N)) / tail_per_feed(mat_feed.comp[j], xP_j, xW_j);
  return numerator / denominator;
};


void pyne_enr::SolveNM()
{
  //This function takes a given initial guess number of enriching and stripping stages 
  //for a given composition of fuel with a given jth key component, knowing the values 
  //that are desired in both Product and Tails streams.  Having this it solves for what 
  //the actual N and M stage numbers are and also what the product and waste streams 
  //compositions are.  It returns precisely these.

  FindNM();

  pyne::comp_map compP;
  pyne::comp_map compW;

  for (pyne::comp_iter i = mat_feed.comp.begin(); i != mat_feed.comp.end(); i++)
  {
    compP[i->first] = xP_i(i->first);
    compW[i->first] = xW_i(i->first);
  };

  mat_prod  = pyne::Material(compP);
  mat_tail = pyne::Material(compW);

  return;
};


void pyne_enr::_norm_comp_secant(double N0, double tolerance)
{
  // This function actually solves the whole system of equations.  It uses SolveNM 
  // to find the roots for the enriching and stripping stage numbers.  It then 
  // checks to see if the product and waste streams meet their target enrichments
  // for the jth component like they should.  If they don't then it trys other values 
  // of N and M varied by Newton's Method.  Rinse and repeat as needed.

  // Is the history of N and M that has been input
  std::vector<double> historyN;
  std::vector<double> historyM;

  // Start iteration Counter
  int counter = 0;

  // Set first two points
  double last_N = N0 + 1.0;
  double last_M = M0 + 1.0;
  double curr_N = N0;
  double curr_M = M0;

  // Initialize 'last_' point
  N = last_N;
  M = last_M;
  SolveNM();
  double last_xP_j = mat_prod.comp[j];
  double last_xW_j = mat_tail.comp[j];
  historyN.push_back(N);
  historyN.push_back(M);

  // Initialize 'curr_ent' point
  N = curr_N;
  M = curr_M;
  SolveNM();
  double curr_xP_j = mat_prod.comp[j];
  double curr_xW_j = mat_tail.comp[j];
  historyN.push_back(N);
  historyN.push_back(M);

  // My guess is that what we are checkin here is that the isotopic compositions
  // make sense with abs(1.0 - masscurr_P) rather than calculatign the 
  // relative product to watse mass streams.
  double tempcurr_N = 0.0;
  double tempcurr_M = 0.0;
  double templast_N = 0.0;
  double templast_M = 0.0;

  while (tolerance < fabs(xP_j - curr_xP_j) || tolerance < fabs(xW_j - curr_xW_j))
  {
    if (1 < bright::verbosity)
      std::cout << "--------------------\n";

    if (tolerance <= fabs(xP_j - curr_xP_j))
    {
      // Make a new guess for N
      tempcurr_N = curr_N;
      templast_N = last_N;
      curr_N = curr_N + (xP_j - curr_xP_j)*((curr_N - last_N)/(curr_xP_j - last_xP_j));
      last_N = tempcurr_N;

      // If the new value of N is less than zero, reset.
      if (curr_N < 0.0)
      {
        curr_N = (tempcurr_N + templast_N)/2.0;
        if (1 < bright::verbosity)
          std::cout << "    N < 0, resetting.\n";
      };
    };

    if (tolerance <= fabs(xW_j - curr_xW_j))
    {
      // Make a new guess for M
      tempcurr_M = curr_M;
      templast_M = last_M;
      curr_M = curr_M + (xW_j - curr_xW_j)*((curr_M - last_M)/(curr_xW_j - last_xW_j));
      last_M = tempcurr_M;

      // If the new value of M is less than zero, reset.
      if (M < 0.0)
      {
        curr_M = (tempcurr_M + templast_M)/2.0;
        if (1 < bright::verbosity)
          std::cout << "    M < 0, resetting.\n";
      };
    };

    // Check for infinite loops
    for (int h = 0; h < historyN.size(); h++)
    {
      if (historyN[h] == curr_N && historyM[h] == curr_M)
      {
        if (1 < bright::verbosity)
          std::cout << "~~~ Infinite loop found and exception thrown! ~~~.\n";
        throw EnrichmentInfiniteLoopError();
      };
    };

    if (150 <= historyN.size())
    {
      historyN.erase(historyN.begin());
      historyM.erase(historyM.begin());
    };
    historyN.push_back(N);
    historyM.push_back(M);

    if (10000 < counter)
    {
      if (1 < bright::verbosity)
        std::cout << "~~~ Secant method counter limit hit! ~~~.\n";
      throw EnrichmentIterationLimit();
    }
    else
    {
      counter = counter + 1;
    };

    // Calculate new isotopics for valid (N, M)        
    last_xP_j = curr_xP_j;
    last_xW_j = curr_xW_j;

    N = curr_N;
    M = curr_M;
    SolveNM();
    curr_xP_j = mat_prod.comp[j];
    curr_xW_j = mat_tail.comp[j];

    if (1 < bright::verbosity)
    {
      std::cout << "Product Mass: " << curr_xP_j << "\tTails Mass: " << curr_xW_j << "\n";
      std::cout << "====================\n";
    };
  };

  return;
};


// I have serious doubts that this works...
void pyne_enr::_norm_comp_other()
{
  // This give the order-of-exactness to which N and M are solved for.
  double ooe = 5.0;
  double tolerance = pow(10.0, -ooe);

  // Is the hisorty of N and M that has been input
  std::vector<double> historyN;
  std::vector<double> historyM;

  // Initial point
  N = N0;
  M = M0;
  SolveNM();
  double massP = mat_prod.mass;
  double massW = mat_tail.mass;

  while (tolerance < fabs(1.0 - massP) && tolerance < fabs(1.0 - massW) )
  {
    if (tolerance <= fabs(1.0 - massP))
      N = N - (1.0 - massP)/(1.0 + massP);

    if (tolerance <= fabs(1.0 - massW))
      M = M + (1.0 - massW)/(1.0 + massW);

    // Note this infinite loop checker does not raise an exception
    // Thus the exception cannot be caught by a try statement and then another
    // root finding Comp2Unity method tried.
    // This simply leaves N and M in whatever their curr_ent state is at the time.
    // This is fine for M* = 235.1 for U, since M* won't be optimized here as an outlier 
    // and since it is looping it is probably signalling around some actual value.

    // Check for infinite loops
    for (int h = 0; h < historyN.size(); h++)
    {
      if (historyN[h] == N && historyM[h] == M)
      {
        //throw EnrichmentInfiniteLoopError(); //Possible future use.
        return;
      };

      if (150 <= historyN.size())
      {
        historyN.erase(historyN.begin());
        historyM.erase(historyM.begin());
      };
    };
    historyN.push_back(N);
    historyM.push_back(M);

    // Calculate new masses
    SolveNM();
    massP = mat_prod.mass;
    massW = mat_tail.mass;
  };

  return;
};


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
    double PoF = prod_per_feed(mat_feed.comp[j], xP_j, xW_j);
    double WoF = tail_per_feed(mat_feed.comp[j], xP_j, xW_j);

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
