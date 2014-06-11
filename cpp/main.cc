#include <iostream>
#include <string>

#include "pyne.h"
#include "data.h"
#include "nucname.h"
#include "rxname.h"

#define SHOW(X) \
  std::cout << __FILE__ << ":" << __LINE__ << ": "#X" = " << X << "\n"

int main(int argc, char* argv[]) {
  pyne::NUC_DATA_PATH = "/home/opotowsky/.local/lib/python2.7/site-packages/pyne/nuc_data.h5";

  double air = pyne::ext_air_dose("Cs137", 0);
//  double soil = pyne::ext_soil_dose("Cs137", 0)
//  double ingest = pyne::ingest_dose("Cs137", 0);
//  double inhale = pyne::inhale_dose("Cs137", 0);
  SHOW(air);
//  SHOW(soil);
//  SHOW(ingest);
//  SHOW(inhale);

  return 0;
}
