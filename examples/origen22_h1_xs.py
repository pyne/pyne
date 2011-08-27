"""This example irradiates 1 kg of water for 100 days in ORIGEN 2.2, increasing the 
capture cross section of Hydrogen-1 by 10% each time.  The Hydrogen-2 concentration 
is then gathered and displayed.
"""
from subprocess import check_call

from pyne import origen22
from pyne.api import Material


# 1 kg of water
water = Material()
water.from_atom_frac({'H1': 2.0, 'O16': 1.0})
water.mass = 1E3

# Make a tape4 file for water
origen22.write_tape4(water)

# Make a tape 5 for this calculation
#   * Just output the concentration tables
#   * The cross-section library numbers must 
#     the library / deck numbers in tape9 
origen22.write_tape5_irradiation("IRF", 1000.0, 4E14, 
                                 xsfpy_nlb=(381, 382, 383), 
                                 out_table_num=[5])

# Grab a base tape9 from which we will overlay new values
# This must be supplied by the user
base_tape9 = origen22.parse_tape9("BASE_TAPE9.INP")

base_h1_xs = base_tape9[381]['sigma_gamma'][10010]

# Init a dumb overlay tape9
overlay_tape9 = {381: {'_type': 'xsfpy', 
                       '_subtype': 'activation_products', 
                       'sigma_gamma': {10010: base_h1_xs}, 
                       }
                }


# Run origen, increasing the cross section each time.
h2_concentration = []
for i in range(11):
    overlay_tape9[381]['sigma_gamma'][10010] = (1.0 + i*0.1) * base_h1_xs

    # Merge the base and overlay, and write out
    new_tape9 = origen22.merge_tape9([overlay_tape9, base_tape9])
    origen22.write_tape9(new_tape9, 'TAPE9.INP')

    # Run and parse origen output
    rtn = check_call(['o2_therm_linux.exe'])
    tape6 = origen22.parse_tape6('TAPE6.OUT')
    h2_concentration.append(tape6['table_5']['summary']['data'][10020][0][-1])

print 
print "H2 Concentration: ", h2_concentration


# Clean up
import os
for f in os.listdir('.'):
    if (f.endswith('.INP') or f.endswith('.OUT')) and f != 'BASE_TAPE9.INP':
        os.remove(f)

# Try to plot the results
try: 
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit

import numpy as np
    
plt.plot((1.0 + np.arange(11)*0.1)*base_h1_xs, h2_concentration, 'k-', label='H2')
plt.xlabel('H1 capture cross section [barns]')
plt.ylabel('H2 concentration [grams]')
plt.legend(loc=0)
plt.show()


