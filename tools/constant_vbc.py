"""
Imposes a constant v_bc on the initial conditions files for a
specified level and magnitude of the v_bc.
"""

import sys
import grafic_tools as grafic

if len(sys.argv) < 3:
    print('Usage: python constant_vbc.py <level> <v_bc>')
    sys.exit()
    
level = int(sys.argv[1])
v = float(sys.argv[2])

grafic.impose_vbc('./', level, v)
