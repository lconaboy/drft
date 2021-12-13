"""
Snap coordinates to the base grid, so that you can guarantee
equal extents.

Tests
=====

Rounding
--------

ce = 0.00390625 + 0.001953126 # centre

print(nco_s)
if (nco - nco_s > 0.5): nco_s += 1
print(nco_s)

>>> 0
>>> 1

ce = 0.00390625 + 0.001953126 # centre

print(nco_s)
if (nco - nco_s > 0.5): nco_s += 1
print(nco_s)

>>> 0
>>> 1

ce = 0.00390625 + 0.001953125 # centre

>>> 0
>>> 1

ce = 0.00390625 + 0.001953124 # centre

>>> 0
>>> 0

"""
# ce, co, ex all in [0..1]
ce = (0.08984, 0.82031, 0.5000) # centre
co = None # corner
ex = 0.0078125  # extent
lmin = 8
nmin = 2 ** lmin
dmin = 1. / nmin

if ce is not None:
    co = (x - (ex / 2.) for x in ce)

nco = (x * nmin for x in co)
nco_s = (int(x) if (x - int(x)) < 0.5 else int(x) + 1 for x in nco)
co_s = (x * dmin for x in nco_s)

print('{0:12.11f}, {1:12.11f}, {2:12.11f}'.format(*co_s))
