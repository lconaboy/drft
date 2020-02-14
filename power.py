import numpy as np
from cosmology import compute_power

def main(fn, boxsize=1.0):
    x = np.loadtxt(fn)
    power, k  = compute_power(x, boxsize)

    np.savetxt(fn+'.power', power)
    np.savetxt(fn+'.bins', k)

    return 0

if __name__=='__main__':
    import sys

    if len(sys.argv) < 2:
        print('Usage: python power.py filename [boxsize]')
        sys.exit()
    if len(sys.argv) < 3:
        print('Defaulting to boxsize = 1.0')
        fn = sys.argv[1]
        boxsize = 1.0
    else:
        fn = sys.argv[1]
        boxsize = sys.argv[2]

    main(fn, boxsize)
    
