import numpy as np
import grafic_tools as grafic
from cosmology import power_spectrum


def main(fn, boxsize=1.0, quickplot=True):
    cube = grafic.load_cube(fn)
    power, k  = power_spectrum(cube.load_box(), boxsize)
    
    np.savetxt(fn+'.power', power)
    np.savetxt(fn+'.bins', k)

    print(np.max(power))

    if quickplot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 6))
        plt.plot(k[1:], power[1:], c='k', label=fn)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'k (h/Mpc)')
        plt.ylabel(r'P(k) (h/Mpc)${{^-3}}$')
        plt.legend()
        plt.tight_layout()
        plt.savefig(fn+'_quickplot.png')

    return 0

if __name__=='__main__':
    import sys

    if len(sys.argv) < 2:
        print('Usage: python power.py filename [boxsize]')
        sys.exit()
    if len(sys.argv) < 3:
        fn = sys.argv[1]
        boxsize = 1.0
    else:
        fn = sys.argv[1]
        boxsize = float(sys.argv[2])

    print('Using boxsize = {0:5.3f} h^-1 Mpc'.format(boxsize))
    main(fn, boxsize)
    
