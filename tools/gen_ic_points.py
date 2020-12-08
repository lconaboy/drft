"""Generates region point files for use by MUSIC.

Usage: python gen_ic_points.py <mode> [<extra args>]
"""

import sys
import numpy as np

def sphere(r, c):
    """Generates a set of points that MUSIC can use to compute a minimum
    bounding sphere. Works in box coordinates [0., 1.].

    r (float) - radius
    c (float) - centre
    
    """
    # More general, but more expensive
    # x = np.linspace(c[0] - 2.*r, c[0] + 2.*r, num=64)
    # y = np.linspace(c[1] - 2.*r, c[1] + 2.*r, num=64)
    # z = np.linspace(c[2] - 2.*r, c[2] + 2.*r, num=64)

    # xg, yg, zg = np.meshgrid(x, y, z)
    # rg = (xg-c[0])**2 + (yg-c[1])**2 + (zg-c[2])**2 
    
    # idxs = np.where(rg < r**2)

    # points = np.array([xg[idxs], yg[idxs], zg[idxs]]).T

    # Easy
    points = np.array([[c[0] - r, c[1], c[2]],
                       [c[0] + r, c[1], c[2]],
                       [c[0], c[1] - r, c[2]],
                       [c[0], c[1] + r, c[2]],
                       [c[0], c[1], c[2] - r],
                       [c[0], c[1], c[2] + r]])

    fn = 'points.dat'
    np.savetxt(fn, points, fmt='%.5f')

    print('---- written output to', fn)


def main(mode):

    mode_dict = {'s':'spherical'}

    print('---- generating a {0} region'.format(mode_dict[mode]))

    if mode == 's':
        if len(sys.argv) < 6:
            print('Usage: python gen_ic_points.py s <rad> <xc> <yc> <zc>')
            sys.exit(1)

        r = float(sys.argv[2])
        c = (float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))

        print('-------- using rad = {0:.5f} cen = {1:.5f}, {2:.5f}, {3:.5f}'.format(r, c[0], c[1], c[2]))
        
        sphere(r, c)

if __name__ == '__main__':

    mode = sys.argv[1]

    main(mode)

