import numpy as np

def cic(x, dx, verbose=False):
    """Takes a particle based quantity and interpolates it onto a grid
    using cloud-in-cell interpolation, so that it can be used in
    conjunction with gridded data. Designed to work with three
    dimensional data.

    :param x: 
        particle quantity to be gridded, can be either (array) of same
        dimension as the gridded data or (scalar) for e.g. mass
    :param dx:
        particle offsets from the grid centres, should be normalised so that a 
        cell width is 1.0, should be of type (array) and the gridded data will
        have shape dx.shape
    :param verbose:
        (bool) want printout?
    :returns:
        the particle values interpolated on to a grid
    :rtype:
        (array)

    """
    assert verbose in [True, False, None], 'verbose flag should be bool or None'
    
    # Get dimensions of data
    n1, n2, n3 = dx.shape
    xg = np.zeros(shape=(n1, n2, n3))

    # Do interpolation
    if verbose: print('Interpolating...')
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                # Allow for periodic boundary conditions
                ip = i % n1
                jp = j % n2
                kp = k % n3

                # Define convenience variables
                dxx = dx[i, j, k]
                dyy = dx[i, j, k]
                dzz = dx[i, j, k]
                txx = 1.0 - dxx
                tyy = 1.0 - dyy
                tzz = 1.0 - dzz

                # Interpolate over eight neighbouring cells
                xg[i, j, k] += x[i, j, k] * txx * tyy * tzz
                xg[ip, j, k] += x[ip, j, k] * dxx * tyy * tzz
                xg[i, jp, k] += x[i, jp, k] * txx * dyy * tzz
                xg[i, j, kp] += x[i, j, kp] * txx * tyy * dzz
                xg[ip, jp, k] += x[ip, jp, k] * dxx * dyy * tzz
                xg[ip, j, kp] += x[ip, j, kp] * dxx * tyy * dzz
                xg[i, jp, kp] += x[i, jp, kp] * txx * dyy * dzz
                xg[ip, jp, kp] += x[ip, jp, kp] * dxx * dyy * dzz

    return xg
