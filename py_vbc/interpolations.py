import numpy as np
from scipy.interpolate import interp1d, interp2d

from py_vbc.constants import *

def interpolate_tf(flag, z, dz=None):
    """Reads in CAMB transfer functions at z=0, z=zstart, z=zstart+dz and
    interpolates w.r.t. k.

    :param flag: (str) 'c' for dark matter, 'b' for baryons, 'g' for radiation,
                       't' for total
    :param z: (float) redshift to interpolate TFs at 
    :param dz: (float) delta(z) of TFs to calculate dTF/dz
    :returns: tfz_spline -- the spline of the TF at z=zstart, 
              dtf_spline -- the spline of the derivative of the TF 
                            around z=zstart
    :rtype: splines, which are used by supplying them with a a range of 
            x-values like y_interp = spline(x)

    """

    assert flag in ['c', 'b', 'g', 't', 'vb', 'vc'], 'Flag should be either "c", "vc", "b", "vb", "g" or "t".'
    
    # Dict to convert flag to a TF index
    i_flag = {'c':1, 'b':2, 'g':3, 't':6, 'vc':10, 'vb':11}
    
    z = int(z)

    if (dz is None) and (z==0): 
        # Determine filename
        fn0 = tf_base.format(camb_base, 0)

        # Load in TFs
        tf0 = np.loadtxt(fn0)

        # Convert from k/h to just k
        k = tf0[:, 0]
        kh = k*hconst

        # Compute the spline for z=0
        tf0_spline = interp1d(kh, tf0[:, i_flag[flag]], kind='cubic')

        return tf0_spline

    else:
        # Calculate the spline for z=zstart
        fnz = tf_base.format(camb_base, z)
        tfz = np.loadtxt(fnz)

        # Convert from k/h to k
        k = tfz[:, 0]
        kh = k*hconst
        tfz_spline = interp1d(kh, tfz[:, i_flag[flag]], kind='cubic')

        # If dz is supplied, calculate the derivative
        if dz is not None:
            fnzm = tf_base.format(camb_base, int(z-dz))
            fnzp = tf_base.format(camb_base, int(z+dz))

            # the first line of TF files be the number of lines in the TF
            # Load in TFs
            tfzm = np.loadtxt(fnzm)
            tfzp = np.loadtxt(fnzp)

            # Check TFs are same length
            if (tfz.shape[0] != tfzm.shape[0] != tfzp.shape[0]):
                raise Exception('Input transfer functions are different lengths.')

            # Calculate the derivative of the TFs
            dtf = (tfzm - tfzp)/(2.0 * dz)
    
            # Compute the splines for dz
            dtf_spline = interp1d(kh, dtf[:, i_flag[flag]], kind='cubic')

            return tfz_spline, dtf_spline

        else:
            return tfz_spline


def interpolate_tf2d(flag, zs):
    """Reads in CAMB transfer functions at z={zs} and interpolates wrt
    to k and z.  Returned spline is a function f(k, z)

    :param flag: (str) 'c' for dark matter, 'b' for baryons, 'g' for
        radiation, 't' for total
    :param zs: (array) redshifts to interpolate TFs at

    :returns: tfz_spline -- the spline of the TF for k and z

    :rtype: splines, which are used by supplying them with a a range
            of x-values like y_interp = spline(x)
    """

    assert flag in ['c', 'b', 'g', 't', 'vb', 'vc'], 'Flag should be either "c", "vc", "b", "vb", "g" or "t".'
    
    # Dict to convert flag to a TF index
    i_flag = {'c':1, 'b':2, 'g':3, 't':6, 'vc':10, 'vb':11}

    fnz = tf_base.format(camb_base, zs[0])
    tfz = np.loadtxt(fnz)
    k = tfz[:, 0]
    kh = k * hconst  # convert from h/Mpc to 1/Mpc
    # Set up array to interpolate over
    nz = len(zs)
    nk = tfz.shape[0]
    tfkz = np.zeros((nz, nk))

    for i, z in enumerate(zs):
        fnz = tf_base.format(camb_base, z)
        tfz = np.loadtxt(fnz)
        tfkz[i, :] = tfz[:, i_flag[flag]]        
    
    # Convert from k/h to k
    # k = tfz[:, 0]
    # kh = k*hconst
    # tfz_spline = interp1d(kh, tfz[:, i_flag[flag]], kind='cubic')

    tfkz_spline = interp2d(kh, zs, tfkz, kind='cubic')
    
    return tfkz_spline


def interpolate_recfast():
    """Interpolates RECFAST data w.r.t. redshift for later use.

    :param test: (bool) use test RECFAST file
    :returns: spline of temperature and electron fraction values
    :rtype: spline

    """

    vals = np.loadtxt(rf_base)
    # Reverse order so redshift is increasing
    vals = vals[::-1]

    # Extract the variables
    z = vals[:, 0]
    xe = vals[:, 1]
    T = vals[:, -1]  # RECFAST++ has extra columns

    # Compute the splines, again why are we doing this? If the spline
    # is the same length as the input array isn't this pointless?
    T_spline = interp1d(z, T, kind='cubic')
    xe_spline = interp1d(z, xe, kind='cubic')
    
    return T_spline, xe_spline
