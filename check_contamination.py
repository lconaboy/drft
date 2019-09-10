import seren3
from pymses.utils import constants as C
import numpy as np
import matplotlib.pyplot as plt

# TODO 
#  - change to cube
#  - produce contamination maps

def check_contam(i_snap, cen=(0.5, 0.5, 0.5), rad=0.025, snap_dir='.'):
    """Prouduces a histogram of the relative masses of particles as a
    function of their distance from the centre of a spherical filter.

    :param i_snap: (int) which output to load
    :param cen: (tuple) centre of the spherical filter, in code units
    :param rad: (float) radius of the spherical filter, in code units
    :param snap_dir: (str) directory of the snapshots
    :returns: None
    :rtype:

    """
    # Load snapshot
    print('Working on snapshot {0}'.format(i_snap))
    s = seren3.load_snapshot(snap_dir, i_snap)

    # Calculate comoving box length
    a = 1.0/(1.0 + s.z)
    h = s.info["H0"]/100.0
    L = s.info["unit_length"].express(C.Mpc)
    cL = L/a * h

    # Centre and radius of spherical filter
    cen = cen
    rad = rad  # convert to code units

    # Apply spherical filter and extract dark matter position and mass
    sub = s[s.get_sphere(cen, rad)]
    sub = sub.d
    val = sub[["pos", "mass"]].flatten()

    # Calculate distance from the centre for each particle
    s_cen = s.array(cen, s.info["unit_length"])
    s_rad = s.array(rad, s.info["unit_length"])
    
    # print('s_rad = {0}'.format(s_rad))

    pos = val["pos"] - s_cen
    pos = np.linalg.norm(pos, axis=1)
    pos_norm = pos/s_rad
    pos_arr = np.asarray(pos_norm.tolist())/float(s_rad.units)

    # Calculate mass of each particle
    mass = val["mass"]
    # Normalise the mass with respect to the smallest (finest)
    mass_norm = mass/np.min(mass)
    mass_norm_ref = np.log2(mass_norm)/3.0
    mass_arr = np.asarray(mass_norm_ref.tolist())

    # Calculate comoving spherical radius from comoving box length
    crad_phys = rad * cL
    
    ### PLOTTING ###

    # # Scatter
    # fig, ax = plt.subplots(figsize=(6, 6))
    # ax.set_xlim([0, 1])
    # ax.set_ylim([0, np.max(mass_arr)])
    # for i in range(int(np.max(mass_arr)+1)):
    #     ax.scatter(pos_arr[mass_arr==i], mass_arr[mass_arr==i])
    # ax.set_ylabel("log$_2$[m$_i$/min(m)]/3")
    # ax.set_xlabel("x/r")
    # ax.legend()
    # fig.tight_layout()
    # fig.savefig("scatter.png")

    # Hist
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim([0, 1])

    bins = np.linspace(0, 1, 100)
    for i in range(int(np.max(mass_arr))+1):
        ax.hist(pos_arr[mass_arr==i], bins=bins, histtype="step", label=r"$\ell - \ell_\mathrm{{max}}$ = {}".format(i))
        # Keep the axis limits fixed
        if i == 0:
            # The initial snapshot will have the most particles
            lims = ax.get_ylim()

        # Lower limit of 1 to work with log scale
        ax.set_ylim([1, 1.1*lims[1]])

        ax.set_xlabel(r"r/r$_0$ (r$_0$ = L/2)".format(crad_phys))
    ax.set_ylabel("N")
    ax.set_yscale("log")
    ax.set_title('{0:3.2f} cMpc/h zoom region, z = {1:3.2f}'.format(crad_phys * 2.0, s.z))
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig("contam_hist_o{0:03d}.png".format(i_snap))


def check_cube(i_snap, cen=(0.5, 0.5, 0.5), l=0.025, snap_dir='.'):
    # Load snapshot
    print('Working on snapshot {0}'.format(i_snap))
    s = seren3.load_snapshot(snap_dir, i_snap)

    # Filter to a cube
    print('Filtering')
    sub = s[s.get_sphere(cen, l/2)]
    # Filter to DM
    sub = sub.d
    # Extract mass and position values
    print('Flattening')
    vals = sub[["pos", "mass"]].flatten()
    pos = vals["pos"]
    mass = vals["mass"]

    # Normalise the position to be between [0, 1]
    print('Normalising')
    x_min = cen[0] - l/2.0
    x_max = cen[0] + l/2.0
    pos = (pos - x_min)/x_max
    
    # Normalise the mass to the minimum
    mass /= np.min(mass)
    # Convert to refinement level
    mass = np.log2(mass)/3.0

    # # Pick out the level wanted
    # pos = pos[mass==level]
    # mass = mass[mass==level]
    
    # Convert both to arrays and return
    pos = np.array(pos.tolist(), dtype=float)
    mass = np.array(mass.tolist(), dtype=float)

    return pos, mass


def cube_to_grid(pos, mass, level, n=256):

    # # Pick out the level wanted
    pos = pos[mass==level]
    mass = mass[mass==level]

    # Digitize
    bins = np.linspace(0.0, 1.0, num=n)
    inds = np.digitize(pos, bins, right=True)

    vals = np.zeros(shape=(n, n, n))

    for i in range(pos.size):
        vals[inds[i]-1] = mass[i]

    # Axis can change dependent on which viewpoint is wanted
    vals = np.sum(vals, axis=2)

    return vals


if __name__ == '__main__':
    # Python 2.7
    import sys
    from mpi4py import MPI
    
    
    # Read in command line arguments
    o_min = int(sys.argv[1])
    o_max = int(sys.argv[2])
    o_step = int(sys.argv[3])
    rad = float(sys.argv[4])

    # Create a range of outputs to read from the args
    o_range = np.arange(o_min, o_max+1, o_step, dtype=int)
    pos, mass = check_cube(o_min, l=rad)

    vals = np.array([pos, mass])
    
    np.savetxt('vals.dat', vals)
    # # Set up MPI
    # comm = MPI.COMM_WORLD
    # size = comm.Get_size()
    # rank = comm.Get_rank()

    # if rank == 0:
    #     chunks = np.array_split(o_range, size)
    # else:
    #     o_range = None
    #     chunks = None

    # # Scatter the chunks
    # tasks = comm.scatter(chunks, root=0)

    # # Now produce the hist
    # for task in tasks:
    #     f = check_contam(task, rad=rad)

    # ff = comm.gather(f, root=0)

    # print('Done!')

    vals = np.loadtxt('vals.dat')
