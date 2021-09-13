"""Based on a script originally written by David Sullivan, which can
found at https://github.com/sully90/seren3/"""

import os
import gc
import sys
import time
import glob
import pickle
import numpy as np
from mpi4py import MPI  

import utils as vbc_utils
import grafic_tools as grafic

class Result(object):
    '''
    Simple wrapper object to contain result of single iteration MPI computation
    '''
    def __init__(self, rank, idx):
        self.rank = rank
        self.idx = idx
        self.result = None

    def __repr__(self):
        return "rank: %d idx: %s result: %s" % (self.rank, self.idx, self.result)

    def __eq__(self, other):
        return self.result == other.result

    def __ne__(self, other):
        return self.result != other.result

    def __hash__(self):
        return hash(self.result)

    
class Patch(object):

    def __init__(self, patch, dx, data, field):
        self.patch = patch
        self.dx = dx
        self.data = data
        self.field = field


def work(path, level, patch_size, levelmin, lin=False, verbose=True, ret_vbc=False):
    """Computes a new set of biased grafic fields by convolving patches of
    the fields with bias factors computed using py_vbc.

    :param path: (str) path to (unbiased) grafic fields
    :param level: (int) level of ICs to work on
    :param patch_size: (float) patch size in comoving Mpc
    :param levelmin: (int) overall minimum level of ICs, used for determining whether
    to use periodic CIC interpolation
    :param lin: (bool) only bias deltac/deltab (True) or deltab/velb* (False)?
    :param verbose: (bool) controls printout
    :param ret_vbc: (bool) testing feature, generate an IC file of the patched vbc
    :returns: 
    :rtype:
    """

    # MPI stuff
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    barrier = comm.Barrier
    finalize = MPI.Finalize

    if (lin and (rank == 0)):
        vbc_utils.msg(rank, 'Biasing both deltab and deltac fields', verbose)
    elif (rank == 0):
        vbc_utils.msg(rank, 'Biasing deltab and velb fields', verbose)

    if (ret_vbc and (rank == 0)):
        vbc_utils.msg(rank, 'TESTING: writing vbc_patch field', verbose)
        
    #mpi.msg("Loading initial conditions")
    vbc_utils.msg(rank, "Loading initial conditions.", verbose)
    
    if rank == 0:        
        # First, make root patches dir if it doesn't exist
        if not os.path.isdir("./patches"):
            os.mkdir("./patches")
        # Make patches dir for level
        if os.path.isdir("./patches/level_{0:03d}".format(level)):
            raise Exception("'patches' directory already exists. Either run in 'write' mode or remove and re-run in 'work' mode.")
        elif not os.path.isdir("./patches/level_{0:03d}".format(level)):
            os.mkdir("./patches/level_{0:03d}".format(level))
            vbc_utils.msg(rank, 'Made patches directory.', verbose)

        # Next, derive any fields we need that don't already
        # exist. Start with vbc. Need to know if we are doing periodic
        # or non-periodic.
        # per = (level == levelmin)
        
        if not os.path.isfile(path+"level_{0:03d}/ic_vbc".format(level)):
            # vbc_utils.msg(rank, 'Deriving ic_vbc.', verbose)
            # grafic.derive_vbc(path, level, per)
            raise Exception("'ic_vbc' doesn't exist. Run cic_vel.f90.")
        # Make sure the deltac field exists, if we need it
        if (lin) and (not os.path.isfile(path+"level_{0:03d}/ic_deltac".format(level))):
            # vbc_utils.msg(rank, 'Deriving ic_deltac.', verbose)
            # grafic.derive_deltac(path, level, per)
            raise Exception("'lin' mode is not currently supported.")
    else:
        # Wait for rank 0 to write the velb, velc and vbc fields
        while not os.path.isfile(path+"level_{0:03d}/ic_vbc".format(level)):
            time.sleep(1.0e-3)

    if lin:
        ics = [grafic.load_snapshot(path, level, field=field) for field in
               ['deltab', 'deltac', 'vbc']]
    else:
        ics = [grafic.load_snapshot(path, level, field=field) for field in
               ['deltab', 'velbx', 'velby', 'velbz', 'vbc']]

    # This only works with cubic zoom regions, so check that and exit
    # if not
    if ((ics[0].n[0] != ics[0].n[1]) or (ics[0].n[0] != ics[0].n[2])):
        raise Exception('Only works with cubic zoom regions. Use force_equal_extent = yes with music to generate cubic zoom regions.')

    # Padding around the patches
    pad = 8
    # print('---- testing with extra padding')
    # pad = 16

    # Calculate the number of cubes for each dimension, we want an
    # equal number of cubes in each dimension and this is most easily
    # achieved if the zoom region is actually cubic
    ncubes = np.zeros(3)
    for i in range(3):
        # div = np.array([float(i) for i in vbc_utils.divisors(ics[0].n[i], mode='yield')])
        div = np.array([float(i) for i in vbc_utils.divisors(ics[0].n[i], mode='yield')
                        if i <= 20])  # rough limit to stop using too small patches
        idx = np.abs((ics[0].n[i] / div) * ics[0].dx - patch_size).argmin()
        ncubes[i] = int(div[idx])

    # Compute cube positions in cell units
    cubes = dx = None
    if rank == 0:
        vbc_utils.msg(rank, "Using {0} cubes per dimension.".format(ncubes), verbose)

        cubes, dx = vbc_utils.cube_positions(ics[0], ncubes, ics[0].N)
        cubes = np.array(cubes)
        # Split the cubes into chunks that can be scattered to each processor 
        chunks = np.array_split(cubes, size)
    else:
        chunks = None
    
    # dx is the same for every cube, so this is broadcast to each processor
    dx = comm.bcast(dx, root=0)

    # Initialise \delta sums for correction
    deltab_tot = 0
    deltab_b_tot = 0

    # Iterate over patch positions in parallel
    # dest = {}
    patches = comm.scatter(chunks, root=0)

    # Now the work branches into two, depending on whether we're using
    # the linear method or not
    if lin:
        # For the linear method, we only bias the delta fields
        for i, patch in enumerate(patches):
            # Always print this
            vbc_utils.msg(rank, '{0}/{1}'.format(i+1, len(patches)))

            # Convert dx to float
            dx = dx.astype(np.float32)

            origin = np.array(patch - dx / 2. - pad, dtype=np.int64)
            dx_eps = dx + float(2 * pad)

            # Convert dx_eps to int
            dx_eps = dx_eps.astype(np.int32)

            # Initialise
            deltab = None
            deltac = None
            vbc = None

            vbc_utils.msg(rank, "Loading patch: {0}.".format(patch), verbose)
            # vbc_utils.msg(rank, "TESTING: origin {0}, dx_eps {1}".format(origin, dx_eps), verbose)
            deltab = ics[0].load_patch(origin, dx_eps)
            deltac = ics[1].load_patch(origin, dx_eps)
            vbc = ics[2].load_patch(origin, dx_eps)

            # Compute the bias
            vbc_utils.msg(rank, "Computing bias.", verbose)
            # Commented the below for testing
            k, b_c, b_b, b_vc, b_vb = vbc_utils.compute_bias(ics[2], vbc)

            # Convolve with field
            vbc_utils.msg(rank, "Performing convolution.", verbose)
            deltab_biased = vbc_utils.apply_density_bias(ics[0], k, b_b, deltab.shape[0], delta_x=deltab)
            deltac_biased = vbc_utils.apply_density_bias(ics[1], k, b_c, deltac.shape[0], delta_x=deltac)

            # print('deltab before', delta_biased)
            # Remove the padded region
            x_shape, y_shape, z_shape = deltab_biased.shape
            deltab_biased = deltab_biased[0 + pad:x_shape - pad,
                                          0 + pad:y_shape - pad,
                                          0 + pad:z_shape - pad]

            x_shape, y_shape, z_shape = deltac_biased.shape
            deltac_biased = deltac_biased[0 + pad:x_shape - pad,
                                          0 + pad:y_shape - pad,
                                          0 + pad:z_shape - pad]

            # Store the patched v_bc that is used to calculate the bias
            if ret_vbc:
                vbc_patch = np.ones_like(deltab_biased) * vbc_utils.vbc_rms(vbc)

                # Store
                biased_patches = [Patch(patch, dx, deltab_biased, 'deltab'),
                                  Patch(patch, dx, deltac_biased, 'deltac'),
                                  Patch(patch, dx, vbc_patch, 'vbc_patch')]
            else:
                # Store
                biased_patches = [Patch(patch, dx, deltab_biased, 'deltab'),
                                  Patch(patch, dx, deltac_biased, 'deltac')]

            for patch in biased_patches:
                with open(r"patches/level_{0:03d}/patch_{1}{2:03d}.p".format(level, patch.field, rank), "ab") as f:
                    pickle.dump(patch, f)

            del vbc
            del deltab
            del deltac
            del deltac_biased
            del deltab_biased
            if ret_vbc: del vbc_patch

            gc.collect()

    else:
        # For the other method, we only bias the baryon field
        # (including the velocity fields)
        for i, patch in enumerate(patches):
            # Always print this
            vbc_utils.msg(rank, '{0}/{1}'.format(i+1, len(patches)))

            # Convert dx to float
            dx = dx.astype(np.float32)

            origin = np.array(patch - dx / 2. - pad, dtype=np.int64)
            dx_eps = dx + float(2 * pad)

            # Convert dx_eps to int
            dx_eps = dx_eps.astype(np.int32)

            # Initialise
            delta = None
            velbx = None
            velby = None
            velbz = None
            vbc = None

            vbc_utils.msg(rank, "Loading patch: {0}.".format(patch), verbose)
            # vbc_utils.msg(rank, "TESTING: origin {0}, dx_eps {1}".format(origin, dx_eps), verbose)
            delta = ics[0].load_patch(origin, dx_eps)
            velbx = ics[1].load_patch(origin, dx_eps)
            velby = ics[2].load_patch(origin, dx_eps)
            velbz = ics[3].load_patch(origin, dx_eps)
            vbc = ics[4].load_patch(origin, dx_eps)

            # TESTING
            # vbc_utils.vbc_patch_dist(vbc, origin)
            # continue
            
            # Add to \delta sum for correction
            x_shape, y_shape, z_shape = delta.shape
            deltab_tmp = np.sum(delta[0 + pad:x_shape - pad,
                                      0 + pad:y_shape - pad,
                                      0 + pad:z_shape - pad] + 1.)
            
            deltab_tot = deltab_tot + deltab_tmp
            
            # Compute the bias
            vbc_utils.msg(rank, "Computing bias.", verbose)
            # Commented the below for testing
            k, b_c, b_b, b_vc, b_vb = vbc_utils.compute_bias(ics[4], vbc)

            # Convolve with field
            vbc_utils.msg(rank, "Performing convolution.", verbose)
            delta_biased = vbc_utils.apply_density_bias(ics[0], k, b_b, delta.shape[0], delta_x=delta)
            velbx_biased = vbc_utils.apply_density_bias(ics[1], k, b_vb, velbx.shape[0], delta_x=velbx)
            velby_biased = vbc_utils.apply_density_bias(ics[2], k, b_vb, velby.shape[0], delta_x=velby)
            velbz_biased = vbc_utils.apply_density_bias(ics[3], k, b_vb, velbz.shape[0], delta_x=velbz)

            # print('deltab before', delta_biased)
            

            # Add to \delta_biased sum for correction
            deltab_b_tmp = np.sum(delta_biased[0 + pad:x_shape - pad,
                                               0 + pad:y_shape - pad,
                                               0 + pad:z_shape - pad] + 1.)
            deltab_b_tot = deltab_b_tot + deltab_b_tmp
            
            # Remove the padded region
            x_shape, y_shape, z_shape = delta_biased.shape
            delta_biased = delta_biased[0 + pad:x_shape - pad,
                                        0 + pad:y_shape - pad,
                                        0 + pad:z_shape - pad]

            x_shape, y_shape, z_shape = velbx_biased.shape
            velbx_biased = velbx_biased[0 + pad:x_shape - pad,
                                        0 + pad:y_shape - pad,
                                        0 + pad:z_shape - pad]

            x_shape, y_shape, z_shape = velby_biased.shape
            velby_biased = velby_biased[0 + pad:x_shape - pad,
                                        0 + pad:y_shape - pad,
                                        0 + pad:z_shape - pad]

            x_shape, y_shape, z_shape = velbz_biased.shape
            velbz_biased = velbz_biased[0 + pad:x_shape - pad,
                                        0 + pad:y_shape - pad,
                                        0 + pad:z_shape - pad]

            if ret_vbc:
                vbc_patch = np.ones_like(delta_biased) * vbc_utils.vbc_rms(vbc)

                # Store
                biased_patches = [Patch(patch, dx, delta_biased, 'deltab'),
                                  Patch(patch, dx, velbx_biased, 'velbx'),
                                  Patch(patch, dx, velby_biased, 'velby'),
                                  Patch(patch, dx, velbz_biased, 'velbz'),
                                  Patch(patch, dx, vbc_patch, 'vbc_patch')]
            else:
                # Store
                biased_patches = [Patch(patch, dx, delta_biased, 'deltab'),
                                  Patch(patch, dx, velbx_biased, 'velbx'),
                                  Patch(patch, dx, velby_biased, 'velby'),
                                  Patch(patch, dx, velbz_biased, 'velbz')]
            
            for patch in biased_patches:
                with open(r"patches/level_{0:03d}/patch_{1}{2:03d}.p".format(level, patch.field, rank), "ab") as f:
                    pickle.dump(patch, f)

            del vbc
            del delta
            del delta_biased
            del velbx_biased
            del velby_biased
            del velbz_biased
            if ret_vbc: del vbc_patch

            gc.collect()

    del biased_patches
    gc.collect()

    vbc_utils.msg(rank, 'Done patches', verbose)
    # Wait until everyone has done patches
    barrier()
    # finalize()

    # Gather all of the \delta sums and reduce
    deltab_tot_sum = comm.reduce(deltab_tot, op=MPI.SUM, root=0)
    deltab_b_tot_sum = comm.reduce(deltab_b_tot, op=MPI.SUM, root=0)
    
    # Then write
    if rank == 0:
        # Save the correction as a small text file, use like:
        # \deltab_biased + correction = \deltab (what about the +1s?)
        np.savetxt("./patches/level_{0:03d}/"
                   "deltab_correction.txt".format(level),
                   np.array([deltab_tot_sum - deltab_b_tot_sum]))
        
        vbc_utils.msg(rank, 'Writing patches', verbose)
        write(path, level, lin, verbose, ret_vbc, comm=comm)

        
def write(path, level, lin, verbose=True, ret_vbc=False, comm=None):
    """Writes the patches computed by 'work' to new grafic IC files

    :param path: (str) path to (unbiased) grafic fields
    :param level: (int) level of ICs to work on
    :param lin: (bool) only write deltac/deltab (True) or deltab/velb* (False)?
    :param verbose: (bool) controls printout
    :returns: 
    :rtype:

    """

    # MPI stuff
    if comm is None:
        comm = MPI.COMM_WORLD
        
    size = comm.Get_size()
    rank = comm.Get_rank()
    barrier = comm.Barrier
    finalize = MPI.Finalize

    if rank == 0:
        # Make patches dir
        if not os.path.isdir("./patches/level_{0:03d}".format(level)):
            raise Exception("'patches' directory does not exist for this level. Run in 'work' mode first.")

        vbc_utils.msg(rank, "Writing fields.", verbose)

        if lin:
            # In the linear case we will only write out the delta
            # fields, and the other fields will have to be written
            # using gen_ics
            fields = ['deltab', 'deltac']
        else:
            # Otherwise write out all the baryon fields
            fields = ['deltab', 'velbx', 'velby', 'velbz']

        if ret_vbc: fields.append('vbc_patch')

        # ics = [grafic.load_snapshot(path, level, field=field) for field in fields]
        ics = [grafic.load_snapshot(path, level, field='deltab')]

        # Loop over fields
        for field in fields:
            # Get all of the patches for each field name, assumes that
            # patches will make up the entirety of the pickled objects
            # in the patches directory
            fns = glob.glob('./patches/level_{0:03d}/*{1}*.p'.format(level,
                                                                     field))
            fns.sort()
            size = len(fns)
            
            # Write new ICs
            output_field = np.zeros((ics[0].n[1], ics[0].n[0], ics[0].n[2]),
                                    dtype=np.float32)

            dest = []
            for i, fn in enumerate(fns):
                # Unpickle
                with open(fn, "rb") as f:
                    vbc_utils.msg(rank, 'Loading {0} pickle [{1}/{2}]'.format(field, i+1, size), verbose)
                    while True:
                        try:
                            dest.append(pickle.load(f))
                        except EOFError:
                            break

            # Load up correction if we're working on the deltab field
            # and if we have a correction
            if field == 'deltab':
                try:
                    c_fn = './patches/level_{0:03d}/deltab_correction.txt'
                    correction = np.loadtxt(c_fn.format(level))
                    # Calculate each cell's contribution to the correction
                    correction = np.float64(correction) / (ics[0].n[1] *
                                                           ics[0].n[0] *
                                                           ics[0].n[2])
                except:
                    correction = 0.0

                vbc_utils.msg(rank, 'Correcting deltab by {0:.3e}/cell'.format(correction))

            # Move from the list of patches to an array
            for item in dest:
                patch = item.patch
                dx = item.dx
                biased = item.data
                field = item.field

                # Bounds of this patch
                x_min, x_max = (int((patch[0]) - (dx[0] / 2.)), int((patch[0]) + (dx[0] / 2.)))
                y_min, y_max = (int((patch[1]) - (dx[1] / 2.)), int((patch[1]) + (dx[1] / 2.)))
                z_min, z_max = (int((patch[2]) - (dx[2] / 2.)), int((patch[2]) + (dx[2] / 2.)))

                # Apply the correction -- since it is going to be
                # quite small, we have to cast the biased array to
                # float64 then apply the correction, then revert to
                # float32 for output. Doing it piecemeal to avoid a
                # sudden jump in memory usage.
                if field == 'deltab':
                    print('TESTING:')
                    print('biased.sum() before correction', biased.sum())
                    biased = biased.astype(np.float64) + correction
                    print('biased.sum() after correction', biased.sum())
                    biased = biased.astype(np.float32)
                    print('biased.sum() after downcasting', biased.sum())

                # Place into output
                output_field[y_min:y_max, x_min:x_max, z_min:z_max] = biased

            # Write the initial conditions
            ics_dir = os.path.join(ics[0].level_dir, "ics_ramses_vbc/")
            if not os.path.isdir(ics_dir):
                os.mkdir(ics_dir)
            out_dir = os.path.join(ics_dir, "level_{0:03d}/".format(level))
            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)

            vbc_utils.msg(rank, 'Writing {0} field.'.format(field), verbose)
            ics[0].write_field(output_field, field, out_dir=out_dir)
            vbc_utils.msg(rank, 'Wrote {0} field.'.format(field), verbose)

            # TESTING
            # if field == 'deltab':
            #     ics[0].write_field(output_field+correction,
            #                        'deltab_corrected', out_dir=out_dir)

            
        # Remove patches/level_xxx dir
        # vbc_utils.clean(level)
        # vbc_utils.msg(rank, 'Cleaned up.')

        if lin:
            vbc_utils.msg(rank, 'Remember to use gen_ics to make the rest of the ICs!')

    # We have to wait until rank 0 has done the final reading and
    # writing, then everything can finish at the same time
    vbc_utils.msg(rank, 'Done!')
    # barrier()
    # finalize()
        

if __name__ == '__main__':
    import sys
    import traceback

    if len(sys.argv) < 4:
        print('Usage: [mpiexec] python bias_ics.py </path/to/ics/ (str)> '
              '<level (int)> <patch size (float)> <mode("work" or "write")>\n'
              '[<lin (bool)> <verbose (bool)>]', flush=True)
        sys.exit()

    path = sys.argv[1]
    level = int(sys.argv[2])
    patch_size = float(sys.argv[3])
    levelmin = int(sys.argv[4])
    mode = str(sys.argv[5])
    lin = False
    verbose = True

    # Optional linear argument
    if len(sys.argv) > 6:
        lin = bool(int(sys.argv[6]))
    # Optional verbose argument
    if len(sys.argv) > 7:
        verbose = bool(int(sys.argv[7]))

    ret_vbc = True

    # Run the main loop
    try:    
        if mode == 'work':
            # In 'work' mode, we do the calculations then write the
            # fields. Not terribly efficient, as the other processes
            # are waiting for rank 0 to write, but the writing is
            # quite quick anyway.
            work(path, level, patch_size, levelmin, lin, verbose, ret_vbc)
            # write(path, level, lin, verbose, ret_vbc)
        elif mode == 'write':
            # In 'write' mode, we just write
            write(path, level, lin, verbose, ret_vbc)
        else:
            raise Exception("mode is {0} -- should be 'work' or 'write'.".format(mode))

    except:
        # Catch exceptions
        from mpi4py import MPI
        
        print(traceback.format_exc(), flush=True)
        MPI.COMM_WORLD.Abort(96)

