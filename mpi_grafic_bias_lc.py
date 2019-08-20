import numpy as np
import os
import gc
import pickle

# VERBOSE = 1  # 0 for all, >0 for just patch, <0 for none
P = False
B = False
C = False

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

    def __init__(self, patch, dx, field):
        self.patch = patch
        self.dx = dx
        self.field = field


def main(path, level, patch_size):
    '''
    Writes a new set of grafIC initial conditions with a drift velocity dependent
    bias in the power spectrum
    '''
    import utils as vbc_utils
    from utils import divisors
    import grafic_ics as grafic_snapshot
    #from seren3.analysis.parallel import mpi
    from mpi4py import MPI
    import gc
    import os
  
    # MPI stuff
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    msg = 'rank {0}: {1}'

    #mpi.msg("Loading initial conditions")
    print(msg.format(rank, "Loading initial conditions"))
    
    # ics = grafic_snapshot.load_snapshot(path, level, sample_fft_spacing=False)
    ics = grafic_snapshot.load_snapshot(path, level) #, field='deltab')

    if rank == 0:
        # Make sure vbc field exists on disk
        if not ics.field_exists_on_disk("vbc"):
            ics.write_field(ics["vbc"], "vbc")
        # Make patches dir
        if os.path.isdir("./patches"):
            raise Exception('Patches already exist. Remove and re-run.')
        if not os.path.isdir("./patches"):
            os.mkdir("./patches")

    pad = 8
    div = np.array([float(i) for i in divisors(ics.N - 2*pad, mode='yield')])
    idx = np.abs(((ics.N - 2*pad) / div) * ics.dx - patch_size).argmin()
    ncubes = int(div[idx])

    # Compute cube positions in cell units
    cubes = dx = None
    if rank == 0:
        print(msg.format(rank, "Using {0} cubes per dimension.".format(ncubes)))
        # Try creating a padded region around the outside of the IC box
        cubes, dx = vbc_utils.cube_positions(ics, ncubes, ics.N - 2*pad)
        cubes = np.array(cubes)
        cubes += pad
        # Split the cubes into chunks that can be scattered to each processor 
        chunks = np.array_split(cubes, size)
    else:
        chunks = None
    
    # dx is the same for every cube, so this is broadcast to each processor
    dx = comm.bcast(dx, root=0)

############################## WORK LOOP ######################################

    # Iterate over patch positions in parallel
    # dest = {}
    patches = comm.scatter(chunks, root=0)
    
    for i, patch in enumerate(patches):

        print(msg.format(rank, '{0}/{1}'.format(i+1, len(patches))))
    
        origin = np.array(patch - float(dx) / 2. - pad, dtype=np.int64)
        dx_eps = float(dx) + float(2 * pad)
        print('dx_eps = {0}'.format(dx_eps))

        delta = vbc = None
        if (P): print(msg.format(rank, "Loading patch: {0}").format(patch))
        delta = ics.load_patch("deltab", origin, int(dx_eps))
        vbc = ics.load_patch("vbc", origin, int(dx_eps))

        # Compute the bias
        if (B): print(msg.format(rank, "Computing bias"))
        k_bias, b_cdm, b_b = vbc_utils.compute_bias_lc(ics, vbc)

        # Convolve with field
        if (C): print(msg.format(rank, "Performing convolution"))
        delta_biased = vbc_utils.apply_density_bias(ics, k_bias, b_b, delta.shape[0], delta_x=delta)

        # Remove the padded region
        x_shape, y_shape, z_shape = delta_biased.shape
        delta_biased = delta_biased[0 + pad:x_shape - pad,
                                    0 + pad:y_shape - pad, 0 + pad:z_shape - pad]

        # Store
        biased_patch = Patch(patch, dx, delta_biased)
        
        with open(r"patches/patch_{0}.p".format(rank), "ab") as f:
            pickle.dump(biased_patch, f)

        del vbc
        del delta
        del delta_biased

        gc.collect()

        # if rank == 0:
            # print('after')
            # objgraph.show_growth(limit=3)
            # print(h.heap())
            # obj = objgraph.by_type('Patch')
            #objgraph.show_backrefs([obj], max_depth=10)


    del biased_patch
    gc.collect()
    
    # print(msg.format(rank, "memory usage = {1:.3f} Mo".format(rank, get_memory_usage())))
    print(msg.format(rank, 'Done patches'))

    # dest = comm.gather(biased_patches, root=0)
    
############################## END OF WORK LOOP ###############################
    if rank == 0:
        import os
        # Write new ICs

        output_field = np.zeros(ics.n)
        # dest = np.zeros(cubes.size, dtype=object)
        # j = 0
        dest = []
        for i in range(size):
            # Unpickle
            with open(r"patches/patch_{0}.p".format(i), "rb") as f:
                while True:
                    try:
                        # dest[j] = pickle.load(f)
                        # j += 1
                        dest.append(pickle.load(f))
                    except EOFError:
                        break

        print('Number of elements {0} (zero: {1})'.format(len(dest), np.count_nonzero(dest==0)))

            
        for item in dest:
            # result = item.result
            # patch = result["patch"]
            # dx = result["dx"]
            # delta_biased = result["field"]

            patch = item.patch
            dx = item.dx
            delta_biased = item.field
            
            # Bounds of this patch
            x_min, x_max = (int((patch[0]) - (dx / 2.)), int((patch[0]) + (dx / 2.)))
            y_min, y_max = (int((patch[1]) - (dx / 2.)), int((patch[1]) + (dx / 2.)))
            z_min, z_max = (int((patch[2]) - (dx / 2.)), int((patch[2]) + (dx / 2.)))

            # Place into output
            output_field[x_min:x_max, y_min:y_max, z_min:z_max] = delta_biased

        # Write the initial conditions
        ics_dir = "{0}/ics_ramses_vbc/".format(ics.level_dir)
        if not os.path.isdir(ics_dir):
            os.mkdir(ics_dir)
        out_dir = "{0}/level_{1:03d}/".format(ics_dir, level)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        ics.write_field(output_field, "deltab", out_dir=out_dir)

if __name__ == "__main__":
    import sys
    import traceback
    
    path = sys.argv[1]
    level = int(sys.argv[2])
    patch_size = float(sys.argv[3])

    #try:
    main(path, level, patch_size)
    # except Exception as e:
    #     from seren3.analysis.parallel import mpi
    #     mpi.msg("Caught exception (message): %s" % e.message)
    #     mpi.msg(traceback.format_exc())
    #     mpi.terminate(500, e=e)

    print("Done!")
