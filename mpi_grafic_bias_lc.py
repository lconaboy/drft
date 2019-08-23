import numpy as np
import os
import gc
import pickle
"""
TODO

- swap ics from a list to a dictionary for more explicit indexing

- figure out padding stuff -- currently there is a region of size pad
around the edge that isn't included, as the periodic access of data
currently isn't implemented
 """

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
    from mpi4py import MPI
    import gc
    import os
    import time
  
    # MPI stuff
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    msg = 'rank {0}: {1}'

    #mpi.msg("Loading initial conditions")
    print(msg.format(rank, "Loading initial conditions"))
    
    if rank == 0:
        # Make sure vbc field exists on disk
        if not os.path.isfile(path+"level_{0:03d}/ic_vbc".format(level)):
            grafic_snapshot.derive_vbc(path, level)
        # Make patches dir
        if os.path.isdir("./patches"):
            raise Exception('Patches already exist. Remove and re-run.')
        elif not os.path.isdir("./patches"):
            os.mkdir("./patches")
    else:
        # Wait for rank 0 to write the velb, velc and vbc fields
        while not os.path.isfile(path+"level_{0:03d}/ic_vbc".format(level)):
            time.sleep(1.0e-3)

        # ics = grafic_snapshot.load_snapshot(path, level, sample_fft_spacing=False)
    ics = [grafic_snapshot.load_snapshot(path, level, field=field) for field in ['deltab', 'vbc']]

    pad = 8
    div = np.array([float(i) for i in divisors(ics[0].N, mode='yield')])
    idx = np.abs((ics[0].N / div) * ics[0].dx - patch_size).argmin()
    ncubes = int(div[idx])

    # Compute cube positions in cell units
    cubes = dx = None
    if rank == 0:
        print(msg.format(rank, "Using {0} cubes per dimension.".format(ncubes)))

        cubes, dx = vbc_utils.cube_positions(ics[0], ncubes, ics[0].N)
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

        delta = vbc = None
        if (P): print(msg.format(rank, "Loading patch: {0}").format(patch))
        delta = ics[0].load_patch(origin, int(dx_eps))
        vbc = ics[1].load_patch(origin, int(dx_eps))

        # Compute the bias
        if (B): print(msg.format(rank, "Computing bias"))
        k_bias, b_cdm, b_b = vbc_utils.compute_bias_lc(ics[1], vbc)

        # Convolve with field
        if (C): print(msg.format(rank, "Performing convolution"))
        delta_biased = vbc_utils.apply_density_bias(ics[0], k_bias, b_b, delta.shape[0], delta_x=delta)

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

    del biased_patch
    gc.collect()
    
    print(msg.format(rank, 'Done patches'))

############################## END OF WORK LOOP ###############################
    if rank == 0:
        import os
        # Write new ICs

        # output_field = np.zeros(ics[0].n)

        # To account for the non-periodic biasing
        output_field = ics[0].load_patch([0, 0, 0], ics[0].N)

        dest = []
        for i in range(size):
            # Unpickle
            with open(r"patches/patch_{0}.p".format(i), "rb") as f:
                print('Loading pickle {0}/{1}'.format(i, size), end='\r')
                while True:
                    try:
                        dest.append(pickle.load(f))
                    except EOFError:
                        break

        # print('Number of elements {0} (zero: {1})'.format(len(dest), np.count_nonzero(dest==0)))

            
        for item in dest:
            patch = item.patch
            dx = item.dx
            delta_biased = item.field
            
            # Bounds of this patch
            x_min, x_max = (int((patch[0]) - (dx / 2.)), int((patch[0]) + (dx / 2.)))
            y_min, y_max = (int((patch[1]) - (dx / 2.)), int((patch[1]) + (dx / 2.)))
            z_min, z_max = (int((patch[2]) - (dx / 2.)), int((patch[2]) + (dx / 2.)))

            print('bounds = [{0}, {1}, {2}]'.format([x_min, x_max], [y_min, y_max], [z_min, z_max]))
            print('delta_biased = {0}'.format(delta_biased))
            # Place into output
            output_field[x_min:x_max, y_min:y_max, z_min:z_max] = delta_biased

        # Write the initial conditions
        ics_dir = "{0}/ics_ramses_vbc/".format(ics[0].level_dir)
        if not os.path.isdir(ics_dir):
            os.mkdir(ics_dir)
        out_dir = "{0}/level_{1:03d}/".format(ics_dir, level)
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        print('np.all(output_field == 0) = {0}'.format(np.all(output_field == 0.0)))
        ics[0].write_field(output_field, "deltab", out_dir=out_dir)

if __name__ == "__main__":
    import sys
    import traceback
    
    path = sys.argv[1]
    level = int(sys.argv[2])
    patch_size = float(sys.argv[3])

    main(path, level, patch_size)

    print("Done!")
