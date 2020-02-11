import os
import gc
import sys
import pickle
import numpy as np
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

    def __init__(self, patch, dx, data, field):
        self.patch = patch
        self.dx = dx
        self.data = data
        self.field = field


def main(path, level, patch_size):
    '''
    Writes a new set of grafIC initial conditions with a drift velocity dependent
    bias in the power spectrum
    '''
    import utils as vbc_utils
    from utils import divisors
    import grafic_tools as grafic
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
            print(msg.format(rank, 'Deriving vbc'))
            grafic.derive_vbc(path, level)
        # Make patches dir
        if os.path.isdir("./patches"):
            raise Exception('Patches already exist. Remove and re-run.')
        elif not os.path.isdir("./patches"):
            os.mkdir("./patches")
            print(msg.format(rank, 'Made patches directory.'))
    else:
        # Wait for rank 0 to write the velb, velc and vbc fields
        while not os.path.isfile(path+"level_{0:03d}/ic_vbc".format(level)):
            time.sleep(1.0e-3)

        # ics = grafic.load_snapshot(path, level, sample_fft_spacing=False)
    ics = [grafic.load_snapshot(path, level, field=field) for field in
           ['deltab', 'velbx', 'velby', 'velbz', 'vbc']]

    # This only works with cubic zoom regions, so check that and exit
    # if not
    if ((ics[0].n[0] != ics[0].n[1]) or (ics[0].n[0] != ics[0].n[2])):
        print('Only works with cubic zoom regions. Use force_equal_extent = yes with music to generate cubic zoom regions.')
        sys.exit()
    
    pad = 8

    # Calculate the number of cubes for each dimension, we want an
    # equal number of cubes in each dimension and this is most easily
    # achieved if the zoom region is actually cubic
    ncubes = np.zeros(3)
    for i in range(3):
        div = np.array([float(i) for i in divisors(ics[0].n[i], mode='yield')])
        idx = np.abs((ics[0].n[i] / div) * ics[0].dx - patch_size).argmin()
        ncubes[i] = int(div[idx])

    # Compute cube positions in cell units
    cubes = dx = None
    if rank == 0:
        print(msg.format(rank, "Using {0} cubes per dimension.".format(ncubes)))

        cubes, dx = vbc_utils.cube_positions(ics[0], ncubes, ics[0].N)
        cubes = np.array(cubes)
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

        if (P): print(msg.format(rank, "Loading patch: {0}").format(patch))
        delta = ics[0].load_patch(origin, dx_eps)
        velbx = ics[1].load_patch(origin, dx_eps)
        velby = ics[2].load_patch(origin, dx_eps)
        velbz = ics[3].load_patch(origin, dx_eps)
        vbc = ics[4].load_patch(origin, dx_eps)

        # Compute the bias
        if (B): print(msg.format(rank, "Computing bias"))
        k, b_c, b_b, b_vc, b_vb = vbc_utils.compute_bias(ics[4], vbc)

        # Convolve with field
        if (C): print(msg.format(rank, "Performing convolution"))
        delta_biased = vbc_utils.apply_density_bias(ics[0], k, b_b, delta.shape[0], delta_x=delta)
        velbx_biased = vbc_utils.apply_density_bias(ics[1], k, b_vb, velbx.shape[0], delta_x=velbx)
        velby_biased = vbc_utils.apply_density_bias(ics[2], k, b_vb, velby.shape[0], delta_x=velby)
        velbz_biased = vbc_utils.apply_density_bias(ics[3], k, b_vb, velbz.shape[0], delta_x=velbz)

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

        
        # Store
        biased_patches = [Patch(patch, dx, delta_biased, 'deltab'),
                          Patch(patch, dx, velbx_biased, 'velbx'),
                          Patch(patch, dx, velby_biased, 'velby'),
                          Patch(patch, dx, velbz_biased, 'velbz')]

        for patch in biased_patches:
            with open(r"patches/patch_{0}{1}.p".format(patch.field, rank), "ab") as f:
                pickle.dump(patch, f)

        del vbc
        del delta
        del delta_biased
        del velbx_biased
        del velby_biased
        del velbz_biased

        gc.collect()

    del biased_patches
    gc.collect()
    
    print(msg.format(rank, 'Done patches'))

############################## END OF WORK LOOP ###############################
    if rank == 0:
        import os
        # Loop over fields
        for field_name in ['deltab', 'velbx', 'velby', 'velbz']:
            # Write new ICs
            output_field = np.zeros(ics[0].n)

            dest = []
            for i in range(size):
                # Unpickle
                with open(r"patches/patch_{0}{1}.p".format(field_name, i), "rb") as f:
                    print(msg.format(rank, 'Loading {0} pickle [{1}/{2}]'.format(field_name, i+1, size)))
                    while True:
                        try:
                            dest.append(pickle.load(f))
                        except EOFError:
                            break

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

                # Place into output
                output_field[x_min:x_max, y_min:y_max, z_min:z_max] = biased

            # Write the initial conditions
            ics_dir = "{0}/ics_ramses_vbc/".format(ics[0].level_dir)
            if not os.path.isdir(ics_dir):
                os.mkdir(ics_dir)
            out_dir = "{0}/level_{1:03d}/".format(ics_dir, level)
            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)

            print(msg.format(rank, 'Writing {0} field'.format(field_name)))
            ics[0].write_field(output_field, field_name, out_dir=out_dir)
            print(msg.format(rank, 'Wrote {0} field'.format(field_name)))

if __name__ == "__main__":
    import sys
    import traceback

    from utils import clean
    
    path = sys.argv[1]
    level = int(sys.argv[2])
    patch_size = float(sys.argv[3])

    main(path, level, patch_size)

    clean()

    print("Done!")
