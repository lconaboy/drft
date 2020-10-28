import numpy as np

"""Tools for reading in grafic files. A Snapshot instance won't
contain any data, but the data can be read in from the instance.

Example

    from grafic_tools import Snapshot

    path = '/path/to/ics/dir/'  # not including level_xxx/
    level = 6
    field = deltab'             # not including 'ic_'
    s = Snapshot(path, level, field)
    data = s.load_box()

TODO 

    - add case where a grafic Cube is being read, but
      scipy.FortranFile is unavailable

"""


class Header:

    def __init__(self, path):
        """Class for reading the header of grafic IC files. Attributes are
        the number of bytes in the header, number of bytes in each slice,
        tuple of number of points in each dimension, tuple of offsets in
        each dimension, spacing between grid points and cosmological
        parameters.

        Example
          h = Header(path)
       
        :param path: 
            (str) path to grafic file
        :returns: 
            None
        :rtype: 
            NoneType
        """

        # Read header
        n, dx, xoff, cosmo, size, area = self.read_header(path)

        # Turn into attributes
        self.n = n
        self.dx = dx
        self.xoff = xoff
        self.cosmo = cosmo
        self.size = size
        self.area = area


    def read_header(self, path):
        """Function for taking the path to a grafic IC file and reading the
        header.

        :param path: 
            (str) path to IC file
        :returns: 
            (3-tuple) number of grid points, (float) initial separation
            of points, (3-tuple) offset of grids (i.e. coordinates), 
            (4-tuple) cosmological parameters, (int) total number of
            bytes in header, (int) number of bytes in a slice
        :rtype:

        """
        
        with open(path, 'r+') as f:
            # Read the first 4 bytes written by Fortran
            d = np.fromfile(f, dtype=np.int32, count=1)

            # The header has the value 44, so we can use this
            # to check we've read to the right position
            assert(d == 44), 'header = {0}'.format(int(d))
    
            # Now read the header info
            (n1, n2, n3) = np.fromfile(f, dtype=np.int32, count=3)
            (dxini, xoff1, xoff2, xoff3, astart, omega_m, omega_l, h0) = np.fromfile(f, dtype=np.float32, count=8)

            # Read the final 4 bytes written by Fortran, again these
            # should have the value 44
            d = np.fromfile(f, dtype=np.int32, count=1)
            assert(d == 44), 'footer = {0}'.format(int(d))

            # # Finally, there is an integer specifying the total number
            # # of bytes in a slice (although I've read that this might
            # # be compiler dependent, so perhaps not a good idea to
            # # rely on this...)
            # area = np.fromfile(f, dtype=np.int32, count=1)
            # assert(area == int(n1 * n2 * 4)), 'area = {0} bytes'.format(area)

            # Wrap up variables conveniently, converting h0 from H0 to
            # little h
            n = (n1, n2, n3)
            dx = dxini
            xoff = (xoff1, xoff2, xoff3)
            cosmo = (astart, omega_m, omega_l, h0/100.0)
            area = int(n1 * n2 * 4)  # Total number of bytes in a slice
            size = 44 + 8            # Total number of bytes in header,
                                     # including initial and end record markers

            return n, dx, xoff, cosmo, size, area

class Snapshot:

    def __init__(self, path, level, field):
        """Class for containing and operating on a grafic IC file. 

        :param path: 
            (str) path to IC directory
        :param level: 
            (int) level to load
        :param field: 
            (str) field to load (i.e. 'deltab')
        :returns: 
            snapshot instance
        :rtype: 
            (Snapshot)

        """

        self.path = path
        self.level = level
        self.field = field
        self.fname = self.ic_fname()

        self.header = Header(self.fname)
        self.area = self.header.area
        self.size = self.header.size
        self.n = self.header.n
        self.N = self.n[0]
        self.xoff = self.header.xoff
        self.dx = self.header.dx
        self.z = 1.0/self.header.cosmo[0] - 1.0

        # h here is little h, i.e. H0 = 100*h km/s/Mpc
        self.cosmo = {"aexp":self.header.cosmo[0],
                      "omega_m":self.header.cosmo[1],
                      "omega_l": self.header.cosmo[2],
                      "h":self.header.cosmo[3], "z":self.z}


        self.boxsize = self.N * self.dx


    def load_patch(self, origin, N):
        """Partially load a grafic file. Allows for periodic boundary
        conditions, i.e. origin + N can be > n.

        :param origin: 
            (seq of ints) coordinates of lower left corner of region to load
        :param N:
            (seq of ints) extent of region
        :returns:
            value of the grafic file for the specified region
        :rtype:
            (array)

        """
        # Retrieve the number of bytes in the header, and each slab
        area = self.area
        size = self.size

        fname = self.ic_fname()
        
        # Retrieve the number of points in each slice
        (n1, n2, n3) = self.n

        # Store the data in a patch, swapping the row and column for
        # loading, see comments passim ad nauseam
        patch = np.zeros(shape=(N[1], N[0], N[2]), dtype=np.float32)

        with open(fname, "rb") as f:
            # Seek past the header block to between the final header
            # record marker and the inital data record marker
            f.seek(size, 0)

            # Seek to the first relevant slab
            z0 = int(origin[2] % n3)
            f.seek(z0 * (area + 8), 1)
            
            # The z-axis changes the slowest
            for i3 in range(origin[2], origin[2] + N[2]):
                # Check that we haven't looped back through the
                # file. First need to normalise f.tell() to n_slabs.
                cur_pos = (f.tell() - (self.size)) / (self.area + 8)
                if cur_pos >= n3:
                    f.seek(size, 0)
            
                # Read the initial record marker, this provides a good
                # check that we're in the right place
                rm = np.fromfile(f, dtype=np.int32, count=1)
                assert rm == self.area, 'For slab {0} expected initial record marker of {1} but got {2}'.format(i3, self.area, rm)
    
                # Pick out the slab, and reshape to (y, x). The slab
                # is transposed because Python is row-major whereas
                # Fortran is column-major, or we can use row-major as
                # long as we're consistent about writing out in
                # row-major
                # slab = np.fromfile(f, dtype=np.float32,
                #                    count=(n1 * n2)).reshape((n1, n2)).transpose()
                slab = np.fromfile(f, dtype=np.float32,
                                    count=(n1 * n2)).reshape((n2, n1), order='F')

                # Allow for periodic boundary conditions and the
                # transposed array
                i1 = np.arange(origin[0], origin[0]+N[0])
                i2 = np.arange(origin[1], origin[1]+N[1])
                x = np.mod(i1, n1, dtype=int)
                y = np.mod(i2, n2, dtype=int)
                xg, yg = np.meshgrid(x, y)

                iz = int(i3 - origin[2])

                # Again, the choice of slab indices depends on whether
                # we are using row-major or column-major. Row-major is
                # fine as long as we are consistent about writing out
                # in row-major.
                # patch[:, :, iz] = slab[xg, yg]
                patch[:, :, iz] = slab[yg, xg]

                # Read the end record marker, this provides a good
                # check that we're in the right place
                rm = np.fromfile(f, dtype=np.int32, count=1)
                assert rm == self.area, 'For slab {0} expected end record marker of {1} but got {2}'.format(i3, self.area, rm)

        return patch


    def load_box(self):
        """Loads the entire grafic file into memory.

        :returns: 
            values for that snapshot
        :rtype: 
            (array)

        """
        
        origin = np.array([0, 0, 0])
        N = self.N

        area = self.area
        size = self.size

        fname = self.ic_fname()
        
        # Retrieve the number of points in each slice
        (n1, n2, n3) = self.n

        # This is a bit inefficient but at least it is consistent
        box = self.load_patch(origin, [n2, n1, n3])

        # Store the data in a patch, n2 comes before n1 because
        # Fortran is column-major where Python is row-major, or if we
        # write out in row-major we can read in in row-major and get
        # the correct answer
        # box = np.zeros(shape=(n2, n1, n3), dtype=np.float32)
        # box = np.zeros(shape=(n1, n2, n3), dtype=np.float32)

        # z0 = origin[2]
        # z = z0

        # with open(fname, "rb") as f:
        #     # Seek past the header block to between the final header
        #     # record marker and the inital data record marker
        #     f.seek(size, 0)

        #     # The z-axis changes the slowest
        #     for iz in range(n3):
        #         # Read the initial record marker, this provides a good
        #         # check that we're in the right place
        #         rm = np.fromfile(f, dtype=np.int32, count=1)
        #         assert rm == self.area, 'For slab {0} expected initial record marker of {1} but got {2}'.format(iz, self.area, rm)
                
        #         # Pick out the plane, and reshape to (y, x) since
        #         # Python is row-major where Fortran is column-major
        #         # or, as long as we are consistent about writing out in
        #         # row-major format, then we can read in in row-major
        #         slab = np.fromfile(f, dtype=np.float32,
        #                            count=(n1 * n2)).reshape((n1, n2))
        #         # slab = np.fromfile(f, dtype=np.float32,
        #         #                    count=(n1 * n2)).reshape(n1, n2).transpose()

        #         # Store the slab in the box
        #         box[:, :, iz] = slab

        #         # Read the final record marker, this provides a good
        #         # check that we're in the right place
        #         rm = np.fromfile(f, dtype=np.int32, count=1)
        #         assert rm == self.area, 'For slab {0} expected end record marker of {1} but got {2}'.format(iz, self.area, rm)

        return box
        
    def plot_slice(self, nslc=0, ret_im=False):
        #from FortranFile import FortranFile
        from scipy.io import FortranFile
        import matplotlib.pyplot as plt
        
        fname = self.ic_fname()

        f = FortranFile(fname)
        # Seek past the header block to between the final header
        # record marker and the inital data record marker
        # f.seek(size, 0)

        # slc = f.read_reals(float).reshape((self.n[0], self.n[1]), order='F')
        f.read_record('i4')  # read header, they aren't all ints but we don't need them anyway

        # Skip first n slices
        for i in range(nslc):
            f.read_record('f4')

        slc = f.read_reals('f4').reshape((self.n[0], self.n[1]), order='F')

        fig, ax = plt.subplots(figsize=(6, 6))
        # imsave('{0}{1}_{2}_slice.png'.format(self.path, self.field, self.level), slc, cmap='gist_stern')
        im = ax.imshow(slc, cmap='gist_stern')
        plt.colorbar(im)
        plt.savefig('{0}{1}_{2}_slice.pdf'.format(self.path, self.field, self.level))
        plt.close('all')
        

        if ret_im: return slc


    def ic_fname(self):
        """Generates the filename string for the IC files, given the path,
        level and field supplied to Snapshot.

        :returns: filename
        :rtype: str

        """
        
        return self.path+'level_{0:03d}/ic_{1}'.format(self.level, self.field)


    def write_header(self, f):
        """Writes a header in the style of grafic headers, i.e. one consistent
        with Fortran unformatted binary files.

        :param f: (file object) the file object to write the header to
        :returns:
        :rtype:

        """
        
        (n1, n2, n3) = self.header.n
        n = np.array([n1, n2, n3], dtype=np.int32)
        dx = self.dx  # Mpc a 
        origin = self.xoff
        cosmo = self.cosmo

        header_bytes = np.array([44], dtype=np.int32)

        # We can use the numpy.ndarray.tofile functionality to easily
        # write the header while keeping control over data types

        header_bytes.tofile(f)
        n.tofile(f)
        np.array([dx, origin[0], origin[1], origin[2], cosmo['aexp'],
                  cosmo['omega_m'], cosmo['omega_l'],
                  cosmo['h'] * 100.0], dtype=np.float32).tofile(f)
        header_bytes.tofile(f)


    def write_field(self, data, field_name, out_dir=None, **kwargs):
        """Function for writing the data as a Fortran unformatted
        binary. Writes a header to mimic grafic files, as well as
        appending and prepending each 'record' (i.e. each [ix, :, :])
        with the number of bytes in that slab, as done by Fortran
        (when compiled with GCC). Assumes the data is from Python,
        i.e. that it is in row-major format. Fortran uses column-major
        format, so each slab will be transposed before writing.

        :param data: 
            (numpy.ndarray) data should be in np.float32 format and 
            have (n2, n1, n3) elements
        :param field_name: 
            (str) field to write as, will be prepended with 'ic_' 
        :param out_dir: 
            (str) directory to write the field to
        :returns: 
            None
        :rtype:

        """
        import os
        # Create filename
        fname = self.field_fname(field_name)
        if out_dir is not None:
            fname = '{0}/ic_{1}'.format(out_dir, field_name)

        # Check if file already exists
        if (os.path.isfile(fname)):
            raise Exception("Refusing to overwrite field {0}.".format(field_name))
        # if not, open the file for writing
        else:
            with open(fname, "w+") as f:
                # Write the data header
                self.write_header(f)

                # Get data shape
                (n1, n2, n3) = data.shape

                # To be consistent with the GCC I/O, each Fortran
                # binary has a header and a footer specifying the
                # number of bytes in the pre-/proceeding of each
                # record
                area = np.array([n1 * n2 * 4], dtype=np.int32)

                for iz in range(n3):
                    # Initial record marker
                    area.tofile(f)
                    
                    # Write data, converting back to column-major as
                    # is used by Fortran
                    #
                    slab = data[:, :, iz].flatten(order='F').astype(np.float32)
                    #
                    # or, as long as we are consistent about reading
                    # in in row-major format, then writing out in
                    # row-major will yield the correct output field
                    # slab = data[:, :, iz].flatten(order='C').astype(np.float32)
                    
                    slab.tofile(f)
                    # End record marker
                    area.tofile(f)

    @property
    def level_dir(self):
        return '{0}/level_{1:03d}/'.format(self.path, self.level)


    def field_exists_on_disk(self, field):
        '''
        Checks if field is written to disk
        '''
        import os
        fname = self.field_fname(field)
        return os.path.isfile(fname)


    def field_fname(self, field):
        return '%s/ic_%s' % (self.level_dir, field)


class Cube:    
    
    def __init__(self, path, grafic=False):
        """Class for reading in the result of using the modified RAMSES
        utils/f90/*2cube routines. I modified the routines to output
        aexp as a header.

        :param path: 
            (str) path to cube
        :param grafic:
            (bool) use grafic file format?
        :returns: 
            snapshot instance
        :rtype: 
            (Snapshot)

        """

        self.path = path
        self.grafic = grafic
        self.aexp = self.read_header()

    def load_box(self):
        """Loads the entire grafic file into memory.

        :returns: 
            values for that snapshot
        :rtype: 
            (array)

        """
        try:
            from scipy.io import FortranFile
            f = FortranFile(self.path, 'r')
            print('---- using FortranFile')
            # Read in grafic format, good for files over 2GB, where
            # the binary record would get split up into sub-records
            if self.grafic:
                print('-------- reading in grafic format')
                h = f.read_record(dtype=np.int32)  # Not bothered about
                                                   # aexp, so read all as
                                                   # ints
                box = np.zeros(shape=(h[0], h[1], h[2]))

                for k in range(h[2]):
                    box[:, :, k] = f.read_record(dtype=np.float32).reshape((h[0], h[1]), order='F')     
                print('---- read in cube of size ({0}, {1}, {2})'.format(h[0], h[1], h[2]))
                
            else:
                print('-------- reading in binary format')
                h = f.read_record(dtype=np.int32)  # Not bothered about
                                                   # aexp, so read all as
                                                   # ints
                print(h)
                print('---- read in cube of size ({0}, {1}, {2})'.format(h[1], h[2], h[3]))
                                             
                box = f.read_record(dtype=np.float32).reshape((h[1], h[2], h[3]), order='F')
            
        except:
            # Quietly fail
            if self.grafic: return
            
            with open(self.path, "rb") as f:
                # Skip past the initial record marker and the aexp header
                hi = np.fromfile(f, dtype=np.uint32, count=1) # f.seek(8, 0)
                aexp = np.fromfile(f, dtype=np.float32, count=1) # f.seek(8, 0)
                # Read in size of cube
                nx, ny, nz = np.fromfile(f, dtype=np.int32, count=3)
                print('Read in cube of size ({0}, {1}, {2})'.format(nx, ny, nz))
                # Skip past the final record marker
                hf = np.fromfile(f, dtype=np.uint32, count=1) # f.seek(8, 0) # f.seek(4, 1)
                hi = np.fromfile(f, dtype=np.uint32, count=1) # f.seek(8, 0) # f.seek(4, 1)     
                box = np.fromfile(f, dtype=np.float32,
                                  count=(nx * ny * nz)).reshape(nx, ny, nz)
                hf = np.fromfile(f, dtype=np.uint32, count=1)
        return box

    def read_header(self):
        try:
            from scipy.io import FortranFile
            f = FortranFile(self.path, 'r')
            h = f.read_record(dtype=np.float32)  # Not bothered about
                                                 # n*, so read all as
                                                 # floats
            if self.grafic:
                aexp = h[7]
            else:
                aexp = h[0]

        except:
            # Quietly fail
            if self.grafic: return
            
            with open(self.path, "rb") as f:
                # Skip past initial record marker
                f.seek(4, 0)
                aexp = np.fromfile(f, dtype=np.float32, count=1)
                print('Output at a = ', aexp)

        return aexp

    

def load_snapshot(path, level, field):
    """Convenience function for loading grafic files.

    :param path: 
        (str) path to directory containing the level_... directories
    :param level: 
        (int) level to load
    :param field: 
        (str) field to load, i.e. 'deltab' to load file 'ic_deltab'
    :returns: 
        the grafic file in a Snapshot class
    :rtype: 
        (Snapshot)

    """
    return Snapshot(path, level, field)


def load_cube(path):
    """Convenience function for loading *2cube outputs

    :param path:
        (str) path to *2cube output

    """
    return Cube(path)


def grid_velc(path, level):
    """Grids the CDM peculiar velocity fields.

    :param path:
        (str) path to ic_velc* fields
    :param level:
        (int) level of the fields
    :returns: 
        None
    :rtype: 
        NoneType

    """
    from interp_part import cic

    # First create the fields
    coords = ['x', 'y', 'z']
    qty = ['pos', 'vel']
    fields = [['{0}c{1}'.format(q, c) for c in coords] for q in qty]
    print(fields)
    
    # Load up each snapshot
    s = [[load_snapshot(path, level, f) for f in field] for field in fields]

    # Get the sizes of the IC files
    N = s[0][0].N
    assert(s[0][0].N == s[0][1].N == s[0][2].N), 'posx, posy and posz files different sizes.'
    assert(s[0][0].N == s[1][1].N == s[1][2].N), 'posx, vely and velz files different sizes.'
    
    # Interpolate
    velg = [cic(s[1][i].load_box(), s[0][i].load_box()) for i in range(3)]

    # Now write this data to a field
    out_dir = s[0][0].level_dir
    for i in range(3):
        s[0][0].write_field(velg[i], 'velcg{0}'.format(coords[i]), out_dir=out_dir)


def derive_vel(path, level, species):
    """Calculates the magnitude of the velocity fields for the specified
    species. Then writes the ic_velb and ic_velcg files.

    :param path: (str) path to directory containing the level_... directories
    :param level: (int) level to load
    :param species: (str) should be either 'b' or 'cg'
    :returns: None
    :rtype:

    """
    # Check that species is correct
    assert species in ['b', 'c'], "Species should be either 'b' or 'c'"
    
    # First create the fields
    fields = ['vel{0}{1}'.format(species, c) for c in ['x', 'y', 'z']]

    # Load up each snapshot
    s = [load_snapshot(path, level, field) for field in fields]

    # Get the sizes of the IC files
    N = s[0].N
    assert(s[0].N == s[1].N == s[2].N), 'velx, vely and velz files different sizes.'
    
    # Now extract the values from each snapshot
    vx = s[0].load_box()
    vy = s[1].load_box()
    vz = s[2].load_box()

    # Calculate velocity
    v = np.sqrt(vx**2 + vy**2 + vz**2)

    # These files may be quite large
    del vx
    del vy
    del vz

    # Now write this data to a field
    out_dir = s[0].level_dir
    s[0].write_field(v, 'vel{0}'.format(species), out_dir=out_dir)
    

def derive_vbc(path, level, per):
    """Calculates the velb, gridded velc and vbc fields using the
    subroutine velc from cic.f95. Make sure that cic.f95 has been
    compiled with f2py before running.

    :param path: 
        (str) path to directory containing the level_... directories
    :param level: 
        (int) level to load
    :returns: 
        None
    :rtype:
        NoneType

    """
    import os
    import cic
    import time
    import warnings

    level_path = path+'level_{0:03d}/'.format(level)

    # Check that none of the fields exist already. cic just overwrites
    # the existing files, so this can be a warning instead.
    if os.path.isfile(level_path+'ic_vbc'):
        warnings.warn(level_path+'ic_vbc exists, not regenerating.')
        return

    # Initialise vbc
    s = load_snapshot(path, level, 'deltab')
    vbc = np.zeros(s.n)

    for c in ['x', 'y', 'z']:
        print('---- working on vel' + c)
        vb = load_snapshot(path, level, 'velb'+c).load_box()
        print('-------- read in velb' + c)
        print('DEBUGGING: min, max', vb.min(), vb.max())
        
        vc = load_snapshot(path, level, 'velc'+c).load_box()
        print('-------- read in velc' + c)
        print('DEBUGGING: min, max', vc.min(), vc.max())
        
        
        vbc += (vb - vc) ** 2.

    vbc = np.sqrt(vbc)
    print('---- writing vbc field')
    s.write_field(vbc, 'vbc')
    
    # # Wait up to 20 seconds for it to be written to disk
    # w = 0.0
    # while (not os.path.isfile(level_path+'ic_vbc')):
    #     if w < 30.0:
    #         time.sleep(0.5) # Wait half a second
    #     else:
    #         raise Exception('Waited 30 seconds for ic_vbc to be written to disk. Quitting.')


def derive_deltac(path, level, per, omega_b=None):
    """Calculates the deltac field by CIC interpolating the CDM particle
    positions, using the delc routine from cic.f95. Make sure that
    cic.f95 has been compiled with f2py before running.

    :param path: 
        (str) path to directory containing the level_... directories
    :param level: 
        (int) level to load
    :param omega_b:
        (float or None) value of baryon density parameter, if None this 
        is read from the py_vbc params
    :returns: 
        None
    :rtype:
        NoneType

    """
    import os
    import cic
    import warnings
    
    level_path = path+'level_{0:03d}/'.format(level)

    print('Deriving deltac field')
    
    # Check that none of the fields exist already, otherwise cic will
    # throw an error
    if os.path.isfile(level_path+'ic_deltac'):
        warnings.warn(level_path+'ic_deltac exists, not regenerating.')
        return

    # Get omega_b
    if omega_b is None:
        import configparser
        config_path, _ = os.path.split(__file__)
        config_fname = 'py_vbc/planck2018_params.ini'
        config = configparser.ConfigParser()
        config.read(config_path + '/' + config_fname)

        omega_b = np.float32(config.get('cosmology', 'omega_b'))
        print('---- using py_vbc params for omega_b: ', omega_b)
    elif type(omega_b) == type(0.1):
        omega_b = np.float32(omega_b)
        print('---- using provided value for omega_b: ', omega_b)        
    else:
        raise Exception('omega_b should be float or None')
        
    cic.gen_delc(level_path, omega_b, per)

    # Wait up to 20 seconds for it to be written to disk
    w = 0.0
    while (not os.path.isfile(level_path+'ic_deltac')):
        if w < 30.0:
            time.sleep(0.5) # Wait half a second
        else:
            raise Exception('Waited 30 seconds for ic_deltac to be written to disk. Quitting.')
    
