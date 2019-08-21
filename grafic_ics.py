import numpy as np

# level = 8
# N = 32
# # base = '/home/lc589/IC_2048z_14Mpc/level_{0:03d}/'.format(level)
# path = '/home/lc589/IC_2048z_14Mpc/'
# field = 'deltab'
# # path = base+field
# dx = 0.5**level

# origin = np.array([1, 1, 1], dtype=int)

class Header:

    def __init__(self, path):
        """Class for reading the header of grafic IC files. Attributes are
        the number of bytes in the header, number of bytes in each slice,
        tuple of number of points in each dimension, tuple of offsets in
        each dimension, spacing between grid points and cosmological
        parameters.

        Example
          h = Header(path)
       
        :param path: (str) path to grafic file
        :returns: None
        :rtype: NoneType
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

        :param path: (str) path to IC file
        :returns: (3-tuple) number of grid points, (float) initial separation
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
            size = 44 + 8 + 4        # Total number of bytes in header,
                                     # including first record header

            return n, dx, xoff, cosmo, size, area

class Snapshot:

    def __init__(self, path, level, field):
        """Class for containing and operating on a grafic IC file. 

        :param path: (str) path to IC directory
        :param level: (int) level to load
        :param field: (str) field to load (i.e. 'deltab')
        :returns: snapshot instance
        :rtype: Snapshot

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
                      "omega_M_0":self.header.cosmo[1],
                      "omega_lambda_0": self.header.cosmo[2],
                      "h":self.header.cosmo[3], "z":self.z}


        self.boxsize = self.N * self.dx


    def load_patch(self, origin, N):

        # Retrieve the number of bytes in each slice
        area = self.area
        size = self.size

        fname = self.ic_fname()
        
        # Retrieve the number of points in each slice
        (n1, n2, n3) = self.n

        # Store the data in a patch
        patch = np.zeros(shape=(N, N, N), dtype=np.float32)

        z0 = origin[2]
        z = z0
        # print('z0 = {0}'.format(z0))
        # print('size = {0}'.format(size))
        # print('area = {0}'.format(area))
        # print('origin = {0}'.format(origin))
        # print('fname = {0}'.format(fname))
        with open(fname, "rb") as f:
            # Seek past the header block to the initial point, then
            # move along another z times to the z-origin
            f.seek((size + (area + 8)*z), 0)

            # The z-axis changes the slowest
            for iz in range(N):
                # Pick out the plane, and reshape to (x, y)
                data = np.fromfile(f, dtype=np.float32,
                                   count=(n1 * n2)).reshape((n1, n2))

                # Originally I had this written as patch[:, :, iz] =
                # ... but in order to get the same results as seren3 I
                # have to use patch[iz, :, :] = ... -- I'm not
                # entirely convinced by this, because the way I see it
                # we are moving along z, not x, but I'll leave it for
                # now to be consistent
                patch[iz, :, :] = data[origin[0]:origin[0]+N,
                                       origin[1]:origin[1]+N]

                # Skip along to the next plane
                z += 1
                f.seek((size + (area + 8)*z), 0)

        return patch

    
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
                  cosmo['omega_M_0'], cosmo['omega_lambda_0'],
                  cosmo['h'] * 100.0], dtype=np.float32).tofile(f)
        header_bytes.tofile(f)


    def write_field(self, data, field_name, out_dir=None, **kwargs):
        """Function for writing the data as a Fortran unformatted
        binary. Writes a header to mimic grafic IC files, as well as
        appending and prepending each 'record' (i.e. each [ix, :, :])
        with the number of bytes in that slice, as done by Fortran
        (when compiled with GCC).

        :param data: (numpy.ndarray) data should be in np.float32 format and 
                     have (n1, n2, n3) elements
        :param field_name: (str) field to write as, will be prepended with 'ic_' 
        :param out_dir: (str) directory to write the field to
        :returns: None
        :rtype:

        """
        import os
        # print('data.shape = {0}'.format(data.shape))
        # print('type(data.shape) = {0}'.format(type(data.shape)))
        # print('self.n = {0}'.format(self.n))
        # print('type(self.n) = {0}'.format(type(self.n)))
        # assert(data.shape == self.n)

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

                for ix in range(n1):
                    # Record header
                    area.tofile(f)
                    # Record
                    tmp = data[ix, :, :].flatten().astype(np.float32)
                    tmp.tofile(f)
                    # Record footer
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


def load_snapshot(path, level, field):
    """Convenience function for loading grafic ICs.

    :param path: (str) path to directory containing the level_... directories
    :param level: (int) level to load
    :param field: (str) field to load, i.e. 'deltab' to load file 'ic_deltab'
    :returns: the grafic IC in a Snapshot class
    :rtype: Snapshot

    """
    return Snapshot(path, level, field)


def derive_vel(path, level, species):
    """Calculates the magnitude of the velocity fields for the specified
    species. Then writes the ic_velb and ic_velc files.

    :param path: (str) path to directory containing the level_... directories
    :param level: (int) level to load
    :param species: (str) should be either 'b' or 'c'
    :returns: None
    :rtype:

    """
    # Check that species is correct
    assert(species == 'b' or species == 'c'), "Species should be either 'b' or 'c'"
    
    # First create the fields
    fields = ['vel{0}{1}'.format(species, c) for c in ['x', 'y', 'z']]

    # Load up each snapshot
    s = [load_snapshot(path, level, field) for field in fields]

    # Get the sizes of the IC files
    N = s[0].N
    assert(s[0].N == s[1].N == s[2].N), 'velx, vely and velz files different sizes.'
    
    # Now extract the values from each snapshot
    vx = s[0].load_patch([0, 0, 0], N)
    vy = s[1].load_patch([0, 0, 0], N)
    vz = s[2].load_patch([0, 0, 0], N)

    # Calculate velocity
    v = np.sqrt(vx**2 + vy**2 + vz**2)

    # These files may be quite large
    del vx
    del vy
    del vz

    # Now write this data to a field
    out_dir = s[0].level_dir
    s[0].write_field(v, 'vel{0}'.format(species), out_dir=out_dir)


def derive_vbc(path, level):
    """Calculates the vbc field from the velb and velc fields. If the velb
    and velc fields don't already exist, this function creates
    them. Then writes the ic_vbc file.

    :param path: (str) path to directory containing the level_... directories
    :param level: (int) level to load
    :returns: None
    :rtype:

    """
    import os
    
    # Check that the velb and velc fields have been written
    if not os.path.isfile(path+'level_{0:03d}/ic_velb'.format(level)):
        derive_vel(path, level, species='b')
    if not os.path.isfile(path+'level_{0:03d}/ic_velc'.format(level)):
        derive_vel(path, level, species='c')

    # if they have, go ahead and load them up
    vb = load_snapshot(path, level, 'velb')
    vc = load_snapshot(path, level, 'velc')
    
    # Check they are the same size
    assert(vb.N == vc.N), 'ic_velb and ic_velc different sizes.'
    N = vb.N
    
    # if so, load up the actual data
    vb_data = vb.load_patch([0, 0, 0], N)
    vc_data = vc.load_patch([0, 0, 0], N)

    # Now calculate v_bc
    vbc = vb_data - vc_data

    # And finally write the field
    out_dir = vb.level_dir
    vb.write_field(vbc, 'vbc', out_dir=out_dir)
