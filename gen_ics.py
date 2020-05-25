def main(path, level, levelmin, cur_dir=False, um=None):
    """Function for generating initial conditions files using the
    ic_deltab, ic_deltac files and various linear approximations,
    which are detailed in linear.py. Writes the IC files and a
    subdirectory of path, called path/level_xxx/ics_linear/level_xxx.

    :param path: 
        (str) path to original IC files (they won't be modified)
    :param level: 
        (int) level of the IC files
    :param levelmin:
        (int) minimum overall level of the IC files, to determine 
        whether periodic CIC interpolation is used
    :param cur_dir:
        (bool) write the IC files to the current dir (i.e. path/level_xxx)
    :param um: 
        (str or None) what to do with unmodified ic_deltab (and 
        possibly ic_refmap and ic_pvar_00001): 'cp' to copy them 
        into the new linear directory, 'sl' to symlink them or None 
        to do nothing
    :returns: 
    :rtype:

    """
    import os
    import utils as vbc_utils
    from linear import generate
    import grafic_tools as grafic
    
    xyz = 'xyz'

    # Periodic CIC?
    per = (level == levelmin)
    
    if cur_dir:
        write_dir = path
        write_path = path + 'level_{0:03d}/'.format(level)
        um = None  # don't do anything with unmodified files
    else:
        write_dir = path + 'level_{0:03d}/ics_linear/'.format(level)
        write_path = path + 'level_{0:03d}/ics_linear/level_{0:03d}/'.format(level)

    print('Generating initial conditon files for level {0} using linear approximations'.format(level))

    # If not reusing the directory, make a new one
    if not cur_dir:
        # Make directories to store linear ICs
        if not os.path.isdir(write_dir):
            os.mkdir(write_dir)
        # Make level directory to store linear ICs
        if not os.path.isdir(write_path):
            os.mkdir(write_path)
        # Check the ics_linear directory is empty
        if len(os.listdir(write_path)) > 0:
            raise Exception(write_path+' not empty')

    # Generate deltac field if it doesn't exist
    if not os.path.isfile(write_path+'ic_deltac'):
        grafic.derive_deltac(path, level, per)

    # Generate the velocity fields for the baryons
    print('---- deriving baryon velocities')

    ic = grafic.load_snapshot(path, level, 'deltab')
    f = generate(ic=ic, qty='vel', params=None, approx=True, real=True)
    # and write them out
    for i in range(3):
        ic.write_field(f[i], 'velb'+xyz[i], out_dir=write_path)

    # Now generate the velocity fields for the CDM
    print('---- deriving dark matter velocities')

    ic = grafic.load_snapshot(path, level, 'deltac')
    f = generate(ic=ic, qty='vel', params=None, approx=True, real=True)
    # and write them out
    for i in range(3):
        ic.write_field(f[i], 'velc'+xyz[i], out_dir=write_path)

    # Npw generate the displacement fields for the CDM
    print('---- deriving dark matter displacements')

    ic = grafic.load_snapshot(path, level, 'deltac')
    f = generate(ic=ic, qty='pos', params=None, approx=True, real=True)
    # and write them out
    for i in range(3):
        ic.write_field(f[i], 'posc'+xyz[i], out_dir=write_path)

    # Finally, the linear ICs need the original deltab, deltac, refmap
    # and pvar files. Make this optional because they may be large
    # files and some might prefer symlinking
    if um == 'sl':
        print('---- symlinking unmodified files')

        os.symlink(path+'level_{0:03d}/ic_deltab'.format(level),
                     write_path+'ic_deltab')

        os.symlink(path+'level_{0:03d}/ic_deltac'.format(level),
                   write_path+'ic_deltac')
        
        # Might not be zoom
        try:
            os.symlink(path+'level_{0:03d}/ic_refmap'.format(level),
                         write_path+'ic_refmap')
        except:
            print('-------- not symlinking refmap')

        # Might not be have passive variables
        try:
            os.symlink(path+'level_{0:03d}/ic_pvar_00001'.format(level),
                         write_path+'ic_pvar_00001')
        except:
            print('-------- not symlinking pvar')

    elif um == 'cp':
        import shutil

        print('---- copying unmodified files')

        shutil.copy2(path+'level_{0:03d}/ic_deltab'.format(level),
                     write_path+'ic_deltab')

        shutil.copy2(path+'level_{0:03d}/ic_deltac'.format(level),
                     write_path+'ic_deltac')
        
        # Might not be zoom
        try:
            shutil.copy2(path+'level_{0:03d}/ic_refmap'.format(level),
                         write_path+'ic_refmap')
        except:
            print('-------- not copying refmap')

        # Might not have passive variables, try up to four just in
        # case. Could do this more efficiently with a glob.
        for ii in range(1, 4):
            try:
                shutil.copy2(path+'level_{0:03d}/ic_pvar_{1:05d}'.format(level, ii),
                             write_path+'ic_pvar_{0:05d}'.format(ii))
            except:
                print('-------- not copying pvar_{0:05d}'.format(ii))

    else:
        print('---- not doing anything with unmodified files (e.g. deltab)!')

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 4:
        print('Usage: python gen_ics.py <path/to/ics/> <current level> <min level> [<write to current dir?> <unmodified files?>]')
        sys.exit()
    
    path = sys.argv[1]
    level = int(sys.argv[2])
    levelmin = int(sys.argv[3])
    cur_dir = False
    um = None
    
    # Optional
    if len(sys.argv) > 4:
        cur_dir = bool(int(sys.argv[4]))
    if len(sys.argv) > 5:
        um = sys.argv[5]

    print('Working on', path)
    print('Working on level', level,'with levelmin', levelmin)
    print('Working in modes: cur_dir', cur_dir, 'um', um)
        
    main(path, level, levelmin, cur_dir, um)
