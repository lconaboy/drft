! Collection of subroutines designed to calculate CIC interpolated
! quanitites. For use with f2py.
!
! L Conaboy, April 2020

subroutine delc(path, omega_b)
  ! --------------------
  !
  ! Subroutine for computing CDM overdensity (\delta_c) from CIC
  ! interpolated particle positions.
  !
  ! Input
  ! -----
  !
  ! path - (char) patch IC level_xxx directory, where the fields will
  ! be written back out to
  !
  ! omega_b - (real) baryon density parameter
  !
  ! --------------------
  
  implicit none
  
  character (len=200), intent(in) :: path
  real, intent(in) :: omega_b
  
  integer :: hblk = 44
  integer :: mblk
  integer :: i, j, k, n1, n2, n3
  integer :: ip, jp, kp
  integer, parameter :: f=50
  real :: mp = 0.0
  real :: dxx, dyy, dzz, tx, ty, tz
  real :: omega_c, boxvol
  real :: dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  real :: rho_cr = 2.775e11  ! h^2 Msol/Mpc^3
  real, allocatable, dimension(:, :, :) :: dx, dy, dz, del_c

  write(6, *) 'Interpolating CDM overdensity'
  
  write(6, *) '---- reading ic_poscx'
  ! Read the x data
  open(f, file=trim(path)//'ic_poscx', form='unformatted')
  rewind f
  ! Read header and print out cosmological parameters, for quick checking
  read(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  write(6, *) '-------- n1, n2, n3', n1, n2, n3
  write(6, *) '-------- astart', astart, ' zstart', 1.0/astart - 1.0
  write(6, *) '-------- omega_m', omega_m
  write(6, *) '-------- omega_l', omega_l
  write(6, *) '-------- H_0', h0
  call flush(6)

  ! We also have to read in omega_b
  ! open(f+1, file='omega_b.dat', form='formatted')
  ! read(f+1, *) omega_b
  ! close(f+1)
  write(6, *) '-------- omega_b', omega_b

  ! Now calculate the particle mass
  rho_cr = rho_cr * (h0/100.)**2
  omega_c = omega_m - omega_b
  boxvol = real(n1*n2*n3) * (dxini**3.)

  ! Calculate mp
  mp = omega_c * rho_cr * boxvol / real(n1*n2*n3)  ! M_sol


  write(6, *) '---- calculated parameters'
  write(6, *) '-------- mp (M_sol)', mp
  write(6, *) '-------- boxvol (Mpc^3)', boxvol
  
  
  ! Now we have n1, n2 and n3, we can allocate arrays
  allocate(dx(n1, n2, n3), dy(n1, n2, n3), dz(n1, n2, n3))
  
  ! and read the data
  do k=1,n3
     read(f) ((dx(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)
  
  ! Next read the y data
  ! f = f + 1
  write(6, *) '---- reading ic_poscy'
  call flush(6)
  open(f, file=trim(path)//'ic_poscy', form='unformatted')
  rewind f
  ! Skip header
  read(f)
  ! Read the data
  do k=1,n3
     read(f) ((dy(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)

  write(6, *) '---- reading ic_poscz'
  call flush(6)
  ! f = f + 1
  open(f, file=trim(path)//'ic_poscz', form='unformatted')
  rewind f
  ! Skip header
  read(f)
  ! Read the data
  do k=1,n3
     read(f) ((dz(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)

  allocate(del_c(n1, n2, n3))
  del_c = 0.0 
  
  ! Interpolate the particles onto a grid
  write(6, *) '---- interpolating'
  call flush(6)
  do i=1,n1
     do j=1,n2
        do k=1,n3

           ! Allow for periodic boundary conditions
           ip = mod(i, n1) + 1
           jp = mod(j, n2) + 1
           kp = mod(k, n3) + 1

           ! Displacement of the particle in the current parent cell
           ! at (i, j, k), converting from cell width of dxini to a
           ! cell width of one
           dxx = dx(i, j, k) / dxini
           dyy = dy(i, j, k) / dxini
           dzz = dz(i, j, k) / dxini

           ! Convenience variables
           tx = 1.0 - dxx
           ty = 1.0 - dyy
           tz = 1.0 - dzz

           ! Interpolate using cloud-in-cell
           del_c(i, j, k) = del_c(i, j, k) + mp*tx*ty*tz
           del_c(ip, j, k) = del_c(ip, j, k) + mp*dxx*ty*tz
           del_c(i, jp, k) = del_c(i, jp, k) + mp*tx*dyy*tz
           del_c(i, j, kp) = del_c(i, j, kp) + mp*tx*ty*dzz
           del_c(ip, jp, k) = del_c(ip, jp, k) + mp*dxx*dyy*tz
           del_c(ip, j, kp) = del_c(ip, j, kp) + mp*dxx*ty*dzz
           del_c(i, jp, kp) = del_c(i, jp, kp) + mp*tx*dyy*dzz
           del_c(ip, jp, kp) = del_c(ip, jp, kp) + mp*dxx*dyy*dzz
           
        end do
     end do
  end do

  write(6, *) '---- done interpolating'
  call flush(6)

  ! Convert del_c to proper density
  del_c = del_c / dxini**3  ! M_sol Mpc^-3

  ! Convert to overdensity
  del_c = del_c/(rho_cr * omega_c) - 1.
  
  deallocate(dx)
  deallocate(dy)
  deallocate(dz)
  
  ! Write the data, including an artefact of grafic2 (perhaps to do
  ! with f77?) whereby an int specifying the number of bytes in the
  ! proceeding block is written before and after that block. For the
  ! header, this is always 44 (11 items multiplyed by 4 bytes). For
  ! the main block, this depends on the number of particles in each
  ! slice, although since they are reals, it is always a multiple of
  ! 4.
  mblk = int(n1 * n2 * 4)
  
  write(6, *) '---- writing deltac'
  call flush(6)
  open(f, file=trim(path)//'ic_deltac', form='unformatted', status='new', access='stream')
  rewind f
  ! Write the header
  write(f) hblk
  write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  write(f) hblk
  ! Write the main block
  do k=1,n3
     write(f) mblk
     write(f) ((del_c(i, j, k), i=1,n1), j=1,n2)
     write(f) mblk
  end do
  close(f)

  deallocate(del_c)

end subroutine delc

subroutine velc(path)
  ! --------------------
  !
  ! Subroutine for computing CIC interpolated gridded velocities.
  !
  ! Input
  ! -----
  !
  ! path - (char) patch IC level_xxx directory, where the fields will
  ! be written back out to
  !
  ! --------------------

  
  character (len=200), intent(in) :: path
  
  integer :: hblk = 44
  integer :: mblk
  integer :: i, j, k, n1, n2, n3
  integer :: ip, jp, kp
  integer, parameter :: f=50
  real :: mp = 1.0
  real :: dxx, dyy, dzz, tx, ty, tz
  real :: dxini, x1off, x2off, x3off, astart, omegam, omegal, omega_b, h0
  real, parameter :: rhoc = 2.775e11  ! h^2 Msol/Mpc^3
  real, allocatable, dimension(:, :, :) :: dx, dy, dz, vcx, vcy, vcz, vcg, vc, vbx, vby, vbz, vb, vbc

  write(6, *) 'Interpolating CDM velocities'
  
  write(6, *) '---- reading ic_poscx'
  ! Read the x data
  open(f, file=trim(path)//'ic_poscx', form='unformatted')
  open(f+1, file=trim(path)//'ic_velcx', form='unformatted')
  rewind f
  ! Read header and print out cosmological parameters, for quick checking
  read(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omegam, omegal, h0
  read(f+1)
  write(6, *) '-------- n1, n2, n3', n1, n2, n3
  write(6, *) '-------- astart', astart, ' zstart', 1.0/astart - 1.0
  write(6, *) '-------- omega_m', omegam
  write(6, *) '-------- omega_l', omegal
  write(6, *) '-------- H_0', h0
  call flush(6)

  ! Now we have n1, n2 and n3, we can allocate arrays
  allocate(dx(n1, n2, n3), dy(n1, n2, n3), dz(n1, n2, n3))
  allocate(vcx(n1, n2, n3))
  allocate(vc(n1, n2, n3))

  ! Initialise vcg to zero
  vc = 0.0
  
  ! and read the data
  do k=1,n3
     read(f) ((dx(i, j, k), i=1,n1), j=1,n2)
     read(f+1) ((vcx(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)
  close(f+1)

  ! Store vcx in vcg so we can deallocate
  do i=1,n1
     do j=1,n2
        do k=1,n3
           vc(i, j, k) = vc(i, j, k) + vcx(i, j, k)**2
        end do
     end do
  end do
  deallocate(vcx)
  
  ! Next read the y data
  allocate(vcy(n1, n2, n3))
  ! f = f + 1
  write(6, *) '---- reading ic_poscy'
  call flush(6)
  open(f, file=trim(path)//'ic_poscy', form='unformatted')
  open(f+1, file=trim(path)//'ic_velcy', form='unformatted')
  rewind f
  ! Skip header
  read(f)
  read(f+1)
  ! Read the data
  do k=1,n3
     read(f) ((dy(i, j, k), i=1,n1), j=1,n2)
     read(f+1) ((vcy(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)
  close(f+1)

  ! Store vcy in vcg so we can deallocate
  do i=1,n1
     do j=1,n2
        do k=1,n3
           vc(i, j, k) = vc(i, j, k) + vcy(i, j, k)**2
        end do
     end do
  end do
  deallocate(vcy)
  
  ! Finally, the z data
  allocate(vcz(n1, n2, n3))
  write(6, *) '---- reading ic_poscz'
  call flush(6)
  ! f = f + 1
  open(f, file=trim(path)//'ic_poscz', form='unformatted')
  open(f+1, file=trim(path)//'ic_velcz', form='unformatted')
  rewind f
  ! Skip header
  read(f)
  read(f+1)
  ! Read the data
  do k=1,n3
     read(f) ((dz(i, j, k), i=1,n1), j=1,n2)
     read(f+1) ((vcz(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)
  close(f+1)

  ! Store vcz in vcg so we can deallocate
  do i=1,n1
     do j=1,n2
        do k=1,n3
           vc(i, j, k) = sqrt(vc(i, j, k) + vcz(i, j, k)**2)
        end do
     end do
  end do
  deallocate(vcz)

  allocate(vcg(n1, n2, n3))
  vcg = 0.0
  
  ! Interpolate the velocities onto a grid
  write(6, *) '---- interpolating'
  call flush(6)
  do i=1,n1
     do j=1,n2
        do k=1,n3

           ! Allow for periodic boundary conditions
           ip = mod(i, n1) + 1
           jp = mod(j, n2) + 1
           kp = mod(k, n3) + 1

           ! Displacement of the particle in the current parent cell
           ! at (i, j, k), converting from cell width of dxini to a
           ! cell width of one
           dxx = dx(i, j, k) / dxini
           dyy = dy(i, j, k) / dxini
           dzz = dz(i, j, k) / dxini

           ! Convenience variables
           tx = 1.0 - dxx
           ty = 1.0 - dyy
           tz = 1.0 - dzz

           vcg(i, j, k) = vcg(i, j, k) + vc(i, j, k)*tx*ty*tz
           vcg(ip, j, k) = vcg(ip, j, k) + vc(ip, j, k)*dxx*ty*tz
           vcg(i, jp, k) = vcg(i, jp, k) + vc(i, jp, k)*tx*dyy*tz
           vcg(i, j, kp) = vcg(i, j, kp) + vc(i, j, kp)*tx*ty*dzz
           vcg(ip, jp, k) = vcg(ip, jp, k) + vc(ip, jp, k)*dxx*dyy*tz
           vcg(ip, j, kp) = vcg(ip, j, kp) + vc(ip, j, kp)*dxx*ty*dzz
           vcg(i, jp, kp) = vcg(i, jp, kp) + vc(i, jp, kp)*tx*dyy*dzz
           vcg(ip, jp, kp) = vcg(ip, jp, kp) + vc(ip, jp, kp)*dxx*dyy*dzz

        end do
     end do
  end do

  write(6, *) '---- done interpolating'
  call flush(6)

  deallocate(dx)
  deallocate(dy)
  deallocate(dz)
  deallocate(vc)
  
  ! Write the data, including an artefact of grafic2 (perhaps to do
  ! with f77?) whereby an int specifying the number of bytes in the
  ! proceeding block is written before and after that block. For the
  ! header, this is always 44 (11 items multiplyed by 4 bytes). For
  ! the main block, this depends on the number of particles in each
  ! slice, although since they are reals, it is always a multiple of
  ! 4.
  mblk = int(n1 * n2 * 4)
  
  write(6, *) '---- writing velcg'
  call flush(6)
  open(f, file=trim(path)//'ic_velcg', form='unformatted', status='new', access='stream')
  rewind f
  ! Write the header
  write(f) hblk
  write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omegam, omegal, h0
  write(f) hblk
  ! Write the main block
  do k=1,n3
     write(f) mblk
     write(f) ((vcg(i, j, k), i=1,n1), j=1,n2)
     write(f) mblk
  end do
  close(f)

  ! Now read in ic_velb* and calculate ic_vbc
  allocate(vbx(n1, n2, n3))
  allocate(vb(n1, n2, n3))
  vb = 0.0
  
  write(6, *) '---- reading ic_velbx'
  call flush(6)
  open(f, file=trim(path)//'ic_velbx', form='unformatted')
  rewind f
  ! Skip header
  read(f)
  ! Read the data
  do k=1,n3
     read(f) ((vbx(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)

  do i=1,n1
     do j=1,n2
        do k=1,n3
           vb(i, j, k) = vb(i, j, k) + vbx(i, j, k) ** 2
        end do
     end do
  end do

  deallocate(vbx)

  allocate(vby(n1, n2, n3))
  write(6, *) '---- reading ic_velby'
  call flush(6)
  open(f, file=trim(path)//'ic_velby', form='unformatted')
  rewind f
  ! Skip header
  read(f)
  ! Read the data
  do k=1,n3
     read(f) ((vby(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)
  
  do i=1,n1
     do j=1,n2
        do k=1,n3
           vb(i, j, k) = vb(i, j, k) + vby(i, j, k) ** 2
        end do
     end do
  end do
  deallocate(vby)

  allocate(vbz(n1, n2, n3))
  write(6, *) '---- reading ic_velbz'
  call flush(6)
  open(f, file=trim(path)//'ic_velbz', form='unformatted')
  rewind f
  ! Skip header
  read(f)
  ! Read the data
  do k=1,n3
     read(f) ((vbz(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)
  
  do i=1,n1
     do j=1,n2
        do k=1,n3
           vb(i, j, k) = sqrt(vb(i, j, k) + vbz(i, j, k) ** 2)
        end do
     end do
  end do
  deallocate(vbz)

  write(6, *) '---- writing velb'
  call flush(6)
  open(f, file=trim(path)//'ic_velb', form='unformatted', status='new', access='stream')
  rewind f
  ! Write the header
  write(f) hblk
  write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omegam, omegal, h0
  write(f) hblk
  ! Write the main block
  do k=1,n3
     write(f) mblk
     write(f) ((vb(i, j, k), i=1,n1), j=1,n2)
     write(f) mblk
  end do

  close(f)

  ! Now calculate v_bc = v_b - v_c
  allocate(vbc(n1, n2, n3))
  vbc = 0.0

  do i=1,n1
     do j=1,n2
        do k=1,n3
           vbc(i, j, k) = vb(i, j, k) - vcg(i, j, k)
        end do
     end do
  end do

  deallocate(vb)
  deallocate(vcg)

  write(6, *) '---- writing vbc'
  call flush(6)
  open(f, file=trim(path)//'ic_vbc', form='unformatted', status='new', access='stream')
  rewind f
  ! Write the header
  write(f) hblk
  write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omegam, omegal, h0
  write(f) hblk
  ! Write the main block
  do k=1,n3
     write(f) mblk
     write(f) ((vbc(i, j, k), i=1,n1), j=1,n2)
     write(f) mblk
  end do
  close(f)

end subroutine velc