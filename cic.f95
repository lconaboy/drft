! Collection of subroutines designed to calculate CIC interpolated
! quanitites. For use with f2py.
!
! Example
!
!     from cic import gen_delc
!
!     path = /path/to/ics/
!     l = int(level)
!     omega_b = float(omega_b)
!
!     gen_delc(path, l, omega_b)
!
! Notes
!
!    - displacements in grafic ic_posc files are in units of comoving
!      Mpc/h, while the cell widths (dxini in the grafic headers) are
!      in units of comoving Mpc
!
! L Conaboy, April 2020

subroutine gen_delc(path, l, omega_b)
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
  integer, intent(in) :: l
  real, intent(in) :: omega_b
  
  integer :: hblk = 44
  integer :: mblk
  integer :: i, j, k, n1, n2, n3
  integer :: ip, jp, kp
  integer :: ii, jj, kk, iip, jjp, kkp
  integer, parameter :: f=50
  real :: mp = 0.0
  real :: dxx, dyy, dzz, tx, ty, tz
  real :: omega_c, boxvol
  real :: dx0, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  real :: rho_cr = 2.775e11  ! h^2 Msol/Mpc^3
  real :: cur_max
  real, allocatable, dimension(:, :, :) :: dx, dy, dz, del_c

  write(6, *) 'Interpolating CDM overdensity'

  dx0 = 0.5 ** real(l)
  
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

  ! allocate(del_c(n1+1, n2+1, n3+1))
  allocate(del_c(n1+2, n2+2, n3+2))
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

           ii = i + 1
           jj = j + 1
           kk = k + 1

           iip = ii + 1
           jjp = jj + 1
           kkp = kk + 1
           
           ! Use non-periodic boundary conditions, we have allocated
           ! del_c to be one bigger in each dimension
           ! ip = i + 1
           ! jp = j + 1
           ! kp = k + 1
           
           ! Displacement of the particle in the current parent cell
           ! at (i, j, k), converting from cell width of dxini to a
           ! cell width of one
           dxx = dx(i, j, k) / dxini! / (h0 / 100.) !* dx0 / dxini
           dyy = dy(i, j, k) / dxini! / (h0 / 100.) !* dx0 / dxini
           dzz = dz(i, j, k) / dxini! / (h0 / 100.) !* dx0 / dxini

           ! Convenience variables
           tx = 1.0 - dxx
           ty = 1.0 - dyy
           tz = 1.0 - dzz

           ! Interpolate using cloud-in-cell
           ! del_c(i, j, k) = del_c(i, j, k) + mp*tx*ty*tz
           ! del_c(ip, j, k) = del_c(ip, j, k) + mp*dxx*ty*tz
           ! del_c(i, jp, k) = del_c(i, jp, k) + mp*tx*dyy*tz
           ! del_c(i, j, kp) = del_c(i, j, kp) + mp*tx*ty*dzz
           ! del_c(ip, jp, k) = del_c(ip, jp, k) + mp*dxx*dyy*tz
           ! del_c(ip, j, kp) = del_c(ip, j, kp) + mp*dxx*ty*dzz
           ! del_c(i, jp, kp) = del_c(i, jp, kp) + mp*tx*dyy*dzz
           ! del_c(ip, jp, kp) = del_c(ip, jp, kp) + mp*dxx*dyy*dzz

           del_c(ii, jj, kk) = del_c(ii, jj, kk) + mp*tx*ty*tz
           del_c(iip, jj, kk) = del_c(iip, jj, kk) + mp*dxx*ty*tz
           del_c(ii, jjp, kk) = del_c(ii, jjp, kk) + mp*tx*dyy*tz
           del_c(ii, jj, kkp) = del_c(ii, jj, kkp) + mp*tx*ty*dzz
           del_c(iip, jjp, kk) = del_c(iip, jjp, kk) + mp*dxx*dyy*tz
           del_c(iip, jj, kkp) = del_c(iip, jj, kkp) + mp*dxx*ty*dzz
           del_c(ii, jjp, kkp) = del_c(ii, jjp, kkp) + mp*tx*dyy*dzz
           del_c(iip, jjp, kkp) = del_c(iip, jjp, kkp) + mp*dxx*dyy*dzz

        end do
     end do
  end do

  write(6, *) '---- done interpolating'
  call flush(6)

  ! Convert del_c to proper density
  del_c = del_c / dxini**3  ! M_sol Mpc^-3

  ! Convert to overdensity
  del_c = del_c/(rho_cr * omega_c) - 1.

  ! Because of how we treated the periodic boundary conds...
  del_c = del_c(2:n1+1, 2:n2+1, 2:n3+1)
  del_c(:, :, 1) = 0.0
  del_c(:, 1, :) = 0.0
  del_c(1, :, :) = 0.0
  del_c(n1, :, :) = 0.0
  del_c(:, n2, :) = 0.0
  del_c(:, :, n3) = 0.0
  
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
  open(f, file=trim(path)//'ic_deltac', form='unformatted') !, access='stream')
  rewind f
  ! Write the header
  ! write(f) hblk
  write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  ! write(f) hblk
  ! Write the main block
  do k=1,n3
     ! write(f) mblk
     write(f) ((del_c(i, j, k), i=1,n1), j=1,n2)
     ! write(f) mblk
  end do
  close(f)

  deallocate(del_c)

  write(6, *) '---- done'
end subroutine gen_delc


subroutine gen_velcg(path, l)
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
  character (len=3) :: xyz

  integer, intent(in) :: l
  integer :: hblk = 44
  integer :: mblk
  integer :: i, j, k, n1, n2, n3
  integer :: ip, jp, kp, ll
  integer :: ii, jj, kk, iip, jjp, kkp
  integer, parameter :: f=50
  real :: mp = 1.0
  real :: dxx, dyy, dzz, tx, ty, tz
  real :: dx0, dxini, x1off, x2off, x3off, astart, omegam, omegal, omega_b, h0
  real :: cur_max, cur_min
  real, parameter :: rhoc = 2.775e11  ! h^2 Msol/Mpc^3
  real, allocatable, dimension(:, :, :) :: dx, dy, dz, vc, vcg

  xyz = 'xyz'

  dx0 = 0.5 ** real(l)  ! Cell spacing for level i
  
  write(6, *) 'Interpolating CDM velocities'

  ! Read the x data
  write(6, *) '---- reading parameters'
  open(f, file=trim(path)//'ic_poscx', form='unformatted')
  rewind f
  ! Read header and print out cosmological parameters, for quick checking
  read(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omegam, omegal, h0
  write(6, *) '-------- n1, n2, n3', n1, n2, n3
  write(6, *) '-------- astart', astart, ' zstart', 1.0/astart - 1.0
  write(6, *) '-------- omega_m', omegam
  write(6, *) '-------- omega_l', omegal
  write(6, *) '-------- H_0', h0
  call flush(6)

  ! Load the offsets
  allocate(dx(n1, n2, n3), dy(n1, n2, n3), dz(n1, n2, n3))
  write(6, *) '---- reading posc'
  open(f+1, file=trim(path)//'ic_poscy', form='unformatted')
  open(f+2, file=trim(path)//'ic_poscz', form='unformatted')
  ! Skip headers
  read(f+1)
  read(f+2)

  ! Read in data
  do k = 1, n3
     read(f) ((dx(i, j, k), i=1,n1), j=1,n2)
     read(f+1) ((dy(i, j, k), i=1,n1), j=1,n2)
     read(f+2) ((dz(i, j, k), i=1,n1), j=1,n2)
  end do

  ! Close files
  close(f)
  close(f+1)
  close(f+2)

  cur_max = 0.0
  cur_min = 0.0
  
  ! Loop over dims and grid each velocity
  do ll = 1, 3
     allocate(vc(n1, n2, n3), vcg(n1+2, n2+2, n3+2))

     ! Read in the the CDM velocity
     write(6, *) '---- reading velc'//xyz(ll:ll)
     open(f, file=trim(path)//'ic_velc'//xyz(ll:ll), form='unformatted')
     read(f)
     do k = 1, n3
        read(f) ((vc(i, j, k), i=1,n1), j=1,n2)
     end do
     close(f)

     ! Now interpolate
     vcg = 0.0
     
     write(6, *) '---- interpolating'
     call flush(6)
     do i=1,n1
        do j=1,n2
           do k=1,n3

              ! Allow for periodic boundary conditions
              ip = mod(i, n1) + 1
              jp = mod(j, n2) + 1
              kp = mod(k, n3) + 1

              ii = i + 1
              jj = j + 1
              kk = k + 1

              iip = ii + 1
              jjp = jj + 1
              kkp = kk + 1
           
              ! Interpolate using cloud-in-cell
              ! del_c(i, j, k) = del_c(i, j, k) + mp*tx*ty*tz
              ! del_c(ip, j, k) = del_c(ip, j, k) + mp*dxx*ty*tz
              ! del_c(i, jp, k) = del_c(i, jp, k) + mp*tx*dyy*tz
              ! del_c(i, j, kp) = del_c(i, j, kp) + mp*tx*ty*dzz
              ! del_c(ip, jp, k) = del_c(ip, jp, k) + mp*dxx*dyy*tz
              ! del_c(ip, j, kp) = del_c(ip, j, kp) + mp*dxx*ty*dzz
              ! del_c(i, jp, kp) = del_c(i, jp, kp) + mp*tx*dyy*dzz
              ! del_c(ip, jp, kp) = del_c(ip, jp, kp) + mp*dxx*dyy*dzz

              ! Displacement of the particle in the current parent
              ! cell at (i, j, k), converting from cell width of dxini
              ! to a cell width of one
              dxx = dx(i, j, k) / dxini! / (h0 / 100.)
              dyy = dy(i, j, k) / dxini! / (h0 / 100.)
              dzz = dz(i, j, k) / dxini! / (h0 / 100.)

              ! Conenience variables
              tx = 1.0 - dxx
              ty = 1.0 - dyy
              tz = 1.0 - dzz

              ! write(6, *), 'tx', tx
              ! write(6, *), 'dx(i,j,k)', dx(i, j, k)
              ! write(6, *), 'dxx', dxx

              if (dxx .gt. cur_max) then
                 cur_max = dxx
              else if (dxx .lt. cur_min) then
                 cur_min = dxx
              end if


              vcg(ii, jj, kk) = vcg(ii, jj, kk) + vc(i, j, k)*tx*ty*tz
              vcg(iip, jj, kk) = vcg(iip, jj, kk) + vc(ip, j, k)*dxx*ty*tz
              vcg(ii, jjp, kk) = vcg(ii, jjp, kk) + vc(i, jp, k)*tx*dyy*tz
              vcg(ii, jj, kkp) = vcg(ii, jj, kkp) + vc(i, j, kp)*tx*ty*dzz
              vcg(iip, jjp, kk) = vcg(iip, jjp, kk) + vc(ip, jp, k)*dxx*dyy*tz
              vcg(iip, jj, kkp) = vcg(iip, jj, kkp) + vc(ip, j, kp)*dxx*ty*dzz
              vcg(ii, jjp, kkp) = vcg(ii, jjp, kkp) + vc(i, jp, kp)*tx*dyy*dzz
              vcg(iip, jjp, kkp) = vcg(iip, jjp, kkp) + vc(ip, jp, kp)*dxx*dyy*dzz
              

              ! vcg(i, j, k) = vcg(i, j, k) + vc(i, j, k)*tx*ty*tz
              ! vcg(ip, j, k) = vcg(ip, j, k) + vc(ip, j, k)*dxx*ty*tz
              ! vcg(i, jp, k) = vcg(i, jp, k) + vc(i, jp, k)*tx*dyy*tz
              ! vcg(i, j, kp) = vcg(i, j, kp) + vc(i, j, kp)*tx*ty*dzz
              ! vcg(ip, jp, k) = vcg(ip, jp, k) + vc(ip, jp, k)*dxx*dyy*tz
              ! vcg(ip, j, kp) = vcg(ip, j, kp) + vc(ip, j, kp)*dxx*ty*dzz
              ! vcg(i, jp, kp) = vcg(i, jp, kp) + vc(i, jp, kp)*tx*dyy*dzz
              ! vcg(ip, jp, kp) = vcg(ip, jp, kp) + vc(ip, jp, kp)*dxx*dyy*dzz

           end do
        end do
     end do


     write(6, *), '-------- max displacement (cell width = 1.0)', cur_max
     write(6, *), '-------- min displacement (cell width = 1.0)', cur_min
                 
     ! Because of how we treated the periodic boundary conditions, we
     ! can either set the edge cells to zero, or use the ungridded
     ! velocity. Both are incorrect, but I think using the ungridded
     ! value is less wrong than setting to zero.
     vcg = vcg(2:n1+1, 2:n2+1, 2:n3+1)
     vcg(:, :, 1) = vc(:, :, 1)
     vcg(:, 1, :) = vc(:, 1, :)
     vcg(1, :, :) = vc(1, :, :)
     vcg(n1, :, :) = vc(n1, :, :)
     vcg(:, n2, :) = vc(:, n2, :)
     vcg(:, :, n3) = vc(:, :, n3)

     
     ! Write out the gridded velocity
     write(6, *) '---- writing velcg'//xyz(ll:ll)
     open(f, file=trim(path)//'ic_velcg'//xyz(ll:ll), form='unformatted')
     write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
     do k = 1, n3
        write(f) ((vcg(i, j, k), i=1,n1), j=1,n2)
     end do
     close(f)

     ! Clean up
     deallocate(vc, vcg)
  end do
  write(6, *) '---- done'
end subroutine gen_velcg


subroutine gen_vbc(path)
  ! --------------------
  !
  ! Subroutine for calculating the value of |_vb_ - _vc_| = vbc, where
  ! _x_ indicates a vector.
  !
  ! Input
  ! -----
  !
  ! path - (char) patch IC level_xxx directory, where the fields will
  ! be written back out to
  !
  ! --------------------

  character (len=200), intent(in) :: path
  character (len=3) :: xyz

  integer :: hblk = 44
  integer :: mblk
  integer :: i, j, k, n1, n2, n3
  integer, parameter :: f=50
  real :: dxini, x1off, x2off, x3off, astart, omegam, omegal, omega_b, h0
  real, allocatable, dimension(:, :, :) :: vbc, vb, vcg
  logical :: per

  xyz = 'xyz'
  per = .false.

  write(6, *) 'Calculating v_bc field'

  ! Read the x data
  write(6, *) '---- reading parameters'
  open(f, file=trim(path)//'ic_poscx', form='unformatted')
  rewind f
  ! Read header and print out cosmological parameters, for quick checking
  read(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omegam, omegal, h0
  write(6, *) '-------- n1, n2, n3', n1, n2, n3
  write(6, *) '-------- astart', astart, ' zstart', 1.0/astart - 1.0
  write(6, *) '-------- omega_m', omegam
  write(6, *) '-------- omega_l', omegal
  write(6, *) '-------- H_0', h0
  call flush(6)

  ! Calculate vbc as the vector difference between vb and vcg
  allocate(vbc(n1, n2, n3))
  vbc = 0.0
  
  do ii = 1, 3
     allocate(vb(n1, n2, n3), vcg(n1, n2, n3))

     ! Read in the the baryon velocity
     write(6, *) '---- reading velb'//xyz(ii:ii)
     open(f, file=trim(path)//'ic_velb'//xyz(ii:ii), form='unformatted')
     read(f)
     do k = 1, n3
        read(f) ((vb(i, j, k), i=1,n1), j=1,n2)
     end do
     close(f)

     ! Read in the gridded CDM velocity
     write(6, *) '---- reading velcg'//xyz(ii:ii)
     open(f, file=trim(path)//'ic_velcg'//xyz(ii:ii), form='unformatted')
     read(f)
     do k = 1, n3
        read(f) ((vcg(i, j, k), i=1,n1), j=1,n2)
     end do
     close(f)

     do k = 1, n3
        do j = 1, n2
           do i = 1, n1
              vbc(i, j, k) = vbc(i, j, k) + (vb(i, j, k) - vcg(i, j, k))**2
           end do
        end do
     end do

     deallocate(vb, vcg)
  end do
  
  ! Write the vbc field
  vbc = sqrt(vbc)

  ! Because of how we handled non-periodic boundaries in gen_velcg, we
  ! don't need to worry about the edge cases here.

  ! if (.not. per) then
  !    vbc(1, :, :) = 0.0
  !    vbc(:, 1, :) = 0.0
  !    vbc(:, :, 1) = 0.0
  !    vbc(:, :, n3) = 0.0
  !    vbc(:, n2, :) = 0.0
  !    vbc(n1, :, :) = 0.0
  ! end if
  
  write(6, *) '---- writing vbc'
  open(f, file=trim(path)//'ic_vbc', form='unformatted')
  write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  do k = 1, n3
     write(f) ((vbc(i, j, k), i=1,n1), j=1,n2)
  end do
  close(f)

  deallocate(vbc)
  write(6, *) '---- done'
end subroutine gen_vbc
