! Grids particle quantities such that they can be differenced with AMR
! quantities in e.g. initial conditions. Follows pm/rho_fine.f90 in
! RAMSES to carry out the gridding.

! TODO
! - read in values to be gridded

program interp
  implicit none

  interface
     subroutine cic_grid(path, field, per)
       logical, intent(in) :: per
       character (len=5), intent(in):: field
       character (len=200), intent(in):: path
     end subroutine cic_grid

     subroutine gen_vbc(path)
       character(len=200), intent(in) :: path
     end subroutine gen_vbc
  end interface

  integer :: i
  logical :: per
  character (len=3) :: xyz
  character (len=4) :: s
  character (len=5) :: field
  character (len=200) :: path

  xyz = 'xyz'
  per = .true.
  s = 'velc'
  path='./'

  do i=1,3
     write(field, '(A4, A1)') s, xyz(i:i)
     call cic_grid(path, field, per)
  end do

  call gen_vbc(path)
  
end program interp

subroutine cic_grid(path, field, per)
  implicit none

  logical, intent(in) :: per
  character (len=5), intent(in) :: field
  character (len=200), intent(in) :: path
  
  integer :: i, j, k, ic, jc, kc, icp, jcp, kcp, n1, n2, n3
  real :: xp, yp, zp, xl, yl, zl, xr, yr, zr
  real :: dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  real, allocatable :: xg(:, :, :), yg(:, :, :), zg(:, :, :), vp(:, :, :)
  real, allocatable :: vg(:, :, :)
  
  write(6, *) '---- reading ic_poscx'
  ! Read the x data
  open(20, file=trim(path)//'ic_poscx', form='unformatted')
  rewind 20
  ! Read header and print out cosmological parameters, for quick checking
  read(20) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  write(6, *) '-------- n1, n2, n3', n1, n2, n3
  write(6, *) '-------- astart', astart, ' zstart', 1.0/astart - 1.0
  write(6, *) '-------- omega_m', omega_m
  write(6, *) '-------- omega_l', omega_l
  write(6, *) '-------- H_0', h0
  call flush(6)

  ! Now we have n1, n2 and n3, we can allocate arrays
  allocate(xg(n1, n2, n3), yg(n1, n2, n3), zg(n1, n2, n3), vp(n1, n2, n3))

  ! If doing non-periodic, add a buffer region where we can dump all
  ! of the particles that go outside
  if (per) then
     allocate(vg(n1, n2, n3))
  else
     allocate(vg(0:n1+1, 0:n2+1, 0:n3+1))
  end if

  ! Initialise
  vg(:, :, :) = 0.0
  
  open(21, file=trim(path)//'ic_poscy', form='unformatted')
  open(22, file=trim(path)//'ic_poscz', form='unformatted')
  ! and the field to be interpolated
  open(23, file=trim(path)//'ic_'//trim(field), form='unformatted')

  ! Skip the headers
  read(21)
  read(22)
  read(23)

  ! Read the data in
  do k=1,n3
     read(20) ((xg(i, j, k), i=1,n1), j=1,n2)
     read(21) ((yg(i, j, k), i=1,n1), j=1,n2)
     read(22) ((zg(i, j, k), i=1,n1), j=1,n2)
     read(23) ((vp(i, j, k), i=1,n1), j=1,n2)
  end do

  close(20)
  close(21)
  close(22)
  write(6, *) '---- read in position data'
  close(23)
  write(6, *) '---- read in '//trim(field)//' data'

  ! Convert displacements from Mpc/h to cell widths
  xg = xg / (dxini * h0 / 100.)
  yg = yg / (dxini * h0 / 100.)
  zg = zg / (dxini * h0 / 100.)

  ! Convert displacements to positions in units of cells
  do k=1,n3
     do j=1,n2
        do i=1,n1
           xg(i, j, k) = xg(i, j, k) + i
           yg(i, j, k) = yg(i, j, k) + j
           zg(i, j, k) = zg(i, j, k) + k
        end do
     end do
  end do

  ! Print out a test slice
  ! open(20, file='x_slice.dat', form='formatted', status='new')
  ! open(21, file='y_slice.dat', form='formatted', status='new')
  ! do i = 1,n1
  !    do j = 1,n2
  !       write(20, *) xg(i, j, 1)
  !       write(21, *) yg(i, j, 1)
  !    end do
  ! end do
  ! close(20)
  ! close(21)

  ! Do CIC
  do k=1,n3
     do j=1,n2
        do i=1,n1
           xp = xg(i, j, k)
           yp = yg(i, j, k)
           zp = zg(i, j, k)

           ! Calculate weights
           xr = xp + 0.5
           ic = int(xr)
           xr = xr - ic
           xl = 1.0 - xr
           ! ic = ic - 1

           !print*, xl, xr
           !if (i .ne. ic) print*, i, ic

           yr = yp + 0.5
           jc = int(yr)
           yr = yr - jc
           yl = 1.0 - yr
           ! jc = jc - 1

           zr = zp + 0.5
           kc = int(zr)
           zr = zr - kc
           zl = 1.0 - zr
           ! kc = kc - 1

           ! Periodic boundary conditions
           if (per) then
              if (ic .lt. 1) ic = ic + n1
              if (ic .gt. n1) ic = ic - n1
              if (jc .lt. 1) jc = jc + n2
              if (jc .gt. n2) jc = jc - n2
              if (kc .lt. 1) kc = kc + n3
              if (kc .gt. n3) kc = kc - n3
           ! Non-periodic
           else
              if (ic .lt. 1) ic = 0
              if (ic .gt. n1) ic = n1 + 1
              if (jc .lt. 1) jc = 0
              if (jc .gt. n2) jc = n2 + 1
              if (kc .lt. 1) kc = 0
              if (kc .gt. n3) kc = n3 + 1
           end if
           
           icp = int(ic) + 1
           jcp = int(jc) + 1
           kcp = int(kc) + 1

           ! Periodic boundary conditions
           if (per) then
              if (icp .lt. 1) icp = icp + n1
              if (icp .gt. n1) icp = icp - n1
              if (jcp .lt. 1) jcp = jcp + n2
              if (jcp .gt. n2) jcp = jcp - n2
              if (kcp .lt. 1) kcp = kcp + n3
              if (kcp .gt. n3) kcp = kcp - n3
           ! Non-periodic
           else
              if (icp .lt. 1) icp = 0
              if (icp .gt. n1) icp = n1 + 1
              if (jcp .lt. 1) jcp = 0
              if (jcp .gt. n2) jcp = n2 + 1
              if (kcp .lt. 1) kcp = 0
              if (kcp .gt. n3) kcp = n3 + 1
           end if
           
           vg(ic, jc, kc) = vg(ic, jc, kc) + vp(i, j, k) * xl * yl * zl
           vg(ic, jc, kcp) = vg(ic, jc, kcp) + vp(i, j, k) * xl * yl * zr
           vg(ic, jcp, kc) = vg(ic, jcp, kc) + vp(i, j, k) * xl * yr * zl
           vg(icp, jc, kc) = vg(icp, jc, kc) + vp(i, j, k) * xr * yl * zl
           vg(ic, jcp, kcp) = vg(ic, jcp, kcp) + vp(i, j, k) * xl * yr * zr
           vg(icp, jc, kcp) = vg(icp, jc, kcp) + vp(i, j, k) * xr * yl * zr
           vg(icp, jcp, kc) = vg(icp, jcp, kc) + vp(i, j, k) * xr * yr * zl
           vg(icp, jcp, kcp) = vg(icp, jcp, kcp) + vp(i, j, k) * xr * yr * zr
        end do
     end do
  end do

  ! deallocate(xg, yg, zg, vp)
  ! deallocate(xg)
  ! deallocate(yg)
  ! deallocate(zg)
  ! deallocate(vp)
  
  ! Write out the data
  open(24, file=trim(path)//'ic_'//trim(field)//'_cic', form='unformatted')
  write(24) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  do k=1,n3
     write(24) ((vg(i, j, k), i=1,n1), j=1,n2)
  end do
  close(24)

  ! deallocate(vg)
  
end subroutine cic_grid

subroutine gen_vbc(path)
  implicit none

  character (len=200), intent(in) :: path

  character (len=3) :: xyz = 'xyz'
  integer :: i, j, k, n1, n2, n3
  real :: dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  real, allocatable :: b_tmp(:, :), c_tmp(:, :)
  real, allocatable :: vbc(:, :, :)
  
  write(6, *) '---- reading ic_velbx'
  ! Read the x data
  open(20, file=trim(path)//'ic_velbx', form='unformatted')
  rewind 20
  ! Read header and print out cosmological parameters, for quick checking
  read(20) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  write(6, *) '-------- n1, n2, n3', n1, n2, n3
  write(6, *) '-------- astart', astart, ' zstart', 1.0/astart - 1.0
  write(6, *) '-------- omega_m', omega_m
  write(6, *) '-------- omega_l', omega_l
  write(6, *) '-------- H_0', h0
  call flush(6)
  close(20)
  
  ! Now we have n1, n2 and n3, we can allocate arrays
  allocate(b_tmp(n1, n2), c_tmp(n1, n2))
  allocate(vbc(n1, n2, n3))

  vbc = 0.0
  
  do i=1,3
     write(6, *) '---- working on vel', xyz(i:i)
     open(20, file=trim(path)//'ic_velb'//xyz(i:i), form='unformatted')
     open(21, file=trim(path)//'ic_velc'//xyz(i:i)//'_cic', form='unformatted')

     ! Skip header
     read(20)
     read(21)
    
     do k=1,n3
        read(20) b_tmp(:, :)
        read(21) c_tmp(:, :)

        vbc(:, :, k) = vbc(:, :, k) + (b_tmp - c_tmp) ** 2.
     end do

     close(20)
     close(21)
    
  end do
 
  vbc = sqrt(vbc)

  ! Write out the data
  open(22, file=trim(path)//'ic_vbc', form='unformatted')
  write(22) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  do k=1,n3
     write(22) ((vbc(i, j, k), i=1,n1), j=1,n2)
  end do
  close(22)
    
end subroutine gen_vbc
