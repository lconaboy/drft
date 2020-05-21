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
!    - Displacements in grafic ic_posc files are taken to be in units
!      of comoving Mpc/h, while the cell widths (dxini in the grafic
!      headers) are in units of comoving Mpc. This is discrepant with
!      the original definition of the grafic2 file (the README
!      available from Edmund Bertschinger's page), which defines
!      displacements in comoving Mpc. After working through the units
!      in pm/init_part.f90 (and the comments therein) as well as
!      looking at other sources (pynbody's grafic reading routines,
!      for instance), it looks as though ramses expects the ic_posc
!      files to be in comoving Mpc/h, so that's what we'll use here.
!
!    - There is some ambiguity in how to handle the edges for the
!    - non-periodic case, which become sharp due to there being no
!    - mass outside either edge to smooth them. For the CDM
!    - overdensity, I chose to set those edge cases (calculated as
!    - anything between the edge and the maximum particle-cell offset)
!    - to zero, and for the velocities I set them as their ungridded
!    - velocities. Obviously, this is not ideal but I think it is a
!    - better result than having a highly underdense border to the
!    - zoom region. My decision is somewhat mitigated by the fact that
!    - the refinement maps usually don't go all the way from 1 -> n,
!    - instead they have a region which is roughly the same size as
!    - the edge that I set to zero (determined by checking). A future
!    - improvement could be to set these edge cases from the previous
!    - level of refinement, which is still not ideal but may offer a
!    - slight advantage. Due to current time constraints I haven't
!    - implemented that.
!
!        + The edge values are now averaged out from the central
!        values, which removes the discontinuity introduced earlier by
!        a sudden change in value.
!
! L Conaboy, April 2020


subroutine gen_delc(path, omega_b, per)
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
  logical, intent(in) :: per
  
  integer :: hblk = 44
  integer :: mblk
  integer :: i, j, k
  integer :: ip, jp, kp, ic, jc, kc, n1, n2, n3
  integer :: icp, jcp, kcp
  integer :: ii, jj, kk, cut, pad!, cnt
  integer, parameter :: f=50
  double precision :: mp, dxx, dyy, dzz, tx, ty, tz, dxx1, dyy1, dzz1
  !real :: mp, dxx, dyy, dzz, tx, ty, tz, dxx1, dyy1, dzz1
  real :: omega_c, boxvol
  real :: dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  double precision :: rho_cr = 2.775d11  ! h^2 Msol/Mpc^3
  ! real :: rho_cr = 2.775e11  ! h^2 Msol/Mpc^3
  double precision :: cur_max, cur_min, dh, ddxini
  ! real :: cur_max, cur_min, old_maxval
  real, allocatable, dimension(:, :, :) :: dx, dy, dz, del_c, cnt

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

  ! Define double convenience variables
  dh = dble(h0) / 100.d0
  ddxini = dble(dxini)
  
  ! Now calculate the particle mass
  rho_cr = rho_cr * dh**2
  omega_c = omega_m - omega_b
  boxvol = real(n1*n2*n3) * (dxini**3)

  ! Calculate mp
  mp = (dble(omega_c) * rho_cr * dble(boxvol)) / dble(n1*n2*n3)  ! M_sol


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
  
  ! Interpolate the particles onto a grid
  write(6, *) '---- interpolating'
  
  if (per) then
     write(6, *) '-------- using periodic boundary conditions'
     call flush(6)
          
     allocate(del_c(n1, n2, n3))
     del_c = 0.0
 
     cur_min = 0.0
     cur_max = 0.0
     
     ! Now loop over particles
     do ip = 1, n1
        do jp = 1, n2
           do kp = 1, n3
              ! Convert the displcaments from comoving Mpc/h to units
              ! of cell width: dxx=1.0 means the particle is offset by
              ! 1 cell from the current cell centre
              dxx = dble(dx(ip, jp, kp)) / ddxini / dh !* dx0 / dxini
              dyy = dble(dy(ip, jp, kp)) / ddxini / dh !* dx0 / dxini
              dzz = dble(dz(ip, jp, kp)) / ddxini / dh !* dx0 / dxini

              ! We can use displacements to figure out which cell the
              ! particle is in, since we've normalised to cell widths
              ic = int(dble(ip) + dxx + 0.5d0)
              jc = int(dble(jp) + dyy + 0.5d0)
              kc = int(dble(kp) + dzz + 0.5d0)

              ! Calcaulate the offset from the new cell centre
              dxx1 = dxx - (dble(ic) - dble(ip) - 0.5d0)  
              dyy1 = dyy - (dble(jc) - dble(jp) - 0.5d0)
              dzz1 = dzz - (dble(kc) - dble(kp) - 0.5d0)

              ! Check if we've wrapped the ic value around
              if (ic .gt. n1) ic = ic - n1
              if (ic .lt. 1) ic = ic + n1
              if (jc .gt. n2) jc = jc - n2
              if (jc .lt. 1) jc = jc + n2
              if (kc .gt. n3) kc = kc - n3
              if (kc .lt. 1) kc = kc + n3
              
              ! Now calculate the +1th cell, and check if we've
              ! wrapped around
              icp = ic + 1
              jcp = jc + 1
              kcp = kc + 1
              if (icp .gt. n1) icp = icp - n1
              if (icp .lt. 1) icp = icp + n1
              if (jcp .gt. n2) jcp = jcp - n2
              if (jcp .lt. 1) jcp = jcp + n2
              if (kcp .gt. n3) kcp = kcp - n3
              if (kcp .lt. 1) kcp = kcp + n3

              ! Check how large the offset is, this gives us an idea
              ! of how much of the edge to remove in the ICs, for the
              ! non-periodic case. Just used for information in the
              ! periodic case.
              if (dxx .lt. cur_min) then
                 cur_min = dxx
              end if
              if (dyy .lt. cur_min) then
                 cur_min = dyy
              end if
              if (dzz .lt. cur_min) then
                 cur_min = dzz
              end if
              
              if (dxx .gt. cur_max) then
                 cur_max = dxx
              end if
              if (dyy .gt. cur_max) then
                 cur_max = dyy
              end if
              if (dzz .gt. cur_max) then
                 cur_max = dzz
              end if
              
              ! Convenience variables
              tx = 1.0 - dxx1
              ty = 1.0 - dyy1
              tz = 1.0 - dzz1

              ! Interpolate using cloud-in-cell
              del_c(ic, jc, kc) = del_c(ic, jc, kc) + real(mp*tx*ty*tz)
              del_c(icp, jc, kc) = del_c(icp, jc, kc) + real(mp*dxx1*ty*tz)
              del_c(ic, jcp, kc) = del_c(ic, jcp, kc) + real(mp*tx*dyy1*tz)
              del_c(ic, jc, kcp) = del_c(ic, jc, kcp) + real(mp*tx*ty*dzz1)
              del_c(icp, jcp, kc) = del_c(icp, jcp, kc) + real(mp*dxx1*dyy1*tz)
              del_c(icp, jc, kcp) = del_c(icp, jc, kcp) + real(mp*dxx1*ty*dzz1)
              del_c(ic, jcp, kcp) = del_c(ic, jcp, kcp) + real(mp*tx*dyy1*dzz1)
              del_c(icp, jcp, kcp) = del_c(icp, jcp, kcp) + real(mp*dxx1*dyy1*dzz1)
           end do
        end do
     end do

     write(6, *) '---- done interpolating'
     call flush(6)


     write(6, *) '-------- min(offset)', cur_min
     write(6, *) '-------- max(offset)', cur_max

     ! Convert del_c to proper density
     del_c = del_c / dxini**3  ! M_sol Mpc^-3

     ! Convert to overdensity
     del_c = del_c/(rho_cr * omega_c) - 1.

     write(6, *) '-------- min(del_c) ', minval(del_c)
     write(6, *) '-------- max(del_c) ', maxval(del_c)
     
     deallocate(dx)
     deallocate(dy)
     deallocate(dz)

  else
     write(6, *) '-------- using non-periodic boundary conditions'
     call flush(6)
     
     allocate(del_c(-1:n1+2, -1:n2+2, -1:n3+2))     
     del_c = 0.0

     cur_min = 0.0
     cur_max = 0.0
     
     ! Now loop over particles
     do ip = 1, n1
        do jp = 1, n2
           do kp = 1, n3
              ! non-periodic CIC, we build up the 
              ! xp - particle indices
              ! xc - cell indices
              ! xcp - cell indices + 1

              ! Convert the displcaments from comoving Mpc/h to units
              ! of cell width: dxx=1.0 means the particle is offset by
              ! 1 cell from the current cell centre
              dxx = dble(dx(ip, jp, kp)) / ddxini / dh !* dx0 / dxini
              dyy = dble(dy(ip, jp, kp)) / ddxini / dh !* dx0 / dxini
              dzz = dble(dz(ip, jp, kp)) / ddxini / dh !* dx0 / dxini


              ! We can use d** to figure out which cell the particle is
              ! in, since we've normalised to cell widths. Might have
              ! gone out of the right hand side, so use mod
              ic = int(dble(ip) + dxx + 0.5d0)
              jc = int(dble(jp) + dyy + 0.5d0)
              kc = int(dble(kp) + dzz + 0.5d0)
           
              ! Calcaulate the offset from the new cell centre
              dxx1 = dxx - (dble(ic) - dble(ip) - 0.5d0)  
              dyy1 = dyy - (dble(jc) - dble(jp) - 0.5d0)
              dzz1 = dzz - (dble(kc) - dble(kp) - 0.5d0)

              ! Check how large the offset is, this gives us an idea
              ! of how much of the edge to remove in the ICs
              if (dxx .lt. cur_min) then
                 cur_min = dxx
              end if
              if (dyy .lt. cur_min) then
                 cur_min = dyy
              end if
              if (dzz .lt. cur_min) then
                 cur_min = dzz
              end if
              
              if (dxx .gt. cur_max) then
                 cur_max = dxx
              end if
              if (dyy .gt. cur_max) then
                 cur_max = dyy
              end if
              if (dzz .gt. cur_max) then
                 cur_max = dzz
              end if

              ! Calculate the values for +1 before shifting
              icp = ic + 1
              jcp = jc + 1
              kcp = kc + 1
              
              ! Might have also gone out of the left hand side, if so
              ! stick all of the mass in the zeroth index, which we
              ! will get rid of anyway
              if (ic .lt. 0) ic = -1 ! n1 + ic
              if (jc .lt. 0) jc = -1 ! n2 + jc
              if (kc .lt. 0) kc = -1 ! n3 + kc
              ! Similarly for the rhs
              if (ic .gt. n1) ic = n1 + 1! ic - n1 - 1
              if (jc .gt. n2) jc = n2 + 1! jc - n2 - 1
              if (kc .gt. n3) kc = n3 + 1! kc - n3 - 1
              
              ! Now we can adjust the +1 values, starting with the
              ! lhs, < 1 corresponds to ic < 0
              if (icp .lt. 1) icp = 0 ! n1 + icp
              if (jcp .lt. 1) jcp = 0 ! n2 + jcp
              if (kcp .lt. 1) kcp = 0 ! n3 + kcp
              ! and now the rhs, again icp > n1+1 corresponds to ic >
              ! n1
              if (icp .gt. n1+1) icp = n1 + 2 ! icp - n1 - 1
              if (jcp .gt. n2+1) jcp = n2 + 2 ! jcp - n2 - 1
              if (kcp .gt. n3+1) kcp = n3 + 2 ! kcp - n3 - 1

              ! Convenience variables
              tx = 1.0 - dxx1
              ty = 1.0 - dyy1
              tz = 1.0 - dzz1

              ! Interpolate using cloud-in-cell
              del_c(ic, jc, kc) = del_c(ic, jc, kc) + real(mp*tx*ty*tz)
              del_c(icp, jc, kc) = del_c(icp, jc, kc) + real(mp*dxx1*ty*tz)
              del_c(ic, jcp, kc) = del_c(ic, jcp, kc) + real(mp*tx*dyy1*tz)
              del_c(ic, jc, kcp) = del_c(ic, jc, kcp) + real(mp*tx*ty*dzz1)
              del_c(icp, jcp, kc) = del_c(icp, jcp, kc) + real(mp*dxx1*dyy1*tz)
              del_c(icp, jc, kcp) = del_c(icp, jc, kcp) + real(mp*dxx1*ty*dzz1)
              del_c(ic, jcp, kcp) = del_c(ic, jcp, kcp) + real(mp*tx*dyy1*dzz1)
              del_c(icp, jcp, kcp) = del_c(icp, jcp, kcp) + real(mp*dxx1*dyy1*dzz1)
           end do
        end do
     end do

     write(6, *) '---- done interpolating'
     call flush(6)

     ! Remove edges
     cut = int(ceiling(max(cur_max, abs(cur_min)))) + 1
     ! del_c = del_c(1:n1, 1:n2, 1:n3)

     write(6, *) '-------- min(del_c) ', minval(del_c)
     write(6, *) '-------- max(del_c) ', maxval(del_c)
     write(6, *) '-------- min(offset)', cur_min
     write(6, *) '-------- max(offset)', cur_max
     write(6, *) '------------ cut', cut

     ! Convert del_c to proper density
     del_c = del_c / dxini**3  ! M_sol Mpc^-3

     ! Convert to overdensity
     del_c = del_c/(rho_cr * omega_c) - 1.

     ! Deal with the sharp edges caused by non-periodic interp. Set to
     ! zero initially.
     del_c(:, :, 1:cut) = 0.0
     del_c(:, 1:cut, :) = 0.0
     del_c(1:cut, :, :) = 0.0
     del_c(n1-cut+1:n1, :, :) = 0.0
     del_c(:, n2-cut+1:n2, :) = 0.0
     del_c(:, :, n3-cut+1:n3) = 0.0

     deallocate(dx)
     deallocate(dy)
     deallocate(dz)

     ! Count for the averages
     allocate(cnt(n1, n2, n3))
     cnt = 1.0
     
     ! Now average outwards in one dimension, from the zeroed cells
     do i = cut, 1, -1
        do j = 1, n2
           do k = 1, n3
             ! Check lhs edge
              del_c(i, j, k) = del_c(i, j, k) + del_c(i+1, j, k)
              cnt(i, j, k) = cnt(i, j, k) + 1.
              
              ! Check rhs edge
              ii = n1 - i + 1
              jj = j
              kk = k

              del_c(ii, jj, kk) =  del_c(ii, jj, kk) + del_c(ii-1, j, k)
              cnt(ii, jj, kk) = cnt(ii, jj, kk) + 1.
           end do
        end do
     end do
    
     do i = 1, n1
        do j = cut, 1, -1
           do k = 1, n3
              del_c(i, j, k) = del_c(i, j, k) + del_c(i, j+1, k)
              cnt(i, j, k) = cnt(i, j, k) + 1.

              ! Check rhs edge
              ii = i
              jj = n2 - j + 1
              kk = k

              del_c(ii, jj, kk) =  del_c(ii, jj, kk) + del_c(ii, jj-1, kk)
              cnt(ii, jj, kk) = cnt(ii, jj, kk) + 1.
              
           end do
        end do
     end do

     do i = 1, n1
        do j = 1, n2
           do k = cut, 1, -1
              ! Check lhs edge
              del_c(i, j, k) = del_c(i, j, k) + del_c(i, j, k+1)
              cnt(i, j, k) = cnt(i, j, k) + 1.

              ! Check rhs edge
              ii = i
              jj = j
              kk = n3 - k + 1

              del_c(ii, jj, kk) =  del_c(ii, jj, kk) + del_c(ii, jj, kk-1)
              cnt(ii, jj, kk) = cnt(ii, jj, kk) + 1.
           end do
        end do
     end do

     ! Average for the edges
     del_c = del_c(1:n1, 1:n2, 1:n3) / cnt
     
     write(6, *) '-------- each cell modified on average', sum(cnt) / real(3*n1*n2*cut)
     write(6, *) '-------- min(del_c) ', minval(del_c)
     write(6, *) '-------- max(del_c) ', maxval(del_c)
     
     
  end if
  
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
  ! Write the header, don't need to do this if we don't use
  ! access='stream'
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

  ! if (.not. per) then
  !    write(6, *) '---- writing cnt'
  !    call flush(6)
  !    open(f, file=trim(path)//'ic_cnt', form='unformatted') !, access='stream')
  !    rewind f
  !    ! Write the header, don't need to do this if we don't use
  !    ! access='stream'
  !    ! write(f) hblk
  !    write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  !    ! write(f) hblk
  !    ! Write the main block
  !    do k=1,n3
  !       ! write(f) mblk
  !       write(f) ((cnt(i, j, k), i=1,n1), j=1,n2)
  !       ! write(f) mblk
  !    end do
  !    close(f)
  !    deallocate(cnt)
  ! end if

  deallocate(del_c)

  write(6, *) '---- done'
  
end subroutine gen_delc


subroutine gen_velcg(path, per)
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

  implicit none
  
  character (len=200), intent(in) :: path
  logical, intent(in) :: per
  
  character (len=3) :: xyz
  integer :: hblk = 44
  integer :: mblk
  integer :: ic, jc, kc, icp, jcp, kcp, n1, n2, n3
  integer :: i, j, k, ip, jp, kp, ll
  integer :: ii, jj, kk, cut, pad!, cnt
  integer, parameter :: f=50
  double precision :: dxx, dyy, dzz, tx, ty, tz, dxx1, dyy1, dzz1
  real :: dxini, x1off, x2off, x3off, astart, omega_m, omega_l, omega_b, h0
  double precision :: cur_max, cur_min, ddxini, dh
  double precision, parameter :: rhoc = 2.775d11  ! h^2 Msol/Mpc^3
  real, allocatable, dimension(:, :, :) :: dx, dy, dz, vc, vcg, cnt

  xyz = 'xyz'
  
  write(6, *) 'Interpolating CDM velocities'

  ! Read the x data
  write(6, *) '---- reading parameters'
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

  ! Define double convenience variables
  dh = dble(h0) / 100.d0
  ddxini = dble(dxini)

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
     allocate(vc(n1, n2, n3))

     ! Read in the the CDM velocity
     write(6, *) '---- reading velc'//xyz(ll:ll)
     open(f, file=trim(path)//'ic_velc'//xyz(ll:ll), form='unformatted')
     read(f)
     do k = 1, n3
        read(f) ((vc(i, j, k), i=1,n1), j=1,n2)
     end do
     close(f)

     ! Now interpolate
     write(6, *) '---- interpolating'
     call flush(6)
     if (per) then
        write(6, *) '-------- using periodic boundary conditions'
        call flush(6)
          
        allocate(vcg(n1, n2, n3))
        vcg = 0.0
 
        cur_min = 0.0
        cur_max = 0.0
     
        ! Now loop over particles
        do ip = 1, n1
           do jp = 1, n2
              do kp = 1, n3
                 ! periodic CIC
                 ! xp - particle indices
                 ! xc - cell indices
                 ! xcp - cell indices + 1

                 ! Convert the displcaments from comoving Mpc/h to
                 ! units of cell width: dxx=1.0 means the particle is
                 ! offset by 1 cell from the current cell centre
                 dxx = dble(dx(ip, jp, kp)) / ddxini / dh !* dx0 / dxini
                 dyy = dble(dy(ip, jp, kp)) / ddxini / dh !* dx0 / dxini
                 dzz = dble(dz(ip, jp, kp)) / ddxini / dh !* dx0 / dxini

                 ! We can use d** to figure out which cell the particle is
                 ! in, since we've normalised to cell widths. Might have
                 ! gone out of the right hand side, so use mod
                 ic = int(dble(ip) + dxx + 0.5d0)
                 jc = int(dble(jp) + dyy + 0.5d0)
                 kc = int(dble(kp) + dzz + 0.5d0)

                 ! Calcaulate the offset from the new cell centre
                 dxx1 = dxx - (dble(ic) - dble(ip) - 0.5d0)  
                 dyy1 = dyy - (dble(jc) - dble(jp) - 0.5d0)
                 dzz1 = dzz - (dble(kc) - dble(kp) - 0.5d0)

                 ! Check if we've wrapped the ic value around
                 if (ic .gt. n1) ic = ic - n1
                 if (ic .lt. 1) ic = ic + n1
                 if (jc .gt. n2) jc = jc - n2
                 if (jc .lt. 1) jc = jc + n2
                 if (kc .gt. n3) kc = kc - n3
                 if (kc .lt. 1) kc = kc + n3
              
                 ! Now calculate the +1th cell, and check if we've
                 ! wrapped around
                 icp = ic + 1
                 jcp = jc + 1
                 kcp = kc + 1
                 if (icp .gt. n1) icp = icp - n1
                 if (icp .lt. 1) icp = icp + n1
                 if (jcp .gt. n2) jcp = jcp - n2
                 if (jcp .lt. 1) jcp = jcp + n2
                 if (kcp .gt. n3) kcp = kcp - n3
                 if (kcp .lt. 1) kcp = kcp + n3
              
                 ! Check how large the offset is, this gives us an
                 ! idea of how much of the edge to remove in the ICs
                 if (dxx .lt. cur_min) then
                    cur_min = dxx
                 end if
                 if (dyy .lt. cur_min) then
                    cur_min = dyy
                 end if
                 if (dzz .lt. cur_min) then
                    cur_min = dzz
                 end if
              
                 if (dxx .gt. cur_max) then
                    cur_max = dxx
                 end if
                 if (dyy .gt. cur_max) then
                    cur_max = dyy
                 end if
                 if (dzz .gt. cur_max) then
                    cur_max = dzz
                 end if
              
                 ! Convenience variables
                 tx = 1.0 - dxx1
                 ty = 1.0 - dyy1
                 tz = 1.0 - dzz1

                 ! Interpolate using cloud-in-cell
                 vcg(ic, jc, kc) = vcg(ic, jc, kc) + vc(ip, jp, kp)*real(tx*ty*tz)
                 vcg(icp, jc, kc) = vcg(icp, jc, kc) + vc(ip, jp, kp)*real(dxx1*ty*tz)
                 vcg(ic, jcp, kc) = vcg(ic, jcp, kc) + vc(ip, jp, kp)*real(tx*dyy1*tz)
                 vcg(ic, jc, kcp) = vcg(ic, jc, kcp) + vc(ip, jp, kp)*real(tx*ty*dzz1)
                 vcg(icp, jcp, kc) = vcg(icp, jcp, kc) + vc(ip, jp, kp)*real(dxx1*dyy1*tz)
                 vcg(icp, jc, kcp) = vcg(icp, jc, kcp) + vc(ip, jp, kp)*real(dxx1*ty*dzz1)
                 vcg(ic, jcp, kcp) = vcg(ic, jcp, kcp) + vc(ip, jp, kp)*real(tx*dyy1*dzz1)
                 vcg(icp, jcp, kcp) = vcg(icp, jcp, kcp) + vc(ip, jp, kp)*real(dxx1*dyy1*dzz1)
              end do
           end do
        end do

        write(6, *) '---- done interpolating'
        call flush(6)


        write(6, *) '-------- min(offset)', cur_min
        write(6, *) '-------- max(offset)', cur_max
        write(6, *) '-------- min(vcg) ', minval(vcg)
        write(6, *) '-------- max(vcg) ', maxval(vcg)

        deallocate(vc, vcg)
     else
        write(6, *) '-------- using non-periodic boundary conditions'
        call flush(6)
     
        allocate(vcg(-1:n1+2, -1:n2+2, -1:n3+2))     
        vcg = 0.0

        cur_min = 0.0
        cur_max = 0.0
     
        ! Now loop over particles
        do ip = 1, n1
           do jp = 1, n2
              do kp = 1, n3
                 ! non - periodic CIC, we build up the 
                 ! xp - particle indices
                 ! xc - cell indices
                 ! xcp - cell indices + 1

                 
                 ! Convert the displcaments from comoving Mpc/h to
                 ! units of cell width: dxx=1.0 means the particle is
                 ! offset by 1 cell from the current cell centre
                 dxx = dble(dx(ip, jp, kp)) / ddxini / dh !* dx0 / dxini
                 dyy = dble(dy(ip, jp, kp)) / ddxini / dh !* dx0 / dxini
                 dzz = dble(dz(ip, jp, kp)) / ddxini / dh !* dx0 / dxini

                 ! We can use d** to figure out which cell the particle is
                 ! in, since we've normalised to cell widths. Might have
                 ! gone out of the right hand side, so use mod
                 ic = int(dble(ip) + dxx + 0.5d0)
                 jc = int(dble(jp) + dyy + 0.5d0)
                 kc = int(dble(kp) + dzz + 0.5d0)
           
                 ! Calcaulate the offset from the new cell centre
                 dxx1 = dxx - (dble(ic) - dble(ip) - 0.5d0)  
                 dyy1 = dyy - (dble(jc) - dble(jp) - 0.5d0)
                 dzz1 = dzz - (dble(kc) - dble(kp) - 0.5d0)

                 ! Check how large the offset is, this gives us an idea
                 ! of how much of the edge to remove in the ICs
                 if (dxx .lt. cur_min) then
                    cur_min = dxx
                 end if
                 if (dyy .lt. cur_min) then
                    cur_min = dyy
                 end if
                 if (dzz .lt. cur_min) then
                    cur_min = dzz
                 end if
              
                 if (dxx .gt. cur_max) then
                    cur_max = dxx
                 end if
                 if (dyy .gt. cur_max) then
                    cur_max = dyy
                 end if
                 if (dzz .gt. cur_max) then
                    cur_max = dzz
                 end if
              
                 ! Calculate the values for +1 before shifting
                 icp = ic + 1
                 jcp = jc + 1
                 kcp = kc + 1
              
                 ! Might have also gone out of the left hand side, if
                 ! so stick all of the mass in the zeroth index, which
                 ! we will get rid of anyway
                 if (ic .lt. 0) ic = -1 ! n1 + ic
                 if (jc .lt. 0) jc = -1 ! n2 + jc
                 if (kc .lt. 0) kc = -1 ! n3 + kc
                 ! Similarly for the rhs
                 if (ic .gt. n1) ic = n1 + 1! ic - n1 - 1
                 if (jc .gt. n2) jc = n2 + 1! jc - n2 - 1
                 if (kc .gt. n3) kc = n3 + 1! kc - n3 - 1
              
                 ! Now we can adjust the +1 values, starting with the
                 ! lhs
                 if (icp .lt. 1) icp = 0 ! n1 + icp
                 if (jcp .lt. 1) jcp = 0 ! n2 + jcp
                 if (kcp .lt. 1) kcp = 0 ! n3 + kcp
                 ! and now the rhs
                 if (icp .gt. n1+1) icp = n1 + 2 ! icp - n1 - 1
                 if (jcp .gt. n2+1) jcp = n2 + 2 ! jcp - n2 - 1
                 if (kcp .gt. n3+1) kcp = n3 + 2 ! kcp - n3 - 1

                 ! Convenience variables
                 tx = 1.0 - dxx1
                 ty = 1.0 - dyy1
                 tz = 1.0 - dzz1

                 ! Interpolate using cloud-in-cell
                 vcg(ic, jc, kc) = vcg(ic, jc, kc) + vc(ip, jp, kp)*real(tx*ty*tz)
                 vcg(icp, jc, kc) = vcg(icp, jc, kc) + vc(ip, jp, kp)*real(dxx1*ty*tz)
                 vcg(ic, jcp, kc) = vcg(ic, jcp, kc) + vc(ip, jp, kp)*real(tx*dyy1*tz)
                 vcg(ic, jc, kcp) = vcg(ic, jc, kcp) + vc(ip, jp, kp)*real(tx*ty*dzz1)
                 vcg(icp, jcp, kc) = vcg(icp, jcp, kc) + vc(ip, jp, kp)*real(dxx1*dyy1*tz)
                 vcg(icp, jc, kcp) = vcg(icp, jc, kcp) + vc(ip, jp, kp)*real(dxx1*ty*dzz1)
                 vcg(ic, jcp, kcp) = vcg(ic, jcp, kcp) + vc(ip, jp, kp)*real(tx*dyy1*dzz1)
                 vcg(icp, jcp, kcp) = vcg(icp, jcp, kcp) + vc(ip, jp, kp)*real(dxx1*dyy1*dzz1)
              end do
           end do
        end do
        write(6, *) '---- done interpolating'
        ! Remove edges
        cut = int(ceiling(max(cur_max, abs(cur_min)))) + 1
        ! vcg = vcg(1:n1, 1:n2, 1:n3)

        write(6, *) '-------- min(vcg) ', minval(vcg)
        write(6, *) '-------- max(vcg) ', maxval(vcg)
        write(6, *) '-------- min(offset)', cur_min
        write(6, *) '-------- max(offset)', cur_max
        write(6, *) '------------ cut', cut


        write(6, *) '-------- max displacement (cell width = 1.0)', cur_max
        write(6, *) '-------- min displacement (cell width = 1.0)', cur_min
                 
        ! Because of how we treated the periodic boundary conditions, we
        ! can either set the edge cells to zero, or use the ungridded
        ! velocity. Both are incorrect, but I think using the ungridded
        ! value is less wrong than setting to zero.
        vcg(:, :, 1:cut) = vc(:, :, 1:cut)
        vcg(:, 1:cut, :) = vc(:, 1:cut, :)
        vcg(1:cut, :, :) = vc(1:cut, :, :)
        vcg(n1-cut+1:n1, :, :) = vc(n1-cut+1:n1, :, :)
        vcg(:, n2-cut+1:n2, :) = vc(:, n2-cut+1:n2, :)
        vcg(:, :, n3-cut+1:n3) = vc(:, :, n3-cut+1:n3)

        deallocate(dx)
        deallocate(dy)
        deallocate(dz)

        ! Count for the averages
        allocate(cnt(n1, n2, n3))
        cnt = 1.0
     
        ! Now average outwards in one dimension, from the zeroed cells
        do i = cut, 1, -1
           do j = 1, n2
              do k = 1, n3
                 ! Check lhs edge
                 vcg(i, j, k) = vcg(i, j, k) + vcg(i+1, j, k)
                 cnt(i, j, k) = cnt(i, j, k) + 1.
              
                 ! Check rhs edge
                 ii = n1 - i + 1
                 jj = j
                 kk = k

                 vcg(ii, jj, kk) =  vcg(ii, jj, kk) + vcg(ii-1, j, k)
                 cnt(ii, jj, kk) = cnt(ii, jj, kk) + 1.
              end do
           end do
        end do
    
        do i = 1, n1
           do j = cut, 1, -1
              do k = 1, n3
                 vcg(i, j, k) = vcg(i, j, k) + vcg(i, j+1, k)
                 cnt(i, j, k) = cnt(i, j, k) + 1.

                 ! Check rhs edge
                 ii = i
                 jj = n2 - j + 1
                 kk = k

                 vcg(ii, jj, kk) =  vcg(ii, jj, kk) + vcg(ii, jj-1, kk)
                 cnt(ii, jj, kk) = cnt(ii, jj, kk) + 1.
              
              end do
           end do
        end do
        
        do i = 1, n1
           do j = 1, n2
              do k = cut, 1, -1
                 ! Check lhs edge
                 vcg(i, j, k) = vcg(i, j, k) + vcg(i, j, k+1)
                 cnt(i, j, k) = cnt(i, j, k) + 1.

                 ! Check rhs edge
                 ii = i
                 jj = j
                 kk = n3 - k + 1

                 vcg(ii, jj, kk) =  vcg(ii, jj, kk) + vcg(ii, jj, kk-1)
                 cnt(ii, jj, kk) = cnt(ii, jj, kk) + 1.
              end do
           end do
        end do

        vcg = vcg(1:n1, 1:n2, 1:n3) / cnt
 
        ! Write out the gridded velocity
        write(6, *) '---- writing velcg'//xyz(ll:ll)
        open(f, file=trim(path)//'ic_velcg'//xyz(ll:ll), form='unformatted')
        write(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
        do k = 1, n3
           write(f) ((vcg(i, j, k), i=1,n1), j=1,n2)
        end do
        close(f)
        deallocate(vc, vcg, cnt)
     end if
  end do
  
  ! Clean up
  deallocate(dx, dy, dz)
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
  real :: dxini, x1off, x2off, x3off, astart, omega_m, omega_l, omega_b, h0
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
  read(f) n1, n2, n3, dxini, x1off, x2off, x3off, astart, omega_m, omega_l, h0
  write(6, *) '-------- n1, n2, n3', n1, n2, n3
  write(6, *) '-------- astart', astart, ' zstart', 1.0/astart - 1.0
  write(6, *) '-------- omega_m', omega_m
  write(6, *) '-------- omega_l', omega_l
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
