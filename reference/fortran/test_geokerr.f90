program test_geokerr
! Dump geokerr camera + geodesic outputs for grtrans-rs validation,
! using the production code path (geodesics module: initialize_pixels +
! initialize_geodesic), exactly as pgrtrans/grtrans_driver use it.
!
! Sections:
!   [camera1]/[geo1]   standard=1 (radial tracing), nup=25, mufill
!   [camera1b]/[geo1b] standard=1, nup=400 (production-like)
!   [camera2]/[geo2]   standard=2 (polar tracing), nup=1
!
! Determinism note: upstream leaves TPMARR unassigned for standard=1 and
! TPRARR unassigned for standard=2; we explicitly zero them here (and the
! Rust port does the same) so the fixture is reproducible.
   use geodesics
   use class_four_vector
   implicit none

   type (geokerr_args) :: gargs
   type (geo) :: g
   integer :: status, i, k, ipix
   integer, parameter :: nro1 = 5, nphi1 = 4, npix1 = nro1 * nphi1
   integer, parameter :: nro2 = 5, nphi2 = 4, npix2 = nro2 * nphi2

   ! ---------------- standard = 1, nup = 25 ----------------
   call initialize_geokerr_args(gargs, npix1)
   call initialize_pixels(gargs, .true., 1, 0.6428d0, -0.5d0, 0.9375d0, &
        0.04d0, 1.d0, 1.d0, 2, -13.d0, 13.d0, -13.d0, 13.d0, &
        nro1, nphi1, 25)
   gargs%t0 = 0.d0
   gargs%tpm = 0
   gargs%tpr = 0
   gargs%tpr = gargs%tpr + 1  ! standard=1 camera sets TPR=1
   write(6,'(A)') '# camera1'
   write(6,'(A,2(ES24.16E3,1X))') '# u0 offset ', gargs%u0, gargs%offset
   do ipix = 1, npix1
      write(6,'(8(ES24.16E3,1X),2I6)') gargs%alpha(ipix), &
           gargs%beta(ipix), gargs%q2(ipix), gargs%l(ipix), &
           gargs%uf(ipix), gargs%muf(ipix), gargs%sm(ipix), &
           gargs%su(ipix), gargs%tpm(ipix), gargs%tpr(ipix)
   end do
   write(6,'(A)') '# geo1'
   do ipix = 1, npix1, 5
      call initialize_geodesic(g, gargs, ipix, status)
      write(6,'(A,I6,1X,ES24.16E3,1X,4I6)') '# pixel ', ipix, &
           g%gk%uf(1), g%gk%tpm(1), g%gk%tpr(1), status, g%npts
      do k = 1, g%npts
         write(6,'(10(ES24.16E3,1X))') g%x(k)%data(1), g%x(k)%data(2), &
              g%x(k)%data(3), g%x(k)%data(4), g%k(k)%data(1), &
              g%k(k)%data(2), g%k(k)%data(3), g%k(k)%data(4), &
              g%lambda(k), dble(g%tpmarr(k) * 1000 + g%tprarr(k))
      end do
      call del_geodesic(g)
   end do
   call del_geokerr_args(gargs)

   ! ---------------- standard = 1, nup = 400 ----------------
   call initialize_geokerr_args(gargs, npix1)
   call initialize_pixels(gargs, .true., 1, 0.6428d0, -0.5d0, 0.9375d0, &
        0.04d0, 1.d0, 1.d0, 2, -13.d0, 13.d0, -13.d0, 13.d0, &
        nro1, nphi1, 400)
   gargs%t0 = 0.d0
   gargs%tpm = 0
   gargs%tpr = 1
   write(6,'(A)') '# camera1b'
   write(6,'(A,2(ES24.16E3,1X))') '# u0 offset ', gargs%u0, gargs%offset
   do ipix = 1, npix1
      write(6,'(8(ES24.16E3,1X),2I6)') gargs%alpha(ipix), &
           gargs%beta(ipix), gargs%q2(ipix), gargs%l(ipix), &
           gargs%uf(ipix), gargs%muf(ipix), gargs%sm(ipix), &
           gargs%su(ipix), gargs%tpm(ipix), gargs%tpr(ipix)
   end do
   write(6,'(A)') '# geo1b'
   do ipix = 1, npix1, 5
      call initialize_geodesic(g, gargs, ipix, status)
      write(6,'(A,I6,1X,ES24.16E3,1X,4I6)') '# pixel ', ipix, &
           g%gk%uf(1), g%gk%tpm(1), g%gk%tpr(1), status, g%npts
      do k = 1, g%npts
         write(6,'(10(ES24.16E3,1X))') g%x(k)%data(1), g%x(k)%data(2), &
              g%x(k)%data(3), g%x(k)%data(4), g%k(k)%data(1), &
              g%k(k)%data(2), g%k(k)%data(3), g%k(k)%data(4), &
              g%lambda(k), dble(g%tpmarr(k) * 1000 + g%tprarr(k))
      end do
      call del_geodesic(g)
   end do
   call del_geokerr_args(gargs)

   ! ---------------- standard = 2, nup = 1 ----------------
   call initialize_geokerr_args(gargs, npix2)
   call initialize_pixels(gargs, .true., 2, 0.26d0, -0.5d0, 0.9d0, &
        0.01d0, 1.d0, 1.d0, 2, -21.d0, 21.d0, -21.d0, 21.d0, &
        nro2, nphi2, 1)
   gargs%t0 = 0.d0
   ! standard=2: TPMARR is set by the camera; TPRARR is left unassigned
   ! upstream, so zero it for determinism.
   gargs%tpr = 0
   write(6,'(A)') '# camera2'
   write(6,'(A,2(ES24.16E3,1X))') '# u0 offset ', gargs%u0, gargs%offset
   do ipix = 1, npix2
      write(6,'(8(ES24.16E3,1X),2I6)') gargs%alpha(ipix), &
           gargs%beta(ipix), gargs%q2(ipix), gargs%l(ipix), &
           gargs%uf(ipix), gargs%muf(ipix), gargs%sm(ipix), &
           gargs%su(ipix), gargs%tpm(ipix), gargs%tpr(ipix)
   end do
   write(6,'(A)') '# geo2'
   do ipix = 1, npix2
      call initialize_geodesic(g, gargs, ipix, status)
      write(6,'(A,I6,1X,ES24.16E3,1X,4I6)') '# pixel ', ipix, &
           g%gk%uf(1), g%gk%tpm(1), g%gk%tpr(1), status, g%npts
      do k = 1, g%npts
         write(6,'(10(ES24.16E3,1X))') g%x(k)%data(1), g%x(k)%data(2), &
              g%x(k)%data(3), g%x(k)%data(4), g%k(k)%data(1), &
              g%k(k)%data(2), g%k(k)%data(3), g%k(k)%data(4), &
              g%lambda(k), dble(g%tpmarr(k) * 1000 + g%tprarr(k))
      end do
      call del_geodesic(g)
   end do
   call del_geokerr_args(gargs)

end program test_geokerr
