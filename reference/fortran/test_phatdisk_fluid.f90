program test_phatdisk_fluid
! Dump PHATDISK tables, interpolated values and the Keplerian fluid velocity
! for validation against the Rust port. The thindisk rin/rout are set
! explicitly (upstream leaves them uninitialized for a standalone PHATDISK).
   use class_four_vector
   use fluid_model_thindisk, only: init_thindisk
   use fluid_model_phatdisk, only: init_phatdisk, phatdisk_vals, fnu_tab, &
        r_tab, om_tab, freq_tab, nfreq_tab, nr
   use fluid_model, only: fluid, initialize_fluid_model, del_fluid_model, get_fluid_vars
   implicit none
   type(four_vector), dimension(:), allocatable :: x0, k0
   type(fluid) :: f
   character(len=20) :: fname, ifile
   real :: a
   integer :: i, k, n, np
   real, dimension(12) :: rgrid, thgrid
   real, dimension(12) :: rvals
   real, dimension(:,:), allocatable :: fnu
   real, dimension(:), allocatable :: omega

   a = 0.9375
   fname = 'PHATDISK'
   ifile = 'phatdisk_test.in'
   open(unit=12, file=ifile, form='formatted', status='replace')
   write(12, '(A)') '&phatdisk'
   write(12, '(A)') ' nw=20, wmin=1e-4, wmax=1e4, nfreq_tab=10, fmin=1e10, fmax=1e15,'
   write(12, '(A)') ' rmax=10000.0, nr=30, sigt=0.4, fcol=1.7,'
   write(12, '(A)') '/'
   write(12, '(A)') '&thindisk'
   write(12, '(A)') ' mbh=620000000.0,'
   write(12, '(A)') ' mdot=1570000000000000.0,'
   write(12, '(A)') '/'
   close(12)
   call init_thindisk(a, mdott=1.57e15, mbht=6.2e8, rint=-1.0, routt=1e8)
   call init_phatdisk(a, ifile)

   write(6, '(A)') '# tables'
   do i = 1, nr, 3
      do k = 1, nfreq_tab, 2
         write(6, '(I4,1X,I4,1X,4(ES24.16E3,1X))') i, k, r_tab(i), freq_tab(k), fnu_tab(i, k), om_tab(i)
      end do
   end do

   np = 12
   do i = 1, np
      rgrid(i) = 2.0 * exp(dble(i - 1) / dble(np - 1) * log(50.0d0 / 2.0d0))
      thgrid(i) = 0.3d0 + 2.5d0 * dble(i - 1) / dble(np - 1)
   end do
   rvals = real(rgrid)
   allocate(fnu(np, nfreq_tab)); allocate(omega(np))
   call phatdisk_vals(rvals, a, fnu, omega)
   write(6, '(A)') '# vals'
   do i = 1, np
      write(6, '(4(ES24.16E3,1X))') rgrid(i), omega(i), fnu(i, 1), fnu(i, nfreq_tab)
   end do

   allocate(x0(np)); allocate(k0(np))
   do i = 1, np
      x0(i)%data(1) = 0.d0
      x0(i)%data(2) = rgrid(i)
      x0(i)%data(3) = thgrid(i)
      x0(i)%data(4) = 0.d0
      k0(i)%data(1) = 1.d0
      k0(i)%data(2) = 0.d0
      k0(i)%data(3) = 0.d0
      k0(i)%data(4) = 0.d0
   end do
   call initialize_fluid_model(f, fname, dble(a), np)
   call get_fluid_vars(x0, k0, dble(a), f)
   write(6, '(A)') '# fluid'
   do i = 1, np
      write(6, '(6(ES24.16E3,1X))') x0(i)%data(2), x0(i)%data(3), f%u(i)%data(1), f%u(i)%data(4), f%fnu(i, 1), freq_tab(1)
   end do
   call del_fluid_model(f)
end program test_phatdisk_fluid
