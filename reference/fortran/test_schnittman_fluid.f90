program test_schnittman_fluid
! Dump SCHNITTMAN spot density, toroidal field and four-velocity for
! validation against the Rust port.
   use class_four_vector
   use fluid_model, only: fluid, initialize_fluid_model, del_fluid_model, get_fluid_vars
   use fluid_model_hotspot_schnittman, only: init_schnittman_hotspot
   implicit none
   type(four_vector), dimension(:), allocatable :: x0, k0
   type(fluid) :: f
   character(len=20) :: fname
   real :: a
   integer :: i, np
   real(8), dimension(12) :: rgrid, thgrid, phigrid

   a = 0.9375
   fname = 'SCHNITTMAN'
   np = 12
   rgrid = (/6.0d0, 6.2d0, 5.8d0, 6.0d0, 6.0d0, 6.5d0, 7.0d0, 8.0d0, 6.0d0, 6.0d0, 3.0d0, 12.0d0/)
   thgrid = (/1.5708d0, 1.5708d0, 1.5708d0, 1.4d0, 1.7d0, 1.5708d0, 1.5708d0, 1.5708d0, 0.8d0, 2.5d0, 1.5708d0, 1.5708d0/)
   phigrid = (/0.0d0, 0.1d0, -0.1d0, 0.2d0, -0.2d0, 0.5d0, 1.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
   call init_schnittman_hotspot(rs=1.5, r0=6.0, n0=1.0e4)
   allocate(x0(np)); allocate(k0(np))
   do i = 1, np
      x0(i)%data(1) = 0.d0
      x0(i)%data(2) = rgrid(i)
      x0(i)%data(3) = thgrid(i)
      x0(i)%data(4) = phigrid(i)
      k0(i)%data(1) = 1.d0
      k0(i)%data(2) = 0.1d0
      k0(i)%data(3) = 0.05d0
      k0(i)%data(4) = 0.2d0
   end do
   call initialize_fluid_model(f, fname, dble(a), np)
   call get_fluid_vars(x0, k0, dble(a), f)
   write(6, '(A)') '# schnittman'
   do i = 1, np
      write(6, '(11(ES24.16E3,1X))') x0(i)%data(2), x0(i)%data(3), x0(i)%data(4), &
           f%rho(i), f%bmag(i), f%u(i)%data(1), f%u(i)%data(4), &
           f%b(i)%data(1), f%b(i)%data(2), f%b(i)%data(3), f%b(i)%data(4)
   end do
   call del_fluid_model(f)
end program test_schnittman_fluid
