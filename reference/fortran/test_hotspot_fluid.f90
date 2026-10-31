program test_hotspot_fluid
! Dump HOTSPOT density, field and four-velocity for all four bl06 field
! models, for validation against the Rust port.
   use class_four_vector
   use fluid_model, only: fluid, initialize_fluid_model, del_fluid_model, get_fluid_vars
   use fluid_model_hotspot, only: init_hotspot
   implicit none
   type(four_vector), dimension(:), allocatable :: x0, k0
   type(fluid) :: f
   character(len=20) :: fname
   real :: a
   integer :: i, bl, np, ib
   integer, dimension(4) :: bls
   real(8), dimension(12) :: rgrid, thgrid, phigrid

   a = 0.9375
   bls = (/1, 0, 2, 3/)
   fname = 'HOTSPOT'
   np = 12
   rgrid = (/6.0d0, 6.2d0, 5.8d0, 6.0d0, 6.0d0, 6.5d0, 7.0d0, 8.0d0, 6.0d0, 6.0d0, 3.0d0, 12.0d0/)
   thgrid = (/1.5708d0, 1.5708d0, 1.5708d0, 1.4d0, 1.7d0, 1.5708d0, 1.5708d0, 1.5708d0, 0.8d0, 2.5d0, 1.5708d0, 1.5708d0/)
   phigrid = (/-1.5708d0, -1.4708d0, -1.6708d0, -1.3708d0, -1.7708d0, -1.0708d0, &
        -0.5708d0, 0.4292d0, -1.5708d0, -1.5708d0, -1.5708d0, -1.5708d0/)
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
   do ib = 1, 4
      bl = bls(ib)
      call init_hotspot(rs=1.5, r0=6.0, n0=1.0e4, bl=bl)
      call initialize_fluid_model(f, fname, dble(a), np)
      call get_fluid_vars(x0, k0, dble(a), f)
      write(6, '(A,I2)') '# bl06=', bl
      do i = 1, np
         write(6, '(11(ES24.16E3,1X))') x0(i)%data(2), x0(i)%data(3), x0(i)%data(4), &
              f%rho(i), f%bmag(i), f%u(i)%data(1), f%u(i)%data(4), &
              f%b(i)%data(1), f%b(i)%data(2), f%b(i)%data(3), f%b(i)%data(4)
      end do
      call del_fluid_model(f)
   end do
end program test_hotspot_fluid
