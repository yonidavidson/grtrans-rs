program test_toy_fluid
! Dump TOY density/angular momentum and four-velocity for validation. The
! upstream aleph/bb are garbage (stale f%u / unassigned bb), so only the
! well-defined toy_vals and f%u are compared.
   use class_four_vector
   use fluid_model, only: fluid, initialize_fluid_model, del_fluid_model, get_fluid_vars
   use fluid_model_toy, only: init_toy, toy_vals
   implicit none
   type(four_vector), dimension(:), allocatable :: x0, k0
   type(fluid) :: f
   character(len=20) :: fname
   real(8) :: a
   integer :: i, np
   real(8), dimension(12) :: rgrid, thgrid, n, l

   a = 0.9375
   fname = 'TOY'
   np = 12
   do i = 1, np
      rgrid(i) = 2.0d0 + 18.0d0 * dble(i - 1) / dble(np - 1)
      thgrid(i) = 0.3d0 + 2.5d0 * dble(i - 1) / dble(np - 1)
   end do
   call init_toy(3.0d7, 1.7d0, 2.0d0, 0)
   call toy_vals(a, cos(thgrid), rgrid, n, l)
   write(6, '(A)') '# toy_vals'
   do i = 1, np
      write(6, '(4(ES24.16E3,1X))') rgrid(i), thgrid(i), n(i), l(i)
   end do

   allocate(x0(np)); allocate(k0(np))
   do i = 1, np
      x0(i)%data(1) = 0.d0
      x0(i)%data(2) = rgrid(i)
      x0(i)%data(3) = thgrid(i)
      x0(i)%data(4) = 0.d0
      k0(i)%data(1) = 1.d0
      k0(i)%data(2) = 0.1d0
      k0(i)%data(3) = 0.05d0
      k0(i)%data(4) = 0.2d0
   end do
   call initialize_fluid_model(f, fname, dble(a), np)
   call get_fluid_vars(x0, k0, dble(a), f)
   write(6, '(A)') '# toy_u'
   do i = 1, np
      write(6, '(4(ES24.16E3,1X))') rgrid(i), thgrid(i), f%u(i)%data(1), f%u(i)%data(4)
   end do
   call del_fluid_model(f)
end program test_toy_fluid
