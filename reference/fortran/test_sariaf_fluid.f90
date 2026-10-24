program test_sariaf_fluid
! Dump SARIAF fluid variables for validation.
   use class_four_vector
   use fluid_model, only: fluid, initialize_fluid_model, del_fluid_model, get_fluid_vars
   use fluid_model_sariaf, only: init_sariaf
   implicit none
   type(four_vector), dimension(:), allocatable :: x0, k0
   type(fluid) :: f
   integer, parameter :: nr = 12, nth = 12
   real(8) :: rgrid(12), thgrid(12)
   real :: a
   integer :: i, k, n
   character(len=20) :: fname

   a = 0.9375
   fname = 'SARIAF'
   call init_sariaf(3.0e7, 1.7e11, 8.0e4, 2.9, 10.0, 0)

   n = nr * nth
   allocate(x0(n)); allocate(k0(n))
   do i = 1, nr
      rgrid(i) = 1.5d0 * exp(dble(i - 1) / dble(nr - 1) * log(40.0d0 / 1.5d0))
   end do
   do k = 1, nth
      thgrid(k) = 0.05d0 + (3.0d0 - 0.05d0) * dble(k - 1) / dble(nth - 1)
   end do
   n = 0
   do k = 1, nth
      do i = 1, nr
         n = n + 1
         x0(n)%data(1) = 0.d0
         x0(n)%data(2) = rgrid(i)
         x0(n)%data(3) = thgrid(k)
         x0(n)%data(4) = 0.d0
         k0(n)%data(1) = 0.d0
         k0(n)%data(2) = 0.d0
         k0(n)%data(3) = 0.d0
         k0(n)%data(4) = 0.d0
      end do
   end do
   call initialize_fluid_model(f, fname, dble(a), n)
   call get_fluid_vars(x0, k0, dble(a), f)
   write(6, '(A)') '# rho p bmag rho2 u0 u1 u2 u3 b0 b1 b2 b3 r th'
   do n = 1, nr * nth
      write(6, '(14(ES24.16E3,1X))') f%rho(n), f%p(n), f%bmag(n), f%rho2(n), &
           f%u(n)%data(1), f%u(n)%data(2), f%u(n)%data(3), f%u(n)%data(4), &
           f%b(n)%data(1), f%b(n)%data(2), f%b(n)%data(3), f%b(n)%data(4), &
           x0(n)%data(2), x0(n)%data(3)
   end do
   call del_fluid_model(f)
end program test_sariaf_fluid
