program test_harm_fluid
! Dump HARM fluid variables (interpolated) for validation.
   use class_four_vector
   use fluid_model_harm, only: initialize_harm_model, harm_vals, del_harm_data
   implicit none
   type (four_vector), dimension(:), allocatable :: x0, b, u
   real, dimension(:), allocatable :: rho, p, bmag
   integer, parameter :: nr = 12, nth = 12
   real(8) :: rr, tt, rgrid(12), thgrid(12)
   integer :: i, k, n
   real :: a
   character(len=20) :: ifile
   character(len=40) :: df, hf

   a = 0.9375
   ifile = 'harm_test.in'
   df = 'dump'
   hf = 'dump040'
   n = nr * nth
   allocate(x0(n)); allocate(b(n)); allocate(u(n))
   allocate(rho(n)); allocate(p(n)); allocate(bmag(n))

   do i = 1, nr
      ! descending: first point at the largest radius so the upstream
      ! time-index logic stays in range (tindx=0) for the single snapshot
      rgrid(i) = 39.0d0 * exp(-dble(i - 1) / dble(nr - 1) * log(39.0d0 / 1.5d0))
   end do
   do k = 1, nth
      thgrid(k) = 0.05d0 + (3.0d0 - 0.05d0) * dble(k - 1) / dble(nth - 1)
   end do

   call initialize_harm_model(dble(a), ifile, df, hf, 1, 40)

   n = 0
   do k = 1, nth
      do i = 1, nr
         n = n + 1
         rr = rgrid(i)
         tt = thgrid(k)
         x0(n)%data(1) = 0.d0
         x0(n)%data(2) = rr
         x0(n)%data(3) = tt
         x0(n)%data(4) = 0.d0
      end do
   end do
   call harm_vals(x0, a, rho, p, b, u, bmag)
   write(6, '(A)') '# rho p bmag bb u0 u1 u2 u3 b0 b1 b2 b3 r th'
   do n = 1, nr * nth
      write(6, '(14(ES24.16E3,1X))') rho(n), p(n), bmag(n), b(n) * b(n), &
           u(n)%data(1), u(n)%data(2), u(n)%data(3), u(n)%data(4), &
           b(n)%data(1), b(n)%data(2), b(n)%data(3), b(n)%data(4), &
           x0(n)%data(2), x0(n)%data(3)
   end do
   call del_harm_data()
end program test_harm_fluid
