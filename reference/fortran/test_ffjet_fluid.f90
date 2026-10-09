program test_ffjet_fluid
! Dump FFJET fluid variables for validation.
   use class_four_vector
   use fluid_model_ffjet, only: initialize_ffjet_model, ffjet_vals, del_ffjet_data
   use kerr, only: calc_rms
   implicit none
   type (four_vector), dimension(:), allocatable :: x0, b, u
   real, dimension(:), allocatable :: rho, p, bmag
   real :: a
   integer, parameter :: nr = 12, nz = 12
   real(8) :: rr, zz, rgrid(12), zgrid(12)
   integer :: i, k, n
   character(len=40) :: df

   df = 'm87bl09rfp10xi5a998fluidvars.bin'
   a = 0.998
   n = nr * nz
   allocate(x0(n)); allocate(b(n)); allocate(u(n))
   allocate(rho(n)); allocate(p(n)); allocate(bmag(n))

   ! log grid in cylindrical (r,z) covering the jet region
   do i = 1, nr
      rgrid(i) = 1.2d0 * exp(dble(i - 1) / dble(nr - 1) * log(120.0d0 / 1.2d0))
   end do
   do k = 1, nz
      zgrid(k) = 0.3d0 * exp(dble(k - 1) / dble(nz - 1) * log(120.0d0 / 0.3d0))
   end do

   call initialize_ffjet_model(dble(a), df=df)

   n = 0
   do k = 1, nz
      do i = 1, nr
         n = n + 1
         rr = rgrid(i)
         zz = zgrid(k)
         x0(n)%data(1) = 0.d0
         x0(n)%data(2) = sqrt(rr * rr + zz * zz)
         x0(n)%data(3) = acos(zz / sqrt(rr * rr + zz * zz))
         x0(n)%data(4) = 0.d0
      end do
   end do
   call ffjet_vals(x0, a, rho, p, b, u, bmag)
   write(6, '(A)') '# rho p bmag bb u0 u1 u2 u3 b0 b1 b2 b3'
   do n = 1, nr * nz
      write(6, '(13(ES24.16E3,1X))') rho(n), p(n), bmag(n), b(n) * b(n), &
           u(n)%data(1), u(n)%data(2), u(n)%data(3), u(n)%data(4), &
           b(n)%data(1), b(n)%data(2), b(n)%data(3), b(n)%data(4), &
           x0(n)%data(2)
   end do
   call del_ffjet_data()
end program test_ffjet_fluid
