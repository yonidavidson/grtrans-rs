program test_bessel
! Dump upstream Bessel function values for grtrans-rs kernel validation.
! Upstream routines take rank-1 arrays, so call with 1-element arrays.
   use bessel
   implicit none
   integer :: i
   real(8) :: xs(1), x
   real(8), parameter :: grid(9) = (/0.1d0, 0.5d0, 1.d0, 1.9d0, 2.d0, &
                                     2.5d0, 3.75d0, 5.d0, 10.d0/)
   do i = 1, 9
      x = grid(i)
      xs(1) = x
      write(6, '(5(ES24.16E3,1X))') x, besseli0(xs), besseli1(xs), &
                                    besselk0(xs), besselk1(xs)
   end do
   write(6, '(A)') '# kn'
   do i = 1, 9
      x = grid(i)
      xs(1) = x
      write(6, '(4(ES24.16E3,1X))') x, besselk(2, xs), besselk(3, xs), &
                                     besselk(5, xs)
   end do
end program
