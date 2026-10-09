program test_polsynch
! Dump polarized power-law synchrotron coefficients for validation.
   use polsynchemis
   implicit none
   integer, parameter :: n = 12
   real(8) :: nu(n), nnth(n), b(n), th(n), p(n), gmin(n), gmax
   real(8) :: e(n, 11)
   integer :: i
   ! values spanning the FFJET regime
   nu   = (/2.7555415006877594d11, 2.683007249434544d11, 2.731003087487894d11, &
            3.45d11, 3.45d11, 3.45d11, 1.0d11, 1.0d12, 2.0d11, 5.0d11, &
            3.45d11, 3.45d11/)
   nnth = (/4.107198715209961d-1, 1.6513969898223877d0, 4.701847076416016d0, &
            1.0d0, 1.0d0, 1.0d0, 0.1d0, 10.0d0, 1.5d0, 0.5d0, 1.0d0, 1.0d0/)
   b    = (/2.221096083521843d0, 3.065144419670105d0, 4.255986772477627d0, &
            3.0d0, 3.0d0, 3.0d0, 0.5d0, 20.0d0, 2.0d0, 5.0d0, 1.0d0, 10.0d0/)
   th   = (/3.368731396732495d-1, 2.6194631579424915d-1, 1.8918080719155692d-1, &
            0.3d0, 1.2d0, 0.05d0, 0.5d0, 1.0d0, 0.8d0, 0.2d0, 1.5d0, 0.1d0/)
   p    = 3.5d0
   gmin = 100.d0
   gmax = 1.0d5
   call initialize_polsynchpl(4)
   call polsynchpl(nu, nnth, b, th, p, gmin, gmax, e)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do
   ! also the unpolarized variant
   write(6, '(A)') '# synchpl'
   call synchpl(nu, nnth, b, th, p, gmin, gmax, e)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do
end program test_polsynch
