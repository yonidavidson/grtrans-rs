program test_maxwell
! Dump MAXJUTT and MAXCOMP (with and without a selection) coefficients for
! validation against the Rust port. Same input grid as test_emis_extra.
   use calc_maxjutt, only: calc_maxjutt_subroutine
   use calc_maxcomp, only: calc_maxcomp_subroutine
   implicit none
   integer, parameter :: n = 12
   real(8) :: nu(n), ne(n), nnth(n), b(n), th(n), tt(n)
   real(8) :: e(n, 11)
   real(8) :: args_j(7), args_c6(8), args_c0(8)
   integer :: i
   nu   = (/2.7555415006877594d11, 2.683007249434544d11, 2.731003087487894d11, &
            3.45d11, 3.45d11, 3.45d11, 1.0d11, 1.0d12, 2.0d11, 5.0d11, &
            3.45d11, 3.45d11/)
   nnth = (/4.107198715209961d-1, 1.6513969898223877d0, 4.701847076416016d0, &
            1.0d0, 1.0d0, 1.0d0, 0.1d0, 10.0d0, 1.5d0, 0.5d0, 1.0d0, 1.0d0/)
   ne = nnth * 1.0d6
   b    = (/2.221096083521843d0, 3.065144419670105d0, 4.255986772477627d0, &
            3.0d0, 3.0d0, 3.0d0, 0.5d0, 20.0d0, 2.0d0, 5.0d0, 1.0d0, 10.0d0/)
   th   = (/3.368731396732495d-1, 2.6194631579424915d-1, 1.8918080719155692d-1, &
            0.3d0, 1.2d0, 0.05d0, 0.5d0, 1.0d0, 0.8d0, 0.2d0, 1.5d0, 0.1d0/)
   tt = (/1.0d10, 3.0d10, 1.0d11, 1.0d9, 3.0d11, 5.0d9, 2.0d10, 1.0d12, 1.0d8, 4.0d10, 1.0d11, 7.0d9/)

   args_j = (/2.0d0, 0.10714295d0, 0.4459917d0, 1.34219997d0, 2.17450422d0, &
              1.36515333d0, 0.52892829d0/)
   args_c6 = (/2.0d0, 6.0d0, 0.10714295d0, 0.4459917d0, 1.34219997d0, &
               2.17450422d0, 1.36515333d0, 0.52892829d0/)
   args_c0 = (/2.0d0, 0.0d0, 0.10714295d0, 0.4459917d0, 1.34219997d0, &
               2.17450422d0, 1.36515333d0, 0.52892829d0/)

   write(6, '(A)') '# maxjutt'
   call calc_maxjutt_subroutine(nu, ne, b, tt, th, args_j, e)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do

   write(6, '(A)') '# maxcomp_sel6'
   call calc_maxcomp_subroutine(nu, ne, b, tt, th, args_c6, e)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do

   write(6, '(A)') '# maxcomp_sel0'
   call calc_maxcomp_subroutine(nu, ne, b, tt, th, args_c0, e)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do
end program test_maxwell
