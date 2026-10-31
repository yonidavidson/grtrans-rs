program test_emis_extra
! Dump BREMS, SYNCHTHBREMS, SYNCHTHAVNOABS and the HYBRID compositions
! for validation against the Rust port. Uses the same input grid as
! test_polsynch, with a thermal density ne = nnth*1e6.
   use polsynchemis
   use emissivity, only: brememisHEROIC
   implicit none
   integer, parameter :: n = 12
   real(8) :: nu(n), ne(n), nnth(n), b(n), th(n), p(n), gmin(n), gmax
   real(8) :: tt(n)
   real(8) :: e(n, 11), kth(n, 11), kpl(n, 11), kbr(n, 11), ksy(n, 11)
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
   p    = 3.5d0
   gmin = 100.d0
   gmax = 1.0d5
   tt = (/1.0d10, 3.0d10, 1.0d11, 1.0d9, 3.0d11, 5.0d9, 2.0d10, 1.0d12, 1.0d8, 4.0d10, 1.0d11, 7.0d9/)
   call initialize_polsynchpl(4)

   write(6, '(A)') '# brems'
   e = 0d0
   call brememisHEROIC(nu, ne, tt, e)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do

   write(6, '(A)') '# synchthbrems'
   call synchemis(nu, ne, b, tt, ksy)
   kbr = 0d0
   call brememisHEROIC(nu, ne, tt, kbr)
   e = ksy + kbr
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do

   write(6, '(A)') '# synchthavnoabs'
   call synchemisnoabs(nu, ne, b, tt, e)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do

   call polsynchth(nu, ne, b, tt, th, kth)
   call polsynchpl(nu, nnth, b, th, p, gmin, gmax, kpl)

   write(6, '(A)') '# hybridthpl'
   e = kth + kpl
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do

   write(6, '(A)') '# hybridth'
   e = kth
   e(:, 5:11) = kth(:, 5:11) + kpl(:, 5:11)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do

   write(6, '(A)') '# hybridpl'
   e = kpl
   e(:, 5:11) = kth(:, 5:11) + kpl(:, 5:11)
   do i = 1, n
      write(6, '(11(ES24.16E3,1X))') e(i, :)
   end do
end program test_emis_extra
