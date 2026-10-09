program test_zroots
! Dump ZROOTS/LAGUER outputs for the polynomials geokerr uses, for
! grtrans-rs validation.
   implicit none
   integer, parameter :: ndeg = 4
   complex*16 :: c(5), roots(4), coefs(7), hroots(6)
   double precision :: a, l, q2, cc, dd, ee
   integer :: i

   ! pixel 1 of the test_geokerr fixture
   a = 0.9375d0
   l = 7.9667540809039652d0
   q2 = 139.39017252377502d0
   cc = a*a - q2 - l*l
   dd = 2.d0*((a-l)**2 + q2)
   ee = -a*a*q2
   c(1) = dcmplx(1.d0, 0.d0)
   c(2) = dcmplx(0.d0, 0.d0)
   c(3) = dcmplx(cc, 0.d0)
   c(4) = dcmplx(dd, 0.d0)
   c(5) = dcmplx(ee, 0.d0)
   write(6,'(A,3(ES24.16E3,1X))') '# cc dd ee ', cc, dd, ee
   write(6,'(A)') '# quartic'
   call zroots(c, 4, roots, .true.)
   do i = 1, 4
      write(6,'(2(ES24.16E3,1X))') dble(roots(i)), dimag(roots(i))
   end do

   coefs(1) = dcmplx(1.d0, 0.d0)
   coefs(2) = dcmplx(-cc/sqrt(ee), 0.d0)
   coefs(3) = dcmplx(-1.d0, 0.d0)
   coefs(4) = dcmplx(sqrt(ee)*(2.d0*cc/ee-(dd/ee)**2), 0.d0)
   coefs(5) = dcmplx(-1.d0, 0.d0)
   coefs(6) = dcmplx(-cc/sqrt(ee), 0.d0)
   coefs(7) = dcmplx(1.d0, 0.d0)
   write(6,'(A)') '# sextic1'
   call zroots(coefs, 6, hroots, .true.)
   do i = 1, 6
      write(6,'(2(ES24.16E3,1X))') dble(hroots(i)), dimag(hroots(i))
   end do

   coefs(1) = ee**(-3.d0)
   coefs(2) = -cc/ee**3.d0
   coefs(3) = -ee**(-2.d0)
   coefs(4) = -ee**(-2.d0)*(dd**2.d0/ee-2.d0*cc)
   coefs(5) = -1.d0/ee
   coefs(6) = -cc/ee
   coefs(7) = 1.d0
   write(6,'(A)') '# sextic2'
   call zroots(coefs, 6, hroots, .true.)
   do i = 1, 6
      write(6,'(2(ES24.16E3,1X))') dble(hroots(i)), dimag(hroots(i))
   end do
end program test_zroots
