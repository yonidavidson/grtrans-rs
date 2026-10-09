program test_elliptic
! Dump RF/RC/RD/RJ over argument grids for grtrans-rs validation.
   implicit none
   integer :: i, j, k
   real(8) :: xs(6) = (/0.01d0, 0.1d0, 0.5d0, 1.d0, 4.d0, 100.d0/)
   real(8) :: rfx(6) = (/0.02d0, 0.2d0, 0.7d0, 1.5d0, 3.d0, 25.d0/)
   real(8) :: rfy(6) = (/0.05d0, 0.3d0, 0.9d0, 2.d0, 5.d0, 50.d0/)
   real(8) :: rfp(6) = (/0.03d0, 0.4d0, 1.1d0, 2.5d0, 6.d0, 200.d0/)
   real(8) :: rf, rc, rd, rj
   write(6,'(A)') '# rf'
   do i = 1, 6
      do j = 1, 6
         write(6,'(4(ES24.16E3,1X))') rfx(i), rfy(i), xs(j), &
              rf(rfx(i), rfy(i), xs(j))
      end do
   end do
   write(6,'(A)') '# rc'
   do i = 1, 6
      do j = 1, 6
         write(6,'(3(ES24.16E3,1X))') xs(i), rfy(j), rc(xs(i), rfy(j))
      end do
   end do
   write(6,'(A)') '# rd'
   do i = 1, 6
      do j = 1, 6
         write(6,'(4(ES24.16E3,1X))') xs(i), rfy(j), rfp(i), &
              rd(xs(i), rfy(j), rfp(i))
      end do
   end do
   write(6,'(A)') '# rj'
   do i = 1, 6
      do j = 1, 6
         write(6,'(5(ES24.16E3,1X))') xs(i), rfy(j), rfx(i), rfp(j), &
              rj(xs(i), rfy(j), rfx(i), rfp(j))
      end do
   end do
end program test_elliptic
