program test_geomu
! Dump GEOMU outputs for pixel 11 of the test_geokerr fixture.
! Self-contained: geokerr_wrapper.f has all needed routines.
   implicit none
   double precision :: u0, uf, mu0, muf, a, l, l2, q2, su, sm
   integer :: iu, tpm, tpr, ncase
   double precision :: h1, u1, u2, u3, u4, rffu0, rffu1, rffmu1, rffmu2
   double precision :: rffmu3, iu0, i1mu, i3mu, iu_d
   logical :: pht, firstpt
   integer :: k
   double precision :: ufs(5)
   double precision :: phimu, tmu, phiu, tu, lambda
   double precision :: rdc, rjc, tu01, tu02, tu03, tu04, tmu1, tmu3
   double precision :: phimu1, phimu3
   interface
      subroutine geomu(u0,uf,mu0,muf,a,l,l2,q2,iu,tpm,tpr,su,sm,ncase,h1, &
           u1,u2,u3,u4,rffu0,rffu1,rffmu1,rffmu2,rffmu3,iu0,i1mu,i3mu,pht,firstpt)
         double precision, intent(in) :: u0, mu0, a, l, l2, q2, su, sm
         double precision, intent(inout) :: uf, muf
         integer, intent(inout) :: tpm, tpr
         double precision, intent(out) :: h1, u1, u2, u3, u4
         double precision, intent(out) :: rffu0, rffu1, rffmu1, rffmu2, rffmu3
         double precision, intent(out) :: iu0, i1mu, i3mu
         integer, intent(out) :: ncase
         double precision, intent(out) :: iu
         logical, intent(in) :: pht, firstpt
      end subroutine geomu
      subroutine geophitime(u0,uf,mu0,muf,a,l,l2,q2,tpm,tpr,su,sm,iu,h1, &
           phimu,tmu,ncase,u1,u2,u3,u4,phiu,tu,lambda,rffu0,rffu1,rffmu1, &
           rffmu2,rffmu3,rdc,rjc,tu01,tu02,tu03,tu04,tmu1,tmu3,phimu1, &
           phimu3,firstpt)
         double precision, intent(in) :: u0, uf, mu0, muf, a, l, l2, q2
         integer, intent(in) :: tpm, tpr, ncase
         double precision, intent(in) :: su, sm, iu, h1, u1, u2, u3, u4
         double precision, intent(in) :: rffu0, rffu1, rffmu1, rffmu2, rffmu3
         double precision, intent(out) :: phimu, tmu, phiu, tu, lambda
         double precision, intent(out) :: rdc, rjc, tu01, tu02, tu03, tu04
         double precision, intent(out) :: tmu1, tmu3, phimu1, phimu3
         logical, intent(in) :: firstpt
      end subroutine geophitime
   end interface

   ! pixel 11 of camera1
   u0 = 1.7506389832288786d-7
   mu0 = 0.6428d0
   a = 0.9375d0
   l = 0.0d0
   l2 = 0.0d0
   q2 = 1.0199343109375000d1
   su = 1.0d0
   sm = -1.0d0
   tpm = 1
   tpr = 0
   pht = .true.
   ! uf values taken from geo1 pixel 11 (u = 1/r at the sampled points)
   ufs = (/1.0d0/1.8505853198045095d1, 1.0d0/1.2178655213497432d1, &
           1.0d0/9.075663069891803d0, 1.0d0/7.232817822530681d0, &
           1.0d0/6.012049970184252d0/)
   write(6,'(A)') '# geomu'
   do k = 1, 5
      uf = ufs(k)
      muf = 0.0d0
      call geomu(u0,uf,mu0,muf,a,l,l2,q2,iu_d,tpm,tpr,su,sm,ncase,h1, &
           u1,u2,u3,u4,rffu0,rffu1,rffmu1,rffmu2,rffmu3,iu0,i1mu,i3mu,pht,.true.)
      write(6,'(A,I2)') '# call ', k
      write(6,'(7(ES24.16E3,1X))') uf, iu_d, muf, u1, u2, u3, u4
      write(6,'(8(ES24.16E3,1X))') rffu0, rffu1, rffmu1, rffmu2, rffmu3, &
           iu0, i1mu, i3mu
      write(6,'(2I6,ES24.16E3,1X)') ncase, tpr, h1
      ! geophitime with the geomu outputs (firstpt=.true.)
      rffu1 = 0.d0
      call geophitime(u0,uf,mu0,muf,a,l,l2,q2,tpm,tpr,su,sm,iu_d,h1, &
           phimu,tmu,ncase,u1,u2,u3,u4,phiu,tu,lambda,rffu0,rffu1, &
           rffmu1,rffmu2,rffmu3,rdc,rjc,tu01,tu02,tu03,tu04,tmu1,tmu3, &
           phimu1,phimu3,.true.)
      write(6,'(A,I2)') '# phitime ', k
      write(6,'(5(ES24.16E3,1X))') phimu, tmu, phiu, tu, lambda
   end do
end program test_geomu
