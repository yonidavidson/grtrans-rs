program test_interp
! Dump interpemis_noabs output for a synthetic table, to validate the
! log-log interpolation and the flattened (point-major, freq-fastest) layout.
   use emissivity, only: interpemis_noabs
   implicit none
   integer, parameter :: np = 5, nf = 4
   real(8) :: nu(np), freqarr(nf), jnu(np*nf), K(np, 11)
   integer :: p, f
   freqarr = (/1.0d10, 1.0d11, 1.0d12, 1.0d13/)
   ! jnu(point, freq) = point_index * 10^(freq_index-1), flattened point-major
   do p = 1, np
      do f = 1, nf
         jnu((p-1)*nf + f) = dble(p) * 10.0d0**(f-1)
      end do
   end do
   ! nu per point: inside brackets, at edges, below and above the grid
   nu = (/5.0d10, 1.0d11, 3.0d12, 1.0d9, 1.0d14/)
   K = 0.d0
   call interpemis_noabs(nu, freqarr, jnu, K)
   do p = 1, np
      write(6, '(ES24.16E3,1X,ES24.16E3)') nu(p), K(p, 1)
   end do
end program test_interp
