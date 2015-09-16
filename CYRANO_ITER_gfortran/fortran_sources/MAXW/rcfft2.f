      subroutine rcfft2(init, idir, npfft, cr, work, cvc)
      
      implicit none
      integer i, init, idir, npfft, ifail, n2
      double precision cr(0:npfft-1), fac
      complex*16 work(3*npfft/2+2), cvc(0:npfft/2)

c     n.b. these dimensions are ok for a single transform at a time.
c     modify if c06fpf to be called with first argument > 1.

c     n.b. here input cr is overwritten!
      
      ifail = 0
      if(init.eq.1)then
      call c06fpf(1, npfft, work(npfft+1), 'i', work, work(npfft+1), ifail)
c      call c06fpf(1, npfft, cr(0), 'i', work, work(npfft+1), ifail)

      else if(init.eq.0)then
      call c06fpf(1, npfft, cr(0), 's', work, work(npfft+1), ifail)
      
c     Special normaliz. for Cray RCFFT2:
      fac = 2.d0*dsqrt(dfloat(npfft))
      n2 = npfft / 2
      
      cvc(0) = dcmplx(fac*cr(0), 0.d0)
      cvc(n2) = dcmplx(fac*cr(n2),0.d0)

        if(idir.gt.0)then
c       case of reverse transform:
          do i = 1, n2 - 1
          cvc(i) = dcmplx(fac*cr(i), -fac*cr(npfft-i))
          end do
        else
c       case of direct transform:
          do i = 1, npfft/2-1
          cvc(i) = dcmplx(fac*cr(i), fac*cr(npfft-i))
          end do
        end if
      end if
            
      return
      end
