      subroutine cfft2(init, idir, npfft, cr, work, cvc)
      
      implicit none
      integer init, idir, npfft
      complex*16 cr(0:npfft-1), cvc(0:npfft-1)
      double precision work(4*npfft)

c     Written to mimic homonymous Cray routine.
c     n.b. these dimensions are ok for a single transform at a time.
c     modify if c06fpf to be called with first argument > 1.

c     n.b. here input cr is overwritten!
c
c     idir gives the sign in the FFT exponent; init=1 initializes (no FFT is
c     performed)
c
c	Input: cr
c     On output: npfft coefficients in cvc(0:...). 
c                Indices 0 to npfft/2 correspond to harmonics 0 to npfft/2
c                Indices npfft/2 to npfft-1 correspond to harmonics -npfft/2 to -1
c
      integer i, ifail, n2, dn, tn, j, k
      double precision fac
	save fac
      
      dn = 2 * npfft + 1
      tn = 3 * npfft + 1
      n2 = npfft / 2

      ifail = 0
      if(init.eq.1)then
c      call c06frf(1, npfft, cr(0), cr(n2), 'i', work, work(dn)
      call c06frf(1, npfft, work(dn), work(tn), 'i', work, work(dn), ifail)
cPL      fac = dsqrt(dfloat(npfft))

      else if(init.eq.0)then
      fac = dsqrt(dfloat(npfft))
        do i = 0, npfft - 1        
        work(dn+i) = dreal(cr(i))
        work(tn+i) = dfloat(-idir) * dimag(cr(i))
        end do
c       Build vectors of real and im parts, stored in complex arrays:
        do i = 0, n2 - 1
        j = dn + 2 * i
        k = tn + 2 * i
        cr(i) = dcmplx(work(j),work(j+1))
        cr(i+n2) = dcmplx(work(k),work(k+1))
        end do
      
      call c06frf(1, npfft, cr(0), cr(n2), 's', work, work(dn), ifail)
      
c     This way, cr and cvc can be equivalenced in calling program:
      call dcopy(npfft, cr(0), 1, work(dn), 1)
      call dcopy(npfft, cr(n2), 1, work(tn), 1)
        if(idir.gt.0)then
c       case of reverse transform:
          do i = 0, npfft - 1
          cvc(i) = dcmplx(fac*work(dn+i), -fac*work(tn+i))
c          do i = 0, n2 - 1
c          cvc(2*i) = dcmplx(fac*dreal(cr(i)), -fac*dreal(cr(n2+i)))
c          cvc(2*i+1) = dcmplx(fac*dimag(cr(i)), -fac*dimag(cr(n2+i)))
          end do
        else
c       case of direct transform:
          do i = 0, npfft - 1
          cvc(i) = dcmplx(fac*work(dn+i), fac*work(tn+i))
c          do i = 0, n2 - 1
c          cvc(2*i) = dcmplx(fac*dreal(cr(i)), fac*dreal(cr(n2+i)))
c          cvc(2*i+1) = dcmplx(fac*dimag(cr(i)), fac*dimag(cr(n2+i)))
          end do
        end if
      end if      
      
      return
      end
