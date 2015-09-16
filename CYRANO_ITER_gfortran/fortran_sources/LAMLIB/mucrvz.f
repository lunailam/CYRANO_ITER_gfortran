      subroutine mucrvz(nra,nca,a,lda,nx,x,ipath,ny,y)

c     multiplies a complex*16 rectangular matrix by a complex*16 vector.
c     based on call to blas ZGEMV.
c     ipath=1: y = a.x
c     ipath=2: y = transpose(a).x
c     else: no action.

      integer nra, nca, lda, nx, ipath, ny             
      complex*16 a(lda,*), x(*), y(*)

      character*1 trans
      complex*16 cone, czero

c      call mucrv(nra,nca,a,lda,nx,x,ipath,ny,y)
      trans = ' '
      cone = dcmplx(1.d0,0.d0)
      czero = dcmplx(0.d0,0.d0)

      if (ipath .eq. 1) then
        trans = 'N'
      elseif (ipath .eq. 2) then
        trans = 'T'
      endif

      call zgemv(trans, nra, nca, cone, a, lda, x, 1, czero, y, 1)

      return
      end
