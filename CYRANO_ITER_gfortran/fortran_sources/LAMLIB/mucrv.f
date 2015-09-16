      subroutine mucrv(nra,nca,a,lda,nx,x,ipath,ny,y)

c     multiplies a complex rectangular matrix by a complex vector.
c     based on call to blas CGEMV.
c     ipath=1: y = a.x
c     ipath=2: y = transpose(a).x
c     else: no action.

      integer nra, nca, lda, nx, ipath, ny             
      complex a(lda,*), x(*), y(*), cun, cze
     
      cun=cmplx(1.e0,0.e0)
      cze=cmplx(0.e0,0.e0)
          if(ipath.eq.1)then
      if(nx.ne.nca)stop 'error mucrv:nx.ne.nca'
      if(ny.ne.nra)stop 'error mucrv:ny.ne.nra'
      call cgemv('No transpose',nra,nca,cun,a,lda,x,1,cze,y,1)
          else if(ipath.eq.2)then
      if(nx.ne.nra)stop 'error mucrv:nx.ne.nra'
      if(ny.ne.nca)stop 'error mucrv:ny.ne.nca'
      call cgemv('Transpose',nra,nca,cun,a,lda,x,1,cze,y,1)
          end if
      
      return
      end
