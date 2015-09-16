      subroutine mcrcr(nra,nca,a,lda,nrb,ncb,b,ldb,nrc,ncc,c,ldc)

c     multiplies two complex rectangular matrices, AB.
c     based on call to blas CGEMM.

      integer nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc             
      complex a(lda,*), b(ldb,*), c(ldc,*), cun, cze
     
      
      if(nca.ne.nrb)stop 'error mcrcr:nca.ne.nrb'
      cun=cmplx(1.e0,0.e0)
      cze=cmplx(0.e0,0.e0)
      call cgemm('No transpose','No transpose',nra,ncb,nca,cun,a,lda
     ;,b,ldb,cze,c,ldc)
      
      return
      end
