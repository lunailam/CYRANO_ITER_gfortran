      subroutine mcrcrz(nra,nca,a,lda,nrb,ncb,b,ldb,nrc,ncc,c,ldc)

c     multiplies two complex*16 rectangular matrices, AB.
c     based on call to blas ZGEMM.

      integer nra, nca, lda, nrb, ncb, ldb, nrc, ncc, ldc             
      complex*16 a(lda,*), b(ldb,*), c(ldc,*), cun, cze
     
      
      if(nca.ne.nrb)stop 'error mcrcr:nca.ne.nrb'
      cun=dcmplx(1.d0,0.d0)
      cze=dcmplx(0.d0,0.d0)
      call zgemm('No transpose','No transpose',nra,ncb,nca,cun,a,lda
     ;,b,ldb,cze,c,ldc)
      
      return
      end
