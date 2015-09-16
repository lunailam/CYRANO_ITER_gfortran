      subroutine lfccg(n,a,lda,fac,ldfac,ipvt,rcond)

c     computes the LU factorization of a general complex matrix and 
c     estimates its L1 condition number.

      integer n,lda,ldfac,ipvt(*),info
      real rcond
      complex a(lda,*), fac(ldfac,*)
c pas fini!      
c      call ccopy
c      call onenorm
      call cgetrf(n,n,fac,ldfac,ipvt,info)
      call cgecon('1',n,fac,ldfac,anorm,rcond,work,rwork,info) 
      
      return
      end
