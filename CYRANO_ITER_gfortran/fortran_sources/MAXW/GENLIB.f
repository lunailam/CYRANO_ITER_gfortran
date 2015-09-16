cERN	NOT USED
      subroutine cumul(a, b, w, n)

      implicit none

      integer n, i, j
      complex*16 a(n,n), b(n,n)
      complex*16 w

      do 1 i=1, n
      do 1 j=1, n
      a(i,j)=a(i,j)+b(i,j)*w
  1   continue
      return
      end
c
c**********************************************
cERN	NOT USED
      subroutine cumuvf(a, b, w, n)

      implicit none

      integer n, i
      complex*16 a(n), b(n)
      complex*16 w

      do 1 i=1, n
      a(i)=a(i)+b(i)*w
  1   continue
      return
      end
c
c**********************************************
cERN	NOT USED
      subroutine cumuv(a, b, n)

      implicit none

      integer n, i
      complex*16 a(n), b(n)
      do 1 i=1,n
      a(i)=a(i)+b(i)
  1   continue
      return
      end
c
c**********************************************
cERN	NOT USED
      subroutine mul3tf(r,a,b,c,fac,n)

      implicit none

      integer n
      double precision a(n,n),c(n,n)
      complex*16 b(n,n),r(n,n)
      
      integer i, j, k
      complex*16 t(12,12),t2(12,12),fac
      complex*16 czero
      data czero/(0.d0,0.d0)/

c     r:=(a transposed)*b*c*fac  (fac scalar)
      do 1 j=1,n
      do 1 i=1,n
      t(i,j)=czero
      t2(i,j)=czero
  1   continue
      do 2 i=1,n
      do 2 j=1,n
      do 2 k=1,n
      t(i,j)=t(i,j)+b(i,k)*c(k,j)
  2   continue
      do 3 i=1,n
      do 3 j=1,n
      do 3 k=1,n
      t2(i,j)=t2(i,j)+fac*a(k,i)*t(k,j)
  3   continue
      do 4 i=1,n
      do 4 j=1,n
      r(i,j)=t2(i,j)
  4   continue
      return
      end
c
c**********************************************
cERN	NOT USED
      subroutine mu3tf2(r,a,b,fac,m,n,nia,nja)

      implicit none

      integer m, n, nia, nja
      complex*16 b(nia,nia),r(nja,nja),t(12,12),t2(12,12),toto(288),fac
      equivalence (t(1,1),toto(1)),(t2(1,1),toto(12*12+1))
      double precision a(nia,nja)

c     n.b.: only columns 1 to n used in a
c     r:=(a transposed)*b*a*fac  (fac scalar)

      integer i, j, k
      complex*16 czero
      czero=(0.d0,0.d0)

      do 1 i=1,288
      toto(i)=czero
  1   continue

      do 2 i=1,m
      do 2 j=1,n
      do 2 k=1,m
      t(i,j)=t(i,j)+b(i,k)*a(k,j)
  2   continue
      do 3 i=1,n
      do 3 j=1,n
      do 3 k=1,m
      t2(i,j)=t2(i,j)+a(k,i)*t(k,j)
  3   continue
      do 4 i=1,n
      do 4 j=1,n
      r(i,j)=t2(i,j)*fac
  4   continue
      return
      end
c
c******************************************************
cERN  USED IN ASPLAS and ASCURL
      subroutine mu3tf2bl(r, al, a, b, nblo, icola)

      implicit none
      integer nblo, icola
      double precision a(6,icola)
      complex*16 b(6,6,nblo), r(12,12,nblo), al
c      double precision a(6,*)
c      complex*16 b(6,6,*), r(12,12,*), al

c     r:= r + al * trans(a) * b * a, where a is a real (6,icola) matrix
c     b is an array of nblo 6*6 blocks
      integer i, j, k, l
      integer n
      complex*16 czero, fac, fac2
      data czero/(0.d0,0.d0)/

      do 1 k = 1, 6
      do 1 i = 1, icola
      fac2 = al * a(k,i)
      do 1 l = 1, 6
      do 1 j = 1, icola
      fac = fac2 * a(l,j)
      do 1 n = 1, nblo
      r(i,j,n) = r(i,j,n) + fac * b(k,l,n)
cPL gives array bound error:
cPL      call zaxpy(nblo, fac, b(k,l,1), 36, r(i,j,1), 144)
  1   continue

      return
      end
c
c******************************************************
cERN  USED IN ASPLAS and ASCURL
      subroutine mu3tf2bl_fast(r, al, a, b, nblo, icola)

      implicit none
      integer nblo, icola
      double precision a(6,icola)
      complex*16 b(6,6,nblo), r(12,12,nblo), al
      integer n
c     r:= r + al * trans(a) * b * a, where a is a real (6,icola) matrix
c     b is an array of nblo 6*6 blocks
	do n = 1, nblo
	   r(:,:,n) = r(:,:,n) + 
     ;              al * matmul( transpose(a), matmul(b(:,:,n),a) )
	end do

      return
      end
c
c**********************************************
cERN  NOT USED
      subroutine mul3t(r,a,b,c,n)

      implicit none


      integer n
      complex*16 b(n,n),r(n,n),t(12,12),t2(12,12)
      double precision a(n,n),c(n,n)

c     r:=a transposed*b*c (n<=12)

      integer i, j, k
      complex*16 czero
      czero=(0.d0,0.d0)

      do 1 j=1,n
      do 1 i=1,n
      t(i,j)=czero
      t2(i,j)=czero
  1   continue

      do 2 i=1,n
      do 2 j=1,n
      do 2 k=1,n
      t(i,j)=t(i,j)+b(i,k)*c(k,j)
  2   continue

      do 3 i=1,n
      do 3 j=1,n
      do 3 k=1,n
      t2(i,j)=t2(i,j)+a(k,i)*t(k,j)
  3   continue

      do 4 i=1,n
      do 4 j=1,n
      r(i,j)=t2(i,j)
  4   continue
      return
      end
c
c**********************************************
cERN  NOT USED
      subroutine mul3t2(r,a,b,m,n,nia,nja)

      implicit none
      integer m, n, nia, nja

      integer i, j, k
      complex*16 b(nia,nia),r(nja,nja),t(12,12),t2(12,12),toto(288)
      equivalence (t,toto),(t2,toto(12*12+1))
      double precision a(nia,nja)

c     ! nia,nja are row and column dimensions declared in calling
c     program. m and n are the work dimensions for this routine.
c     m,n,nia,nja <=12.
c     r:=(a transposed) * b * a

      complex*16 czero
      czero=(0.d0,0.d0)

      do 1 i=1,288
      toto(i)=czero
  1   continue

      do 2 i=1,m
      do 2 j=1,n
      do 2 k=1,m
      t(i,j)=t(i,j)+b(i,k)*a(k,j)
  2   continue

      do 3 i=1,n
      do 3 j=1,n
      do 3 k=1,m
      t2(i,j)=t2(i,j)+a(k,i)*t(k,j)
  3   continue

      do 4 i=1,n
      do 4 j=1,n
      r(i,j)=t2(i,j)
  4   continue
      return
      end
c
c**********************************************
cERN  NOT USED
      subroutine mul2hc(r,a,b,n)

      implicit none


      integer n
      complex*16 b(n,n),r(n,n),t(12,12),a(n,n)

c     r:=(a transp.conj.)*b

      integer i, j, k
      complex*16 czero
      czero=      (0.d0,0.d0)

      do 1 j=1,n
      do 1 i=1,n
      t(i,j)=czero
  1   continue

      do 2 i=1,n
      do 2 j=1,n
      do 2 k=1,n
      t(i,j)=t(i,j)+conjg(a(k,i))*b(k,j)
  2   continue

      do 3 i=1,n
      do 3 j=1,n
      r(i,j)=t(i,j)
  3   continue

      return
      end
c
c**********************************************
cERN  USED in INTBOU
      subroutine mul2c(r,a,b,n)

      implicit none
      integer n
      complex*16 b(n,n),r(n,n),t(12,12),a(n,n)

c           r:=a*b

      integer i, j, k
      complex*16 czero
      czero=(0.d0,0.d0)

      do 1 j=1,n
      do 1 i=1,n
      t(i,j)=czero
  1   continue

      do 2 j=1,n
      do 2 k=1,n
      do 2 i=1,n
      t(i,j)=t(i,j)+a(i,k)*b(k,j)
  2   continue

      do 3 j=1,n
      do 3 i=1,n
      r(i,j)=t(i,j)
  3   continue

      return
      end
c
c*****************************************
cERN  USED in INTBOU
      subroutine mul2c_fast(r,a,b,n)

      implicit none
      integer n
      complex*16 b(n,n),r(n,n),a(n,n)

c           r:=a*b
	r = (0.d0,0.d0)
	r = matmul(a,b)

      return
      end
c
c*****************************************
cERN  USED in INTBOU
      subroutine mul2c2(r,a,b,ma,na,nb)

      implicit none
      integer ma, na, nb

      complex*16 r(ma,nb), a(ma,na), b(na,nb), t(12,12)

c           r:=a*b

      integer i, j, k
      complex*16 czero, temp
cPL      czero=(0.d0,0.d0)
	data czero/(0.d0,0.d0)/
	save czero

      do 1 j=1,nb
      do 1 i=1,ma
      t(i,j)=czero
  1   continue

      do 2 j=1,nb
      do 2 k=1,na
	temp = b(k,j)
      do 2 i=1,ma
      t(i,j)=t(i,j)+a(i,k)*temp
  2   continue

      do 3 j=1,nb
      do 3 i=1,ma
      r(i,j)=t(i,j)
  3   continue

      return
      end
c
c*****************************************
cERN  USED in INTBOU
      subroutine mul2c2_fast(r,a,b,ma,na,nb)

      implicit none
      integer ma, na, nb

      complex*16 r(ma,nb), a(ma,na), b(na,nb)
c           r:=a*b
	r = (0.d0,0.d0)
	r = matmul(a,b)

      return
      end
c
c*****************************************
cERN  USED in INTBOU
      subroutine mul2c3(r,a,b,ma,na,nb,ldr,lda,ldb)

      implicit none

      integer ma, na, nb, ldr, lda, ldb
      complex*16 r(ldr,nb),a(lda,na),b(ldb,nb),t(12,12)
c      complex*16 r(ldr,*),a(lda,*),b(ldb,*),t(12,12)

c           r:=a*b         (ma,na)*(na,nb), max 12*12

      integer i, j, k
      complex*16 czero
cPL      czero=(0.d0,0.d0)
	data czero/(0.d0,0.d0)/
	save czero

      do 1 i=1,ma
      do 1 j=1,nb
      t(i,j)=czero
  1   continue

      do 2 j=1,nb
      do 2 k=1,na
      do 2 i=1,ma
      t(i,j)=t(i,j)+a(i,k)*b(k,j)
  2   continue

      do 3 j=1,nb
      do 3 i=1,ma
      r(i,j)=t(i,j)
  3   continue

      return
      end
c
c*****************************************
cERN  USED in INTBOU
      subroutine mul2c3_fast(r,a,b,ma,na,nb,ldr,lda,ldb)

      implicit none

      integer ma, na, nb, ldr, lda, ldb
      complex*16 r(ldr,nb),a(lda,na),b(ldb,nb)

c     r = a*b         (ma,na)*(na,nb), max 12*12
	r = (0.d0,0.d0)
cPL25/8/04 BUG
c	r = matmul(a,b)
	r(1:ma,1:nb) = matmul(a(1:ma,1:na), b(1:na,1:nb))

      return
      end
c
c*****************************************
cERN  NOT USED
      subroutine add(r,a,b,n)

      implicit none
      integer n
      complex*16 r(n,n),a(n,n),b(n,n)

      integer i, j
      do 1 i=1,n
      do 1 j=1,n
      r(i,j)=a(i,j)+b(i,j)
  1   continue

      return
      end
c
c**************************************
cERN  NOT USED 
      double precision function cmodul(z,n,ipores)

      implicit none
      integer n, ipores
      complex*16 z(n),zdotc
      double precision t
      external zdotc

c     integer i
      t= dreal( zdotc(n, z, 1, z, 1) )
c     t = 0.d0
c     do 1 i=1,n
c     t=t+abs(z(i))**2
c 1   continue
c     ipores .ne.1: result is modulus squared:
      cmodul=t
c     ipores =1: result is modulus:
      if(ipores.eq.1)cmodul=dsqrt(t)
      return
      end
c
c******************************
cERN	USED in OUTRFF
      double precision function rveno(z,n,ipores)

      implicit none
      integer n, ipores
      double precision ddot
      double precision z(n),t

c     integer i
      
      t = ddot(n,z,1,z,1)
c     t=0.d0
c     do 1 i=1,n
c     t=t+z(i)**2
c 1   continue
c     ipores .ne.1: result is norm squared:

      rveno = t
c     ipores =1: result is norm:
      if(ipores.eq.1)rveno = dsqrt(t)
      return
      end
c
c******************************
cERN	USED in OUBIT, OUTPOW and OUTRFF
      double precision function trapeze(nabs,absi,f,incf,ipow)

      implicit none
      integer nabs, incf, ipow
      double precision absi(nabs)
      double precision f(*),t

c     trapezoidal integration of f(x)*x**ipow over [absi(1),absi(nabs)]
c     f tabulated at abscissae absi(i) and stored with increment incf:
c     f(1+(i-1)*incf) is value of f at absi(i).

      integer nabsf, if, i
      
      nabsf = 1 + (nabs-1)*incf

      if(ipow.eq.0)then
c     -----------------
      t = 0.5d0 * (
     ; f(1)*(absi(2)-absi(1)) + f(nabsf)*(absi(nabs)-absi(nabs-1))  )

      if = 1
      do 1 i = 2, nabs-1
      if = if + incf
      t = t + 0.5d0 * f(if) * ( absi(i+1) - absi(i-1) )
  1   continue
      else if(ipow.eq.1)then
c     ----------------------
      t = 0.5d0 * (
     ;   absi(1) * f(1) * (absi(2)-absi(1))
     ; + absi(nabs) * f(nabsf) * (absi(nabs)-absi(nabs-1))
     ;            )

      if = 1
      do 3 i = 2, nabs-1
      if = if + incf
      t = t + 0.5d0 * absi(i) * f(if) * ( absi(i+1) - absi(i-1) )
  3   continue
      else
c     ----
      t = 0.5d0 * (
     ;   absi(1)**ipow * f(1) * ( absi(2) - absi(1) )
     ; + absi(nabs)**ipow * f(nabsf) * ( absi(nabs) - absi(nabs-1) )
     ;            )

      if = 1
      do 2 i = 2, nabs-1
      if = if + incf
      t = t + 0.5d0 * absi(i)**ipow * f(if) * (absi(i+1)-absi(i-1))
  2   continue
      end if
c     ------

      trapeze=t
      return
      end
c
c*****************************************
cERN	USED in OUBIT, OUTPOW, OUTRFF and TABLES

      double precision function gauss(linter, absi, f, incf, ipow)

      implicit none
      integer incf, ipow
      double precision absi(*), linter
      double precision f(*)

      include 'pardim.copy'
      include 'comin2.copy'

      integer if, i
      double precision t
      

c     gauss integration of f(x)*x**ipow over [absi(1),absi(ngauss)]
c     f tabulated at gauss points and stored with increment incf:
c     f(1+(i-1)*incf) is value of f at absi(i).
c     linter is interval length.

      if = 1
      t = 0.d0
      if(ipow.eq.0)then
c     -----------------
      do 1 i = 1, ngauss
      t = t + wga(i) * f(if)
      if = if + incf
  1   continue
      else if(ipow.eq.1)then
c     ----------------------
      do 3 i = 1, ngauss
      t = t + wga(i) * f(if) * absi(i)
      if = if + incf
  3   continue
      else
c     ----
      do 2 i = 1, ngauss
      t = t + wga(i) * absi(i)**ipow * f(if)
      if = if + incf
  2   continue
      end if
c     ------

      gauss = linter * t
      return
      end
c
c*****************************************
cERN  USED in VACMBG2
      subroutine mul3(r,a,b,c,i1,i2)

      implicit none
      integer i1, i2
      complex*16 a(i1,i2),b(i2,i2),c(i2,i1),r(i1,i1)

c     r:=a*b*c , where r,a,b,c are complex matrices with
c                respective dimensions:
c     i1*i1 / i1*i2,i2*i2,i2*i1
c     this version: i1,i2< 12,  r,a,b,c complex*16 ||

      integer i, j, k
      complex*16 t(12,12),t2(12,12),cdzero
      cdzero= (0.d0,0.d0)

      do 3 j=1,12
      do 3 i=1,12
      t(i,j)=cdzero
      t2(i,j)=cdzero
  3   continue

      do 1 i=1,i2
      do 1 j=1,i1
      do 1 k=1,i2
      t(i,j)=t(i,j)+b(i,k)*c(k,j)
  1   continue

      do 2 i=1,i1
      do 2 j=1,i1
      do 2 k=1,i2
      t2(i,j)=t2(i,j)+a(i,k)*t(k,j)
  2   continue
      do 4 i=1,i1
      do 4 j=1,i1
      r(i,j)=t2(i,j)
  4   continue

      return
      end
c
c*******************************************
cERN  USED in VACMBG2
      subroutine mul3_fast(r,a,b,c,i1,i2)

      implicit none
      integer i1, i2
      complex*16 a(i1,i2),b(i2,i2),c(i2,i1),r(i1,i1)

c     r:=a*b*c , where r,a,b,c are complex matrices with
c                respective dimensions: i1*i1,i1*i2,i2*i2,i2*i1

	r = matmul(a, matmul(b,c))

      return
      end
c
c*******************************************
cERN  USED in ASPLAS (circ = true)
      subroutine mul3bl(a, b, c, nblo, tra, ntra)

      implicit none
      integer nblo, ntra
      complex*16 a(6,6), b(6,6,nblo), c(6,6), tra(6,6,ntra)

c     b:=a*b*c , where a,b,c are complex 6*6 matrices
c     b is an array of nblo blocks
c     tra is a (6,6,ntra) complex work area.
c     new version (26/9/91) for ntra possibly <nblo.

      integer ns, n, n1, nb, i, j, k, l
      complex*16 czero,fac,fac2
      data czero/(0.d0,0.d0)/

      if(nblo .le. 0)return
        if(ntra.le.0)then
        write(6,*)'mul3bl: ntra =',ntra,'; no work space available!'
        stop
        end if
      ns = nblo / ntra
      if(ns * ntra .ne. nblo)ns = ns + 1
      n = 1

      do n1 = 1, ns
c     -------------
      nb = min0(ntra, nblo-n+1)
      call zset(6*6*nb, czero, tra, 1)

        do i = 1, 6
          do l = 1, 6
          fac2 = a(i,l)
            do j = 1, 6
              do k = 1, 6
              fac = fac2 * c(k,j)
              call zaxpy(nb, fac, b(l,k,n), 36, tra(i,j,1), 36)
              end do
            end do
          end do
        end do
      call zcopy(6*6*nb, tra, 1, b(1,1,n), 1)
      n = n + nb

      end do
c     ------

      return
      end
c
c*******************************************
cERN  USED in ASPLAS (circ = true)
      subroutine mul3bl_fast(a, b, c, nblo)
cPL      subroutine mul3bl_fast(a, b, c, nblo, tra, ntra)

      implicit none
      integer nblo
cPL      integer nblo, ntra
      complex*16 a(6,6), b(6,6,nblo), c(6,6)
cPL      complex*16 a(6,6), b(6,6,nblo), c(6,6), tra(6,6,ntra)

c     b:=a*b*c , where a,b,c are complex 6*6 matrices
c     b is an array of nblo blocks

cPL25/8/04 Buggy section!!
c      integer ns, n, n1, nb
c
c      n = 1
c
c      do n1 = 1, ns
cc     -------------
c         nb = min0(ntra, nblo-n+1)
c	   b(:,:,n) = matmul(a, matmul(b(:,:,n),c))
c	   n = n + nb
c
c      end do
cc     ------
      integer n

      do n = 1, nblo
	b(:,:,n) = matmul(a, matmul(b(:,:,n),c))
      end do

      return
      end
c
c*******************************************
cERN	USED in OUTRF2
      subroutine ccrosp(a, b, c)

      implicit none
      complex*16 a(3), b(3), c(3)

c     conjugate cross product of 2 complex vectors (3 components)
c     a := conj(b) ^ c

      a(1) = conjg( b(2) ) * c(3) - conjg( b(3) ) * c(2)
      a(2) = conjg( b(3) ) * c(1) - conjg( b(1) ) * c(3)
      a(3) = conjg( b(1) ) * c(2) - conjg( b(2) ) * c(1)
      return
      end
c
c*****************************************************
cERN	NOT USED
      subroutine ccvema(v, a, b, n)

      implicit none
      integer n
      complex*16 v(n), a(n), b(n,n)

c     v := v + conjtrans(a)*b
      integer i, j
      complex*16 tra

      do 1 i = 1, n
      tra = conjg(a(i))
      do 1 j = 1, n
      v(j) = v(j) + tra * b(i,j)
   1  continue

      return
      end
c
c************************************
c
