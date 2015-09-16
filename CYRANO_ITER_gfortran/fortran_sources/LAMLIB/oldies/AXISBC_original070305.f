      subroutine axisbc(c, ldc, nbc, nc)

ccc      USE numerical_libraries  ! PL (for IMSL compilation checks)
      implicit none
      
      integer ldc, nbc, nc
      
      complex*16 c(ldc,*)

c     Builds matrix of inner boundary constraints (magn. axis or metal wall)
c
c     Input:
c     LDC, leading dimension of matrix C.
c 
c     Output:
c     C(nbc,nc) the matrix of constraints
c     NBC, its number of rows, the number of boundary conditions
c     NC, its number of columns
c 
c     Current version: no plasma FLR effects

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'comant.copy'
      include 'commod.copy'
      include 'comsub.copy'
      include 'comreg.copy'
      include 'comswe.copy'
      include 'comfou.copy'
      include 'comrot.copy'
      include 'comro2.copy'
      include 'commag.copy'
      include 'compla.copy'
      include 'comfin.copy'
      include 'comfic.copy'
      include 'comin2.copy'
      include 'comequ.copy' 
      include 'comphy.copy'
      include 'comfft.copy'

      character*3 eletyp

      integer icolo, icpblo, idllo, ibulo, i, j, jd, m2
c     ;, k1, k2

c      double precision al, be
      
      double precision polst, bca(npfft+1), bcai(npfft+1), t1, cr(npfft+1,2)
      complex*16 foutra(-npft2:npft2,2)
           
      eletyp = styp(1)
      icolo = iconn(1)
      ibulo = ibub(1)
      idllo = iddl(1)
      icpblo = icolo + ibulo

c     Number of columns concerned by the boundary conditions:
      nc = icolo * nmode(0)
c      nc = lblock(1)
      if(.not.crown .and. eletyp.eq.'M23' .and. nabcap)
     ;nc = lblock(1) + icolo * nmode(1)
c     ;nc = lblock(1) + lblock(2)
      call zset(nc*ldc, czero, c, 1)

      if(crown)then
c     -------------
	print *, 'B.C. at inner metal shell'
	write(nofile,*) 'B.C. at inner metal shell'
      i = 0
c     Metallic conditions at inner wall
        if(eletyp.eq.'HEC')then
c       E+ - E- = 0; E// = 0:
          do mr = 1, nmode(0)
          jd = (mr-1) * icolo
          i = i + 1
          c(i,jd+1) = cun
          c(i,jd+3) = - cun
          i = i + 1
          c(i,jd+5) = cun
          end do
        else if(eletyp.eq.'M23')then
c       Etheta = 0; Ephi = 0:
          do mr = 1, nmode(0)
          jd = (mr-1) * icolo
          i = i + 1
          c(i,jd+2) = cun
          i = i + 1
          c(i,jd+4) = cun
          end do
        else if(eletyp.eq.'CAR')then
c       Etheta = 0; Ephi = 0:
        print *, 'New CAR element type: must implement convolution for metallic condition in AXISBC'
	  stop
        end if

      else  ! Magnetic axis:
c     ----
        if(circ .or. (dshape.and.kappa.eq.1.d0) .or. eletyp.eq.'CAR')then
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       In circular concentric geometry, thetabar = theta, this deals with const. k// coordinates as well
	  print *, 'CIRCULAR B.C. at axis'
        i = 0
c       NB: in D shape, elongation is the only parameter influencing axis b.c.
          if(eletyp.eq.'HEC')then
c           Essential conditions: thesis (4.23)
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            jd = (mr-1) * icolo
              if(m.ne.-1)then
              i = i + 1
              c(i,jd+1) = cun
              end if
              if(m.ne.1)then
              i = i + 1
              c(i,jd+3) = cun
              end if
              if(m.ne.0)then
              i = i + 1
              c(i,jd+5) = cun
              end if
            end do
            
            if(nabcap)then
c           Natural conditions: thesis (4.28)
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            m2 = m / 2
            jd = (mr-1) * icolo
              if(2*m2.ne.m)then
              i = i + 1
              c(i,jd+2) = cun
              i = i + 1
              c(i,jd+4) = cun
              else if(m.ne.0)then
              i = i + 1
              c(i,jd+2) = dfloat(m2) + 1.d0
              c(i,jd+4) = dfloat(m2) - 1.d0
              end if
              if(iabs(m).ne.1)then
              i = i + 1
              c(i,jd+6) = cun
              end if
            end do
            end if

          else if(eletyp.eq.'M23')then
c           Essential conditions: thesis (4.22)
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            jd = (mr-1) * icolo
              if(iabs(m).ne.1)then
              i = i + 1
c             Erho = 0:
              c(i,jd+1) = cun
              i = i + 1
c             Etheta = 0:
              c(i,jd+2) = cun
              else
              i = i + 1
c             i*m*Erho - Etheta = 0:
              c(i,jd+1) = dcmplx(0.d0,dfloat(m))
              c(i,jd+2) = - cun              
              end if
              if(m.ne.0)then
              i = i + 1
c             Ephi = 0:
              c(i,jd+4) = cun
              end if
            end do

            if(nabcap)then
c           Case natural conditions imposed a priori: thesis (4.27)
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            m2 = m / 2
            jd = (mr-1) * icolo
              if(2*m2.ne.m)then
c  Couldn't work with former treatment of first block (this constraint involves Erho at right node):
              i = i + 1
c             Odd modes: dErho/drho = 0. The following 3 lines express this derivative for the quadratic basis functions.
              c(i,jd+1) = dcmplx(-3.d0,0.d0)
              c(i,icolo*nmode(0)+(mr-1)*ibulo+1) = dcmplx(4.d0,0.d0)
              c(i,icpblo*nmode(0)+jd+1) = dcmplx(-1.d0,0.d0)
              i = i + 1
c             Odd modes: dEtheta/drho = 0
              c(i,jd+3) = cun
              else if(m.ne.0)then
c  Couldn't work with former treatment of first block (this constraint involves Erho at right node):
              i = i + 1
c             Even nonzero modes: (im/2)*dErho/drho - dEtheta/drho = 0
              c(i,jd+1) = dcmplx(0.d0,-3.d0*dfloat(m2))
              c(i,icolo*nmode(0)+(mr-1)*ibulo+1) = dcmplx(0.d0,4.*dfloat(m2))
              c(i,jd+1+icpblo*nmode(0)) = dcmplx(0.d0,-dfloat(m2))
              c(i,jd+3) = dcmplx(-1.d0,0.d0)
              end if
              if(iabs(m).ne.1)then
              i = i + 1
c             All modes except +1 and -1: dEphi/drho = 0
              c(i,jd+5) = cun
              end if
            end do
            end if

          else if(eletyp.eq.'CAR')then
c           Essential conditions: ER, EY, Ephi = 0 for m nonzero
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            jd = (mr-1) * icolo
              if(m.ne.0)then
              i = i + 1
              c(i,jd+1) = cun
              i = i + 1
              c(i,jd+3) = cun
              i = i + 1
              c(i,jd+5) = cun
              end if
            end do
            
            if(nabcap)then
c           Natural conditions: even modes have zero radial derivative
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            m2 = m / 2
            jd = (mr-1) * icolo
              if(2*m2.eq.m)then
              i = i + 1
              c(i,jd+2) = cun
              i = i + 1
              c(i,jd+4) = cun
              i = i + 1
              c(i,jd+6) = cun
	        end if
            end do
            end if

          end if

        else  ! All noncircular cases (not using Cartesian coordinates)
c       ~~~~

c Assumes poloidal mode spectrum includes -1,0,1.
c a completer pour cond. nat. sur H quand nabcap = .t.
c        al = 0.5d0 * (kappa + 1.d0 / kappa)
c        be = 0.25d0 * (kappa - 1.d0 / kappa)

c@     NEW (May 2000):
      polst = twopi / dfloat(npfft)
c     bca receives poloidal integral of coefficient 'newmu' at axis:
c     NB: i) in D-shape, newmu = 1/(lambda thesis eq.4.30)
c        ii) the following approach is valid for arbitrary geometry

      bca(1) = 0.d0
        do i = 2, npfft+1
        bca(i) = bca(i-1) + 0.5 * (eqt(1,i-1,10) + eqt(1,i,10)) *  polst
        end do
          if(cokpco)then
c         Interpolate bca onto uniform Thetabar grid:
          call interp2(ckt(1,1,3), nabplo, bca, 1, npfft+1, polang, bcai, npfft+1)
            do i = 1, npfft+1
            bca(i) = bcai(i)
            end do
	    end if
c      write(nofile,*)'Angle psi for axis boundary conditions: i,theta,psi:'
c        do i = 1, npfft + 1
c        write(nofile,*)i, polang(i), bca(i)
c        end do
      
c     Poloidal FFT of cos and sin(bca):
c     Normalization for RCFFT2: 
c      t1 = 1.d0 / dfloat(2 * npfft)
c     Normalization for IMSL: 
      t1 = 1.d0 / dfloat(npfft)
        do i = 1, npfft
        cr(i,1) = t1 * dcos(bca(i))
        cr(i,2) = t1 * dsin(bca(i))
        end do

        do j = 1, 2
c       fft,n=2**m case: voir cas sym, antisym;
c                                                ci-dessous general en reel.
cccccccccccccccccccccccc ERNESTO cccccccccccccccccccccccccccc

c     WARNING: For FFT of REAL functions, the index of the Fourier coefs. 
c     in CRAY and IMSL routines are different:
cERN    CALL RCFFT2(0, -1, NPFFT, CR(1,j), WORK, foutra(0,j))
	  call df2trf(npfft, cr(1,j), cr(1,j), work2) ! IMSL: forward exp(-i...)

c       Storage of IMSL transform in complex array:
c     NB here we want integral of exp(-i k theta), so no complex conjugate!
cPL25Nov04 CHECK THIS!
	  foutra(0,j) = dcmplx(cr(1,j), 0.d0)
	    do i = 1, npft2-1
	    foutra(i,j) = dcmplx(cr(2*i,j), cr(2*i+1,j))
          end do
	  foutra(npft2,j) = dcmplx(cr(npfft,j), 0.d0)
          do i = 1, npft2
          foutra(-i,j) = dconjg(foutra(i,j))
          end do       
	  end do  ! j

c      do i = 0,5
c         print *, foutra(i,1)
c	end do
c	print *, 'STOP in AXISBC Line 250'
c	stop     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write(nofile,*)'Fourier coeffs. for axis boundary conditions: j,cos,sin'
        do i = 0, npft2
        write(nofile,1000)i, foutra(i,1), foutra(i,2)
        end do
 1000 format(1h , i4, 2x, 2(1h(, g15.7, 1h,, g15.7, 1h), 2x))
c
        i = 0
          if(eletyp.eq.'HEC')then
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            m2 = m / 2
            jd = (mr - 1) * icolo
c             D-shape: we know even modes of E+, E- all vanish at axis:
              if(.not.cokpco .and. (2 * m2 .eq. m))then
              i = i + 1
              c(i,jd+1) = cun
              i = i + 1
              c(i,jd+3) = cun
              else
c             Equations for odd modes, or for all in constant k// coordinates:
c              k1 = max0(-ncrot, minf(0)-m)
c              k2 = min0(ncrot, msup(0)-m)
c             E+: impose mode amplitude relative to m=-1:
                if(m .ne. -1)then
                i = i + 1
C               Coeff.of E+(m):
                c(i,jd+1) = foutra(-1,1) - ci * foutra(-1,2)
C               Coeff.of E+(-1):
                j = (-minf(0) - 1) * icolo
                c(i,j+1) = -(foutra(m,1) - ci * foutra(m,2))
C Old:
C                C(i,jd+1) = dcmplx(m*al+1.d0, 0.d0)
C                  if(k1.le.-2)then
C                  j = jd - 2 * icolo
C                  C(i,j+1) = dcmplx(be*(m-2), 0.d0)
C                  end if
C                  if(k2.ge.2)then
C                  j = jd + 2 * icolo
C                  C(i,j+1) = dcmplx(be*(m+2), 0.d0)
C                  end if
                end if
c               E-: impose mode amplitude relative to m=+1:
                if(m .ne. 1)then
                i = i + 1
                j = (-minf(0) + 1) * icolo
C               Coeff.of E-(m):
                c(i,jd+3) = foutra(1,1) + ci * foutra(1,2)
C               Coeff.of E-(+1):
                c(i,j+3) = -(foutra(m,1) + ci * foutra(m,2))
C Old:
C                C(i,jd+3) = dcmplx(m*al-1.d0, 0.d0)
C                  if(k1.le.-2)then
C                  j = jd - 2 * icolo
C                  C(i,j+3) = dcmplx(be*(m-2), 0.d0)
C                  end if
C                  if(k2.ge.2)then
C                  j = jd + 2 * icolo
C                  C(i,j+3) = dcmplx(be*(m+2), 0.d0)
C                  end if
                end if
              end if
              
              if(m.ne.0)then
c             E//m = 0
              i = i + 1
              c(i,jd+5) = cun
              end if
            end do

            if(nabcap)then
c           Thesis (4.26) in +-//: E is C1; not yet implemented H is C0
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            m2 = m / 2
            jd = (mr - 1) * icolo
              if(2 * m2 .ne. m)then
              i = i + 1
              c(i,jd+2) = cun
              i = i + 1
              c(i,jd+4) = cun
              else 
              i = i + 1
              c(i,jd+6) = cun
              end if
            end do
            end if

          else if(eletyp.eq.'M23')then
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            jd = (mr - 1) * icolo
            m2 = m / 2
              if(.not.cokpco .and. (2 * m2 .eq. m))then
c             D-shape: 
c             Theta-phi coordinates: we know even modes of Erho, Etheta
c             all vanish at axis:
              i = i + 1
              c(i,jd+1) = cun
              i = i + 1
              c(i,jd+2) = cun
              else
c             Equations for odd modes, or all modes in constant k// coordinates:
c              k1 = max0(-ncrot, minf(0)-m)
c              k2 = min0(ncrot, msup(0)-m)
C              i = i + 1
c               E+: impose mode amplitude relative to m=-1:
                if(m .ne. -1)then
                i = i + 1
                j = (-minf(0) - 1) * icolo
C               Coeff.of Erho(m):
                c(i,jd+1) = foutra(-1,1) - ci * foutra(-1,2)
C               Coeff.of Etheta(m):
                c(i,jd+2) = ci * c(i,jd+1)
C               Coeff.of Erho(-1):
                c(i,j+1) = -(foutra(m,1) - ci * foutra(m,2))
C               Coeff.of Etheta(-1):
                c(i,j+2) = ci * c(i,j+1)
                end if
c               E-: impose mode amplitude relative to m=+1:
                if(m .ne. 1)then
                i = i + 1
                j = (-minf(0) + 1) * icolo
C               Coeff.of Erho(m):
                c(i,jd+1) = foutra(1,1) + ci * foutra(1,2)
C               Coeff.of Etheta(m):
                c(i,jd+2) = - ci * c(i,jd+1)
C               Coeff.of Erho(+1):
                c(i,j+1) = -(foutra(m,1) + ci * foutra(m,2))
C               Coeff.of Etheta(+1):
                c(i,j+2) = - ci * c(i,j+1)
                end if
C Old:
C              C(i,jd+1) = cun
C              C(i,jd+2) = dcmplx(0.d0, m*al)
C                if(k1.le.-2)then
C                j = jd - 2 * icolo
C                C(i,j+2) = dcmplx(0.d0,(m-2)*be)
C                end if
C                if(k2.ge.2)then
C                j = jd + 2 * icolo
C                C(i,j+2) = dcmplx(0.d0,(m+2)*be)
C                end if
C                if(iabs(m) .ne. 1)then
C                i = i + 1
C                C(i,jd+1) = dcmplx(0.d0, -m*al)
C                C(i,jd+2) = cun
C                  if(k1.le.-2)then
C                  j = jd - 2 * icolo
C                  C(i,j+1) = dcmplx(0.d0,-(m-2)*be)
C                  end if
C                  if(k2.ge.2)then
C                  j = jd + 2 * icolo
C                  C(i,j+1) = dcmplx(0.d0,-(m+2)*be)
C                  end if
C                end if
              end if

              if(m .ne. 0)then
              i = i + 1
              c(i,jd+4) = cun
              end if
            end do

            if(nabcap)then
c           Thesis (4.26): E is C1; not yet implemented H is C0
            do mr = 1, nmode(0)
            m = minf(0) + mr - 1
            m2 = m / 2
            jd = (mr - 1) * icolo
              if(2 * m2 .ne. m)then
c             Erho'=0 (m odd):
              i = i + 1
              c(i,jd+1) = dcmplx(-3.d0, 0.d0)
              c(i,icolo*nmode(0)+(mr-1)*ibulo+1) = dcmplx(4.d0, 0.d0)
              c(i,icpblo*nmode(0)+jd+1) = dcmplx(-1.d0, 0.d0)
c             Etheta'=0 (m odd):
              i = i + 1
              c(i,jd+3) = cun
              else 
c             Ephi'=0 (m even):
              i = i + 1
              c(i,jd+5) = cun
              end if
            end do
            end if

          end if
c2004        else
c2004       ~~~~
c2004        write(nofile,*)'axisbc: not yet for geneq=.true. I stop.'
c2004	  stop
        end if
c       ~~~~~~

      end if
c     ------

      nbc = i
      return
      end
