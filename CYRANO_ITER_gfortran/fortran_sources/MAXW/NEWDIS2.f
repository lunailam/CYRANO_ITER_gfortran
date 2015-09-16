      subroutine disper(ipath, krho2, nrootp, wridt)

      implicit none
      integer ipath, nrootp
      logical wridt
     
c     Studies local dispersion relation
c     ipath=1: compute approx. starting values for roots
c              otherwise use previous values in krho2
c     nrootp: number of roots to search
c     Species with a general equilibrium distribution are considered Maxwellian
c     here!
c     imeth = 0: old numerical root search
c           = 1: solve quadratic or cubic equ. for kperp**2
c     wridt: true to write some dielectric coefficients to 
c            unit 81 ('Tensor elements'), otherwise false
cERN	Modified in July 2004 to include the case 
cERN	of constant k// coordinates 

      include 'pardim.copy'
      include 'comgeo.copy'
      include 'commod.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'complo.copy'
      include 'complp.copy'
      include 'comnud.copy'
      include 'comswe.copy'
      include 'comphy.copy'

      complex*16 krho2(3,maxpom)
 
      integer i, idum, nkn, nnew, ifail
      
      real*8 kper2, x(2), fvec(2), wa(20)
     ;, npar, majR
      
      complex*16 
     ;  krv, krf, krs, krb
     ;, kperp2(3,maxpom)
     ;, sigsd(3,3), sigod(3,3)
     ;, detphy
     ;, aa(4), cro(3), akpe2(3), akr2(3), a(3), p3, q3, r3, s3
     ;, bp, bm, c, apl, ami, pp, d, pl, pr, pu, pv, del2, del
     
      external fcn, detphy
c      intrinsic dmax1
 
      errrel = 1.d-8
      itmax = 30
 
      if(cyl)then
      kphilo = kphi
      else
      kphilo = kphi * r0orta(intab,intabp)
      end if

      yinv = abscni(intab)
      rhoinv = yinv * rnori
      co = eqt(intab,intabp,14)
      si = eqt(intab,intabp,15)
	majR = eqt(intab,intabp,1)  ! ERN: major radius R(rho,theta)
      
         do 3 m = mstud1, mstud2, mstust
c        -------------------------------
      mr = (m - mstud1) / mstust + 1
      kthe = m * rhoinv
cERN	Include case of constant k// coordinates
	if(cokpco) then
		kpar = hachi(intab) * (m + n*qfactor(intab) ) 	! k//=const. value
	else
		kpar = kthe * si + kphilo * co  ! original value
	end if
      npar = kpar / k0
      npar2 = npar * npar
    
cERN	Include case of constant k// coordinates
	if(cokpco) then
		kper = co/dmax1(si,1.d-8) * (kpar - n/majR/co)
	else
	    kper = kthe * co - kphilo * si
      end if
	kper2 = kper * kper
      krv = k02 - kthe ** 2 - kphilo ** 2
 
      if(vacuum(ireg))then
c     --------------------
      krho2(1,mr) = krv
 
      else
c     ----
c     stupid and vital for HOTTEN:
c     INTAB, INTABP must be passed by COMSWE!
c      write(6,*)'calling hotten from disper'
      call hotten(kpar, sigsd, sigod, .false., idum, .true., .false.)
 
      lstix = sigsd(1,1) + cun
      rstix = sigsd(2,2) + cun
      sstix = 0.5 * (rstix + lstix)
      dstix = 0.5 * (rstix - lstix)
      pstix = sigsd(3,3) + cun
 
      if(flrops(ireg))then
      ghot11 = 0.5d0 + k0rn2 * sigod(1,1)
      ghot22 = 0.5d0 + k0rn2 * sigod(2,2)
      ghotpm = cun - 2.d0 * k0rn2 * sigod(1,2)
      ghot13 = ci * sqrt2 * k02 * rnorm * sigod(1,3)
      ghot23 = ci * sqrt2 * k02 * rnorm * sigod(2,3)
      ghot33 = cun + k0rn2 * sigod(3,3)
      else
      ghot11 = 0.5d0
      ghot22 = 0.5d0
      ghotpm = cun
      ghot13 = czero
      ghot23 = czero
      ghot33 = cun
      end if
      
      if(wridt)write(81,201)eqt(intab,intabp,1), dfloat(m)
     ;, dreal(sstix), dimag(sstix), dreal(dstix), dimag(dstix)
     ;, dreal(lstix), dimag(lstix), dreal(rstix), dimag(rstix)
     ;, dreal(pstix), dimag(pstix)

c      if(imeth.eq.0)then
c     ------------------
c     write(6,*)'disper:'
c     write(6,*)'ktheta= ',kthe,' ;  kphi= ',kphi
c     write(6,200)lstix,rstix,pstix,ghot11,ghot22,ghotpm
c    ;,ghot13,ghot23,ghot33
c 200 format(1h ,(4(g15.7,2x)))
 
c     local physical dispersion:
 
c     matnor(1) = k0 * sqrt(abs(lstix))
c     matnor(2) = k0 * sqrt(abs(rstix))
c     matnor(3) = k0 * sqrt(abs(pstix))
 
c     initial guess = 2 or 3 approx. plasma roots (k**2):
c     search kperp2; deduce krho2 afterwards
      krf = k02 * (npar2 - rstix) * (npar2 - lstix) / (sstix - npar2)
      krs = k02 * pstix * (1.d0-npar2/sstix)
      nguess = 2
 
      detnor=1.d0
 
      if( nrootp .gt. 2 )then
      krb = k02 * 2. * (npar2-sstix)/(0.5d0-ghot11)
      nguess = 3
      end if
 
      if(ipath.eq.1)then
      kperp2(1,mr) = krf
      kperp2(2,mr) = krs
      if(nrootp.gt.2)kperp2(3,mr) = krb
      end if
 
c      write(6,*)'starting root search in newdis2'
c      I = 0
c   1  CONTINUE
c      I = I + 1
c      ERRABS = 1.D-3
CC    ERRABS = 1.D-8 * ABS(DETPHY(1.001 * KPERP2(I,MR)))
CC    WRITE(6,*)'ERRABS= ',ERRABS,' ERRREL= ',ERRREL
c      NKN = I - 1
c      NGUESS = 1
c      NNEW = 1
CC     NKN = 0
CC     NGUESS = NROOTP
CC     NNEW = NROOTP
Cc      CALL ZANLY(DETPHY,ERRABS,ERRREL,NKN,NNEW,NGUESS,KPERP2(1,MR),
Cc     ; ITMAX,KPERP2(1,MR),INFER)

c     Initial guesses:
c      x(1) = dreal(kperp2(i,mr))
c      x(2) = dimag(kperp2(i,mr))
c      ifail = 0
c      call c05nbf(fcn, 2, x, fvec, errrel, wa, 20, ifail)
c      if(ifail.ne.0)write(6,*)'newdis: c05nbf ifail=',ifail
c      kperp2(i,mr) = dcmplx(x(1),x(2))
c      KRHO2(I,MR) = KPERP2(I,MR) - KPER2
c      IF(I.LT.NROOTP)GOTO 1
C     WRITE(6,*)'ITER:',(INFER(I),I=1,NROOTP)
c      WRITE(6,*)'roots by iteration: ',(KRHO2(I,MR),I=1,NROOTP)

c      else if(imeth.eq.1)then
c     -----------------------
        if(nrootp.lt.3)then
c     0th order flr: quadratic equation for nperp**2
      pl = lstix - npar2
      pr = rstix - npar2
      aa(1) = pstix * pl * pr / sstix
      aa(2) = npar2 - (pstix * (sstix - npar2) + rstix * lstix) / sstix
      aa(3) = cun
      aa(4) = czero
      del2 = aa(2)**2 - 4. * aa(1)
c     Manage to have first root the one with smallest real part:
      del = cdsqrt(del2) * dsign(1.d0,dreal(aa(2)))
      cro(1) = (- aa(2) + del) * 0.5
      cro(2) = (- aa(2) - del) * 0.5
      
        else
c     2nd order flr: cubic equation for nperp**2
      bp = - ghot11
      bm = - ghot22
      c = ghotpm
      apl = - ghot13 - npar
      ami =   ghot23 + npar
      pp = pstix
      d = - ghot33
      pl = lstix - npar2
      pr = rstix - npar2
      pu = bp * bm - 0.25 * c * c
      pv = bp * pr + bm * pl

      aa(1) = pp * pl * pr
      aa(2) = pp * pv + d * pl * pr - 0.5*(pl*ami**2+pr*apl**2)
      aa(3) = pp * pu + d * pv - 0.5*(c*apl*ami+bp*ami**2+bm*apl**2)
      aa(4) = d * pu
      cro(3) = 1.d40

        if(aa(1).eq.czero)then
c       One root is zero, remains quadratic equation:
        cro(3) = czero
        aa(1) = aa(2)
        aa(2) = aa(3)
        aa(3) = aa(4)
        aa(4) = czero
        end if

        if(aa(4).ne.czero)then
        a(1) = aa(1) / aa(4)
        a(2) = aa(2) / aa(4)
        a(3) = aa(3) / aa(4)

        call cubeq(a(3), a(2), a(1), cro)

c       Finds one root of cubic, then deflation and solve quadratic 
c       (based on Harwell routine)
c        a(1) = aa(1) / aa(4)
c        a(2) = aa(2) / (3.d0*aa(4))
c        a(3) = aa(3) / (3.d0*aa(4))
c        p3 = a(2) - a(3)*a(3)
c        q3 = 0.5 * (a(3)*(a(2)+2.*p3)-a(1))
c        s3 = 0.5 * cdsqrt(4.*p3**3+q3*q3)
c        r3 = (q3 + s3)**(1./3.)
c        cro(3) = r3 - p3 / r3 - a(3)
c       
c        aa(3) = cun
c        aa(1) = - a(1) / cro(3)
c        aa(2) = (aa(1) - 3. * a(2)) / cro(3)
c        end if

        else
c       Quadratic cases

        if(aa(3).eq.czero)then
        cro(2) = 1.d40
          if(aa(2).eq.czero)then
          cro(1) = 1.d40
          else
          cro(1) = - aa(1) / aa(2)
          end if
        else
        aa(1) = aa(1) / aa(3)
        aa(2) = aa(2) / aa(3)
        del2 = aa(2) * aa(2) - 4. * aa(1)
        del = cdsqrt(del2) * dsign(1.d0,dreal(aa(2)))
        cro(1) = (- aa(2) + del) * 0.5
        cro(2) = (- aa(2) - del) * 0.5
        end if

        end if

        end if

        do i = 1, nrootp
c        akpe2(i) = k02 * cro(i)
c        akr2(I) = akpe2(i) - KPER2
        kperp2(i,mr) = k02 * cro(i)
        krho2(i,mr) = kperp2(i,mr) - KPER2
		end do
c      WRITE(6,*)'analytical roots  : ',(akr2(I),I=1,NROOTP)

c      end if
c     ------
      END IF
C     ------
   3  CONTINUE
C     -------
 
C200  FORMAT(1H ,3(2(E16.8,1X)/1H ),4(4(E16.8,1X)/1H ))
 201  FORMAT(1H ,12(E12.4,1x))
C 201  FORMAT(1H ,E12.4,1x,I4,1x,10(E12.4,1x))
      RETURN
      END
C
C*************
C
      subroutine fcn(neq, x, fvec, iflag)

      IMPLICIT NONE
      integer neq, iflag
      double precision x(2), fvec(2)
c
c     auxiliary subroutine to solve dispersion with Nag routine c05nbf
c     similar to complex function detphy, but now solution of a real nonlinear
c     system of two equations.
c     on entry, x(1) and x(2) are the real and imaginary parts of kperp**2
c

      include 'pardim.copy'
      include 'comnud.copy'
      include 'commod.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comswe.copy'
      include 'comphy.copy'
      INTEGER IPVT(3), LDAP, LDFAC, NMAT, i, j, info

      DOUBLE PRECISION sidet

      COMPLEX*16 
     ;  kperp2, detphy
     ;, KRHO, KRHO2, KPL, KMI, TRA, KP2, NPER2
     ;, MAT(3,3), FAC(3,3)
 
      kperp2 = dcmplx(x(1),x(2))
      NPER2 = KPERP2 / K02
 
      IF(VACUUM(IREG))THEN
      DETPHY = K02 - KPAR**2

      ELSE IF(.NOT.FLROPS(IREG))THEN
      DETPHY = NPER2 * (NPER2 - (PSTIX * (1.d0 - NPAR2 / SSTIX) - NPAR2
     ; + RSTIX * LSTIX / SSTIX))
     ; + PSTIX * (NPAR2 - LSTIX) * (NPAR2 - RSTIX) / SSTIX
c     ; + PSTIX * ((NPAR2 * NPAR2 + RSTIX * LSTIX)/SSTIX -2. * NPAR2)

      ELSE
c     RESTE A FAIRE!!! FAUX MAINTENANT:
c      KPL = (KRHO + CI * KPER) * SQRT2I
c      KMI = (KRHO - CI * KPER) * SQRT2I
      KPL = cdsqrt(kperp2) * SQRT2I
      KMI = dconjg(kpl)
c      TRA = KRHO2 + KPER**2
      KP2 = KPAR**2
 
      MAT(1,1) = K02 * LSTIX - kperp2 * GHOT11 - KP2
      MAT(1,2) = GHOTPM * KPL**2
      MAT(1,3) = KPL * (KPAR + GHOT13)
      MAT(2,1) = GHOTPM * KMI**2
      MAT(2,2) = K02 * RSTIX - kperp2 * GHOT22 - KP2
      MAT(2,3) = KMI * (KPAR + GHOT23)
      MAT(3,1) = KMI * (KPAR + GHOT13)
      MAT(3,2) = KPL * (KPAR + GHOT23)
      MAT(3,3) = K02 * PSTIX - kperp2 * GHOT33
 
      LDAP = 3
      LDFAC = 3
      NMAT = 3
C     WRITE(6,*)MAT
      do j = 1, nmat
      do i = 1, nmat
      fac(i,j) = mat(i,j)
      end do
      end do
      call zgetrf(nmat, nmat, fac, ldfac, IPVT, info)

      sidet = 1.
      if(ipvt(1).ne.1)sidet = -sidet
      if(ipvt(2).ne.2)sidet = -sidet
      if(ipvt(3).ne.3)sidet = -sidet
      detphy = detnor * sidet * fac(1,1) * fac(2,2) * fac(3,3)
 
      END IF
      fvec(1) = dreal(detphy)
      fvec(2) = dimag(detphy)
 
      RETURN
      END
C
C************* to disappear: *******************************
C
      COMPLEX*16 FUNCTION DETPHY(KPERP2)

      IMPLICIT NONE
      COMPLEX*16 KPERP2

      include 'pardim.copy'
      include 'comnud.copy'
      include 'commod.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comswe.copy'
      include 'comphy.copy'

      INTEGER IPVT(3), LDAP, LDFAC, NMAT
     ;, I, J, info
      DOUBLE PRECISION sidet
      COMPLEX*16 KRHO, KRHO2, KPL, KMI, TRA, KP2, NPER2
      COMPLEX*16 MAT(3,3), FAC(3,3)
 
      NPER2 = KPERP2/K02
 
      IF(VACUUM(IREG))THEN
      DETPHY = K02 - KPAR**2
      ELSE IF(.NOT.FLROPS(IREG))THEN
      DETPHY = NPER2*(NPER2 - (PSTIX*(1.-NPAR2/SSTIX) - NPAR2
     ; + RSTIX*LSTIX/SSTIX))
     ; + PSTIX*((NPAR2*NPAR2 + RSTIX*LSTIX)/SSTIX -2.*NPAR2)
      ELSE
C     RESTE A FAIRE!!! FAUX MAINTENANT:
      KPL=(KRHO + CI*KPER)*SQRT2I
      KMI=(KRHO-CI*KPER)*SQRT2I
      TRA = KRHO2 + KPER**2
      KP2=KPAR**2
 
      MAT(1,1)= K02*LSTIX-TRA*GHOT11-KP2
      MAT(1,2)=GHOTPM*KPL**2
      MAT(1,3)=KPL*(KPAR + GHOT13)
      MAT(2,1)=GHOTPM*KMI**2
      MAT(2,2)=K02*RSTIX-TRA*GHOT22-KP2
      MAT(2,3)=KMI*(KPAR + GHOT23)
      MAT(3,1)=KMI*(KPAR + GHOT13)
      MAT(3,2)=KPL*(KPAR + GHOT23)
      MAT(3,3)=K02*PSTIX-TRA*GHOT33
 
C     DO 2 I=1,3
C     DO 2 J=1,3
C     MAT(I,J)=MAT(I,J)/MATNOR(I)
C 2   CONTINUE
 
C     DO 3 J=1,3
C     DO 3 I=1,3
C     MAT(I,J)=MAT(I,J)/MATNOR(J)
C 3   CONTINUE
C
      LDAP=3
      LDFAC=3
      NMAT=3
C     WRITE(6,*)MAT
      do j = 1, nmat
      do i = 1, nmat
      fac(i,j) = mat(i,j)
      end do
      end do
      call zgetrf(nmat, nmat, fac, ldfac, IPVT, info)
c      CALL LFTCG(NMAT,MAT,LDAP,FAC,LDFAC,IPVT)
C     WRITE(6,*) FAC,IPVT
 
c      CALL LFDCG(NMAT,FAC,LDFAC,IPVT,DET1,DET2)
c      DETPHY=DET1*10.D0**DET2*DETNOR
      sidet = 1.
      if(ipvt(1).ne.1)sidet=-sidet
      if(ipvt(2).ne.2)sidet=-sidet
      if(ipvt(3).ne.3)sidet=-sidet
      detphy = detnor * sidet * fac(1,1)*fac(2,2)*fac(3,3)
 
      END IF
 
      RETURN
      END
