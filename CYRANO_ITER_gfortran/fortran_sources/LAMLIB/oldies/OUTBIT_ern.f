      SUBROUTINE OUTBIT_ern

      IMPLICIT NONE

C     BUILDS and writes some 1D tables:
c	antenna poloidal and toroidal spectra, wall current
c     Also writes output for QLFP
      
      include 'PARDIM.COPY'
      include 'DYNOU2.COPY'
      include 'COMGEO.COPY'
      include 'COMANT.COPY'
      include 'COMMOD.COPY'
      include 'COMREG.COPY'
      include 'COMSWE.COPY'
      include 'COMFOU.COPY'
      include 'COMPLO.COPY'
      include 'complp.copy'
      include 'COMPLA.COPY'
      include 'COMFIN.COPY'
      include 'COMFIC.COPY'
      include 'COMIN2.COPY'
      include 'COMGDR.COPY'
      include 'COMPHY.COPY'
       
 
      INTEGER 
     ;  i, ia, iae, ir
     ;, iplom1, ipol, ii
     ;, ic, im
     ;, ima, ist 
     ;, k1, k2, ki, imk, nsum

      DOUBLE PRECISION 
     ;  jthwa(maxthp+1,2)
     ;, wallos, norfac(3)
     ;, trapeze, gauss
     ;, second
 
      COMPLEX*16
     ;  c1
     ;, zdotu
     ;, tra1, tra2
     ;, tofac, toi
     ;, sum1, sum2

      EXTERNAL 
     ;  trapeze, gauss, nsum
 
     
      EXTERNAL tofac, SECOND
 
c      CABS2(TRA1) = DREAL(TRA1 * DCONJG(TRA1))
 
C**********************************************************************

c     Output for QLFP:
      include 'OPENFI.COPY' 
c     NB: if nspgdr=0, nraddr and the abscissae in raqlfp
c     must be given in NAMPLA!
      write(undirk,*)nraddr
        do ir = 1, nraddr
c        i = isrchfge(ist11, absout, 1, raqlfp(ir))
c     Assumes absout and abscis identical:
        i = irafp(ir)
        write(nofile,*)'writing rf fields for QLFP on unit ', undirk
        write(nofile,*)'for minor radius rho= ', absout(i)
c       Field normalization factor, assumes power RFPOW into 1 toroidal mode:
c       sqrt(2) comes from Dirk's convention on E+, E_.
        ii = 1
        if(flrops(1))ii = 2
 
        call dset(3, 1.d0, norfac, 1)
          if(parpow(nspec+1,ii).ne.0.d0)then
			norfac(3) = dsqrt(rfpow / parpow(nspec+1,ii))
			norfac(1) = sqrt2 * norfac(3)
			norfac(2) = norfac(1)
          end if
c       Write index of radius and radius:
        write(undirk,*)ir, absout(i)
c       Pol. mode interval, assumed the same for all elements:
        write(undirk,*)minf(2), msup(2)
c       Write +, -, // components of each mode:
          do im = 1, nmode(2)
          write(undirk,2000)(norfac(ic) * xpmp(2*ic-1,i,im), ic = 1, 3)
          end do
        end do
 
      iplom1 = nploth + 1

C     COMPUTE L2NORM OF SOLUTION:
      DO II = 1, NCOMP
      EL2NO(II) = CZERO
      IST = 1
        DO IREG = 1, NREG
cERN          IF(OUTGAU)THEN
          IST = IST + 1
            DO IEL = IFIEL(IREG), ILAEL(IREG)
            EL2NO(II) = EL2NO(II)
     ;       + GAUSS( RNORM*FL(IEL), ABSOUT(ist), ENORM(IST,II), 1, 0 )
            IST = IST + NGAUSS + 1
c            IST = IST + NGAUSS
            END DO
          IST = IST + 1
cERN           ELSE
cERN           EL2NO(II) = EL2NO(II)
cERN      ;       + TRAPEZE( ISTP(IREG), ABSOUT(IST), ENORM(IST,II), 1 ,0 )
cERN          IST = IST + ISTP(IREG)
cERN          END IF
        END DO
cERN	What normalization is this??
      EL2NO(II) = DSQRT( EL2NO(II) / (RX0M(NREG) - RX0M(0)) )
      END DO
 
C     Moduli of antenna current norms:
        do im = 1, nmoant
          tra1 = czero
          tra2 = czero
          do ia = 1, ntoant
             if(samant)then
               iae = 1
             else
               iae = ia
             end if
             toi = tofac(ia) * tancur(ia)
             tra1 = tra1 + cura(im,iae)*toi
             tra2 = tra2 + chaa(im,iae)*toi
          end do
          ancplo(im,1) = cdabs(tra1)
          ancplo(im,2) = cdabs(tra2) / eps0
      end do
 
cERN  cccccccccccc Write Antenna current poloidal modes ccccccccccccc
	open (UNIT = 701, FILE = paplda // 'Antenna_polspec.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		  write(701,*),"Antenna current and charge/Eps0 modes moduli"
            do i = 1, nmoant
		     write(701,"(i5,g16.6,g16.6)") moant(i),ancplo(i,1),ancplo(i,2)
		  end do
	close(701)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cERN  cccccccccccc Write Antenna current toroidal spectrum cccccccccc
	open (UNIT = 701, FILE = paplda // 'Antenna_torspec.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		  write(701,*),"Antenna current toroidal spectrum"
            write(701,"(i5)"), int(motoan(1))
		  do i = 1, nptosp
		     write(701,"(i5,g16.6)") int(abstom(i)),antosp(i)
		  end do
	close(701)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c     Compute wall surface current and wall losses as perturbation:
c     -------------------------------------------------------------
      iel = nele
      ielm = iel
      if(nmode(iel-1) .gt. nmode(iel))ielm = iel - 1
      ima = minf(ielm) + 1 - minf(imax)
c     hphi for 'm23', heta for 'hec', yield respectively jtheta, j//:
        do ipol = 1, nploth
        c1 = zdotu(nmode(ielm), hrtp(3,nabsci,1), 3*nabplo
     ;, tabexp(ima,ipol), 1)
        jthwa(ipol,1) = dreal(c1)
        jthwa(ipol,2) = dimag(c1)
        end do
      jthwa(iplom1,1) = jthwa(1,1)
      jthwa(iplom1,2) = jthwa(1,2)

cERN  cccccccccccc Write wall poloidal current  ccccccccccccccccccccc
	open (UNIT = 701, FILE = paplda // 'Wall_current.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
		  write(701,*),"Wall poloidal current"
            do i = 1, iplom1
		     write(701,"(g16.6,g16.6,g16.6)"),polplo(i),jthwa(i,1),jthwa(i,2)
		  end do
	close(701)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
c     Wall losses:
      Y = RX0(NREG)
      IEL = NELE
      isubr = nsum(nreg, ns)
      i = IST11
      INTAB = INTAPLOT(i)
      SUM1 = CZERO
        if(circ)then
        CALL FOUCO
          DO M = MINF(IEL), MSUP(IEL)
          IM = M + 1 - MINF(IEL)
          SUM2 = 2.D0 * KPRN
     ;                * DIMAG( DCONJG( XRTP(5,I,IM) ) * XRTP(1,I,IM) )
            IF(CYL)THEN
            K1 = 0
            K2 = 0
            ELSE
            K1 = max0(- NCROT, minf(iel)-m)
            K2 = min0(NCROT, msup(iel)-m)
            END IF
            DO K = K1, K2
            KI = IABS(K)
            IMK = IM + K
            SUM2 = SUM2 + AFOU(KI) * (
     ;         CONJG( XRTP(5,I,IM) ) * XRTP(5,I,IMK)
     ;       + CONJG( XRTP(3,I,IM) ) * XRTP(3,I,IMK)
     ;       + M*(M+K) / Y**2 * CONJG( XRTP(1,I,IM) ) * XRTP(1,I,IMK)
     ;       + 2.D0*(M+K) / Y * DIMAG( DCONJG( XRTP(3,I,IM)) * XRTP(1,I,IMK) )
     ;       )
     ;       + BFOU(KI) * KPRN**2
     ;                  * CONJG( XRTP(1,I,IM) ) * XRTP(1,I,IMK)
            END DO
          SUM1 = SUM1 + SUM2
          END DO

        else
c        CALL FOUCOG(2)
c...to do!
        end if

        if(sum1 .ne. czero)then
        write(nofile,*)'Loss term: sum1=', sum1
        wallos = dreal(sum1) * 0.5d0 * walres / skiwal
     ;           * y / rnorm * (twopi / omegag * mu0)**2
          if(cyl)then
          write(nofile,*) 'Losses in wall=', WALLOS, ' Watt/m'
          else
          wallos = wallos * ra
          write(nofile,*) 'Losses in wall=', WALLOS, ' Watt'
          end if
        if(.not.glovac .and. .not.losles)
     ;  write(nofile,*)'Losses in wall=', 100*WALLOS/TINPOW, '% of input power'
        end if
      
C      WRITE(NOFILE,*)'Plot legend table:'
C      WRITE(NOFILE,*)'------------------'
C      WRITE(NOFILE,1002) (I,PLOLEG(I),I=1,IGLPLO)
C      WRITE(NOFILE,*)'EXIT OUTPUT ; time=',SECOND()-TIMIN
 
      return
 
 1000 format(1h ,'Plot #',i4,'  ',a50)
 1002 format(1h ,i4,'  ',a50)
 1003 format(10g23.15)
 1488 format((' ',4(1x,1p,d11.3)))    
 1588 format(' ',3(1x,i4))          
 2000 format(1h ,21(2x,g13.5))

      end