      SUBROUTINE OUTBIT

      IMPLICIT NONE

C     BUILDS and writes OUTPUT TABLES, to be plotted elsewhere.
C     New, shortened version (1998) of OUTPUT, part 5.
c     Also writes output for QLFP
      
      
      include 'pardim.copy'
      include 'dynou2.copy'
      include 'compow.copy'
      include 'comgeo.copy'
      include 'comant.copy'
      include 'commod.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comswe.copy'
      include 'com3di.copy'
      include 'comfou.copy'
      include 'comro2.copy'
      include 'complo.copy'
      include 'complp.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comfin.copy'
      include 'comfic.copy'
      include 'comequ.copy'
      include 'comin2.copy'
      include 'comgdr.copy'
      include 'comphy.copy'

      character*9  charfi(2), charco(6), charpa(3)
 
      integer 
     ;  i, j, ia, iae, ir
     ;, iplom1, ipol, ii
     ;, ic, im
     ;, ima, ist 
     ;, k1, k2, ki, imk
     ;, nsum, i1st, ceilq
     ;, idamin, idmin, idmax, isrchfge

      double precision 
     ;  jthwa(maxthp+1,2), jphwa(maxthp+1,2), absmon(maxpom)
     ;, wallos, norfac(3)
     ;, cabs2, rveno
     ;, trapeze, gauss
     ;, second
 
      complex*16
     ;  c1
     ;, zdotc, zdotu
     ;, tra1, tra2
     ;, tofac, toi, cdeid
     ;, sum1, sum2

      external 
     ;  nsum, i1st, ceilq, trapeze, gauss, rveno
     ;, zdotc, zdotu, zset, dset, dsum, mucrvz
     ;, cdeid
     ;, idamin, idmin, idmax, isrchfge
 
c      character*7 jobnam
c      character*10 texnum, tex1
c      character*33 ch30
      character*48 text(14), title, yla(maxspe)
 
     
      external tofac, second
c
      data charfi/'Electric ','Magnetic '/
     ;    ,charco/'radial ','poloidal ','toroidal ',
     ;            'L.H.(+) ','R.H.(-) ','parallel '/
     ;    ,charpa/'real','imag','modulus'/
 
      cabs2(tra1) = dreal(tra1 * dconjg(tra1))
 
c**********************************************************************

c     Output for QLFP:
      include 'openfi.copy' 
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
          IF(PARPOW(NSPEC+1,ii).NE.0.D0)THEN
          NORFAC(3) = DSQRT(RFPOW / PARPOW(NSPEC+1,ii))
          NORFAC(1) = SQRT2 * NORFAC(3)
          NORFAC(2) = NORFAC(1)
          END IF
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
cERN	     IF(OUTGAU)THEN  
        IST = IST + 1
          DO IEL = IFIEL(IREG), ILAEL(IREG)
          EL2NO(II) = EL2NO(II) + GAUSS( RNORM*FL(IEL), ABSOUT(ist), ENORM(IST,II), 1, 0 )
          IST = IST + NGAUSS + 1
c          IST = IST + NGAUSS
          END DO
        IST = IST + 1
cERN          ELSE
cERN            EL2NO(II) = EL2NO(II)
cERN     ;      + TRAPEZE( ISTP(IREG), ABSOUT(IST), ENORM(IST,II), 1 ,0 )
cERN            IST = IST + ISTP(IREG)
cERN          END IF ! OUTGAU = T
        END DO
      EL2NO(II) = DSQRT( EL2NO(II) / (RX0M(NREG) - RX0M(0)) )
      END DO
 
C     Moduli of antenna current norms:
      WRITE(75,*)'Antenna current and charge/Eps0 modes moduli:'
      WRITE(75,*)nmoant, 3
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
          tra1 = tra1 + cura(im,iae) * toi
          tra2 = tra2 + chaa(im,iae) * toi
          end do
        ancplo(im,1) = cdabs(tra1)
        ancplo(im,2) = cdabs(tra2) / eps0
        write(75,*)float(moant(im)), (ancplo(im,j), j = 1, 2)
        end do
 
cERN  cccccccccccc Write Antenna current poloidal modes ccccccccccccc
	open (UNIT = 701, FILE = paplda // 'Antenna_polspec.dat', STATUS = "REPLACE", ACTION = "WRITE")
      write(701,*),"Antenna current and charge/Eps0 modes moduli"
        do i = 1, nmoant
        write(701,"(i5,g16.6,g16.6)") moant(i),ancplo(i,1),ancplo(i,2)
        end do
	close(701)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Plot poloidal mode norms, toroidal spectrum, charge modes:
c        IF( .NOT. OUTEZ(3*nficom+3) )THEN
c          DO I = MINF(IMAX), MSUP(IMAX)
c          IR = I + 1 - MINF(IMAX)
c          ABSMON(IR) = I
c          END DO
c        TEXT(9)='                                    VOLT/M'
c                    IF(CYL)THEN
c      CALL GRFTOC(KPHI,TEX1,II)
c      TEX3 = ' KPHI='  // TEX1(1:II)//' M**-1'
c                    ELSE
c      CALL GRFTOC(N,TEX1,II)
c      TEX3 = ' N=' // TEX1(1:II)
c                    END IF
c      title = 'Poloidal modes L2 norms                         '
c      IGLPLO = IGLPLO + 1
c      WRITE(NOFILE,1000) IGLPLO, title
c      PLOLEG(IGLPLO) = title
c      call wrp1d(nofile, IGLPLO, nmode(imax), absmon, monorm
c     ;, 1, 1, 1
c     ;, title
c     ;, 'Index     ', yla
c     ;, '          ', '                                    Volt/m'
c     ;, 0, TEXT,.false.)
c       END IF
 
c     Antenna current modes:
        do ir = 1, nmoant
        absmon(ir) = moant(ir)
        end do
      title = 'Antenna current toroidal spectrum               '
      iglplo = iglplo + 1
      write(nofile,1000)iglplo, title
      ploleg(iglplo) = title
      call wrp1d(75, iglplo, nptosp, abstom, antosp, nptosp, 1, 1, title
     ;, '          ', yla
     ;, '          ', yla
     ;, 0, text,.false.)

cERN  cccccccccccc Write Antenna current toroidal spectrum cccccccccc
	open(UNIT = 701, FILE = paplda // 'Antenna_torspec.dat', STATUS = "REPLACE", ACTION = "WRITE")
      write(701,*),"Antenna current toroidal spectrum"
      write(701,"(i5)"), int(motoan(1))
        do i = 1, nptosp
        write(701,"(i5,g16.6)") int(abstom(i)),antosp(i)
	  end do
	close(701)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      title = 'Antenna current and (charge/eps0) mode ampl.    '
c      IGLPLO = IGLPLO + 1
c      WRITE(NOFILE,1000) IGLPLO, title
c      PLOLEG(IGLPLO) = title
c      call wrp1d(75, IGLPLO, nmoant, absmon, ancplo
c     ;, maxpom, 1, 2
c     ;, title
c     ;, '          ', yla
c     ;, '          ', yla
c     ;, 0, TEXT,.false.)
 
c     Compute wall surface current, and wall losses as perturbation:
c     -------------------------------------------------------------
      iel = nele
      ielm = iel
      if(nmode(iel-1) .gt. nmode(iel))ielm = iel - 1
      ima = minf(ielm) + 1 - minf(imax)
cPL24/01/2005: must make this consistent!
c     In general geometry: Hphi for 'M23', Heta for 'HEC', yield respectively Jtheta, J//:
c     In circular geometry: Hphi for both element types, yields Jtheta:
        if(.not.cokpco)then
          do ipol = 1, nploth
          c1 = zdotu(nmode(ielm), hrtp(3,nabsci,1), 3*nabplo, tabexp(ima,ipol), 1)  ! cern ist11->nabsci
          jthwa(ipol,1) = dreal(c1)
          jthwa(ipol,2) = dimag(c1)
          end do
          do ipol = 1, nploth
          c1 = zdotu(nmode(ielm), hrtp(2,nabsci,1), 3*nabplo, tabexp(ima,ipol), 1)  ! cern ist11->nabsci
          jphwa(ipol,1) = - dreal(c1)
          jphwa(ipol,2) = - dimag(c1)
          end do
	  else ! cokpco case, circular (NB: R, Z matrices must be output for thetabar mesh in noncircular)
cPL24/01/05 only correct when OUTNPFT = .true; otherwise must interpolate to grid of nploth points, different from npfft
          do ipol = 1, nploth
c          c1 = zdotu(nmode(ielm), hrtp(3,nabsci,1), 3*nabplo, tabexp(ima,ipol), 1) * cdeid(n*ckt(nabsci,ipol,4)) 
cPL24/01/05 Using new array einphb:
          c1 = zdotu(nmode(ielm), hrtp(3,nabsci,1), 3*nabplo, tabexp(ima,ipol), 1) * einphb(nabsci,ipol)
          jthwa(ipol,1) = dreal(c1)
          jthwa(ipol,2) = dimag(c1)
          end do
          do ipol = 1, nploth
          c1 = zdotu(nmode(ielm), hrtp(2,nabsci,1), 3*nabplo, tabexp(ima,ipol), 1) * einphb(nabsci,ipol)
          jphwa(ipol,1) = - dreal(c1)
          jphwa(ipol,2) = - dimag(c1)
          end do
	  end if
      jthwa(iplom1,1) = jthwa(1,1)
      jthwa(iplom1,2) = jthwa(1,2)
      jphwa(iplom1,1) = jphwa(1,1)
      jphwa(iplom1,2) = jphwa(1,2)
c      TEXT(4)= ' '
c      TEXT(11)=' '
c      TEXT(12)=' '
c      CH30='Wall poloidal current'
c      TEXT(9)='                                  AMPERE/M'
c      TEXT(2)='          Poloidal angle (radian)'
c      CALL GRFTOC(FREGAG*1.D-6,TEXNUM,II)
c      TEXT(13)= CH30 // '    ' // TEXNUM(1:II) // 'MHZ'
c      TEXT(14)= TEXMOD
c      TEXT(1)='#Real part'
c      TEXT(3)='#Imag.part'
 
      title = 'Wall poloidal current (Re, Im)                  '
      iglplo = iglplo + 1
      write(nofile,1000)iglplo, title
      ploleg(iglplo) = title
      call wrp1d(75, iglplo, iplom1, polplo, jthwa, maxthp+1, 1, 2
     ;, title, 'Radian    ', 'Poloidal angle                                  ', '          ', yla, 0, TEXT,.false.)

cERN  cccccccccccc Write wall poloidal current  ccccccccccccccccccccc
	open (UNIT = 701, FILE = paplda // 'Wall_current.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
      write(701,*),"Wall poloidal current"
        do i = 1, iplom1
        write(701,"(g16.6,g16.6,g16.6)"),polplo(i),jthwa(i,1),jthwa(i,2)
	  end do
	close(701)
cPL  cccccccccccc Write wall toroidal current  cccccccccccccccccccccc
	open (UNIT = 701, FILE = paplda // 'Wall_toroidal_current.dat', 
     ;      STATUS = "REPLACE", ACTION = "WRITE")
      write(701,*),"Wall toroidal current"
        do i = 1, iplom1
        write(701,"(g16.6,g16.6,g16.6)"),polplo(i),jphwa(i,1),jphwa(i,2)
	  end do
	close(701)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Wall losses:
      y = rx0(nreg)
      iel = nele
      isubr = nsum(nreg, ns)
      i = ist11
      intab = intaplot(i)
      sum1 = czero
cPL25/1/2005 use wall current to compute wall dissipation:
        if(.not.cyl)then
          do i = 1, nploth
	    sum1 = sum1 + (jthwa(i,1)**2  + jthwa(i,2)**2 + jphwa(i,1)**2 + jphwa(i,2)**2) * eqt(intab,i,1) * eqt(intab,i,7)
          end do
c       Integral of (surface current squared * R * Ntheta) over outer wall area:
	  sum1 = sum1 * rnorm * y * twopi**2 / dfloat(nploth)
	  else
          do i = 1, nploth
	    sum1 = sum1 + (jthwa(i,1)**2  + jthwa(i,2)**2 + jphwa(i,1)**2 + jphwa(i,2)**2) * eqt(intab,i,7)
          end do
c       Integral of (surface current squared * R * Ntheta) over outer wall perimeter:
	  sum1 = sum1 * rnorm * y * twopi / dfloat(nploth)
	  end if
c     Wall dissipation:
        if(sum1 .ne. czero)then
        write(nofile,*)'Wall loss term: sum1=', sum1
        write(nofile,*)'Wall resistivity: ', walres, ' Ohm.m'
        write(nofile,*)'Wall skin depth: ', skiwal, 'm'
        wallos = sum1 * 0.5d0 * walres / skiwal
          if(cyl)then
          write(nofile,*) 'Losses in wall=', WALLOS, ' Watt/m'
          else
          write(nofile,*) 'Losses in wall=', WALLOS, ' Watt'
          end if
        if(.not.glovac .and. .not.losles)write(nofile,*)'Losses in wall=', 100*WALLOS/TINPOW, '% of input power'
        end if
      
c     Old stuff:
c        if(circ)then
c        call fouco
c          do m = minf(iel), msup(iel)
c          im = m + 1 - minf(iel)
c          sum2 = 2.d0 * kprn * dimag(dconjg(xrtp(5,i,im)) * xrtp(1,i,im))
c            if(cyl)then
c            k1 = 0
c            k2 = 0
c            else
c            k1 = max0(- ncrot, minf(iel)-m)
c            k2 = min0(ncrot, msup(iel)-m)
c            end if
c            do k = k1, k2
c            ki = iabs(k)
c            imk = im + k
c            sum2 = sum2 + afou(ki) * (
c     ;         conjg( xrtp(5,i,im) ) * xrtp(5,i,imk)
c     ;       + conjg( xrtp(3,i,im) ) * xrtp(3,i,imk)
c     ;       + m*(m+k) / y**2 * conjg( xrtp(1,i,im) ) * xrtp(1,i,imk)
c     ;       + 2.d0*(m+k) / y * dimag( dconjg( xrtp(3,i,im)) * xrtp(1,i,imk) )
c     ;       )
c     ;       + bfou(ki) * kprn**2
c     ;                  * conjg( xrtp(1,i,im) ) * xrtp(1,i,imk)
c            end do
c          sum1 = sum1 + sum2
c          end do
c
c        else
c        call foucog(2)
cc...to do!
c        end if
c        if(sum1 .ne. czero)then
c        write(nofile,*)'Loss term: sum1=', sum1
c        wallos = dreal(sum1) * 0.5d0 * walres / skiwal * y / rnorm * (twopi / omegag * mu0)**2
c          if(cyl)then
c          write(nofile,*) 'Losses in wall=', WALLOS, ' Watt/m'
c          else
c          wallos = wallos * ra
c          write(nofile,*) 'Losses in wall=', WALLOS, ' Watt'
c          end if
c        if(.not.glovac .and. .not.losles)write(nofile,*)'Losses in wall=', 100*WALLOS/TINPOW, '% of input power'
c        end if
      
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
