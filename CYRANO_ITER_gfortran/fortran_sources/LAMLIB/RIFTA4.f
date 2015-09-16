      subroutine riftab(nabs, zarg, zkll, rift, tr, tr2, itr, onlyab)

      implicit none

      logical onlyab

      integer nabs, itr(nabs)

      real*8 zarg(nabs), zkll(nabs), tr(nabs), tr2(nabs)

      complex*16 rift(nabs)
c
c plasma dispersion function for real values of k// and omega
c version of the zif4 function further optimized for cray
c with tabulation on uniform mesh and quadratic lagrange interpolation.
c cyrano module - p.lamalle 1994
c 28 jan 2000: option onlyab added
c
c input:
c                 onlyab:        if .true., only compute imaginary parts.
c                 nabs           number of function evaluations
c   raymond's convention on vll:
c   arrays        zarg(nabs)     arguments (omega-n*omc)/sqrt(2)*kll*vll
c   (in cyrano, vth = sqrt(2)*vll )
c                 zkll(nabs)     associated k//
c   work arrays   tr(nabs),tr2(nabs) : real
c                 itr(nabs) : integer
c
c output:         rift(nabs)     associated values of plasma disp.func.
 
C     Table dimensions: these values give max.absolute error 4.D-9
C     with respect to ZIF4 function.
      integer nptab, nptaba

      parameter (nptab=2001, nptaba=251)

      integer 
     ;  ini, i, ipoi1(2), ipoi2(2), npt(2), nabs1, nabs2, nal, ipoi, ireg
     ;, itab, itab1, itab2 

      real*8 x, hs, ha, h(2), hi(2), brkpts, xabs, fi, p, pp1, pm1, a1, a2, a3

      complex*16 
     ;  table(0:nptab), tabla(0:nptaba), zif4, tab(0:nptab,2)

      equivalence (tab(0,1),table(0)), (tab(0,2),tabla(0))

      external zif4
 
      data ini/0/
      save ini
      data brkpts/4.d0/
 
        if(ini.eq.0) then
        ini=1
c       table of nptab points on [0,brkpts]:
        hs = brkpts / dfloat(nptab-1)
          do 301 i = 1, nptab-1
          table(i) = zif4(dfloat(i-1)*hs, 1.d0, .false.)
301       continue
        table(0) = table(1)
        table(nptab) = zif4(brkpts, 1.d0, .false.)
 
c       Table of nptaba points on [0,1/brkpts]for variable 1/x:
        ha = 1.d0/(brkpts*dfloat(nptaba-1))
        tabla(1) = 0.
          do 303 i = 2, nptaba-1
          tabla(i) = zif4(1.d0/(dfloat(i-1)*ha), 1.d0, .false.)
303       continue
        tabla(0) = tabla(1)
        tabla(nptaba) = table(nptab)
 
        h(1)=hs
        h(2)=ha
        hi(1) = 1./hs
        hi(2) = 1./ha
        npt(1)=nptab
        npt(2)=nptaba
        ipoi1(1)=1
        end if
 
      call zset(nabs, (0.d0,0.d0), rift, 1)
 
      nabs1 = 0
        do 2 ipoi = 1, nabs
        xabs = dabs(zarg(ipoi))
          if(xabs.le.brkpts)then
          nabs1 = nabs1 + 1
          itr(nabs1) = ipoi
          tr(nabs1) = zarg(ipoi)
          end if
  2     continue
      nabs2 = 0
        if(nabs1.lt.nabs)then
          do 7 ipoi = 1, nabs
          xabs = dabs(zarg(ipoi))
            if(xabs.gt.brkpts)then
            nabs2 = nabs2 + 1
            nal = nabs1 + nabs2
            itr(nal) = ipoi
            tr(nal) = 1.d0 / zarg(ipoi)
            end if
  7       continue
        end if
 
        do 3 ipoi = 1, nabs
        tr2(ipoi) = zkll(itr(ipoi))
  3     continue
c
      ipoi1(2)=nabs1+1
      ipoi2(1)=nabs1
      ipoi2(2)=nabs
      
      if(onlyab)then
c     --------------
      do ireg = 1, 2
      do ipoi = ipoi1(ireg), ipoi2(ireg)
        if(tr2(ipoi).ne.0.d0) then
        x = tr(ipoi)
        xabs = dabs(x)
c       Interpolate in table: 2 point Lagrange
c       Equidistant case:
        fi = xabs * hi(ireg)
        itab = min0( 2 + int(fi), npt(ireg)-1 )
        itab1 = itab - 1
        itab2 = itab + 1
        p = fi - dfloat(itab1)
        pp1 = p + 1.
        pm1 = p - 1.
        a1 =  p*pm1*0.5
        a2 = -pm1*pp1
        a3 =  p*pp1*0.5
        rift(itr(ipoi)) = dcmplx(0.d0, 
     ;   a1*dimag(tab(itab1,ireg)) 
     ; + a2*dimag(tab(itab,ireg))
     ; + a3*dimag(tab(itab2,ireg))
     ; )
        end if
      end do
      end do
      else
c     ----
      do 1 ireg = 1, 2
      do 1 ipoi = ipoi1(ireg), ipoi2(ireg)
c     ====================================
      if( tr2(ipoi).ne.0.d0) then
      x = tr(ipoi)
      xabs = dabs(x)
c     Interpolate in table: 2 point Lagrange
c     Equidistant case:
      fi = xabs * hi(ireg)
      itab = min0( 2 + int(fi), npt(ireg)-1 )
      itab1 = itab - 1
      itab2 = itab + 1
      p = fi - dfloat(itab1)
      pp1 = p + 1.
      pm1 = p - 1.
      a1 =  p*pm1*0.5
      a2 = -pm1*pp1
      a3 =  p*pp1*0.5
      rift(itr(ipoi)) =
     ;    a1*tab(itab1,ireg) + a2*tab(itab,ireg) + a3*tab(itab2,ireg)
      end if
  1   continue
c     ========
      end if
c     ------
 
      do 5 ipoi = 1, nabs
      if(zkll(ipoi).lt.0.d0)rift(ipoi) = conjg(rift(ipoi))
  5   continue

      do 6 ipoi = 1, nabs
      if(zarg(ipoi).lt.0.d0)rift(ipoi) =
     ;dcmplx(-dreal(rift(ipoi)),dimag(rift(ipoi)))
c    ;- dconjg(rift(ipoi))
  6   continue
 
      return
      end

