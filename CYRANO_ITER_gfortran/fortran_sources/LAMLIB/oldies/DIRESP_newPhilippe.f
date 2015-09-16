      SUBROUTINE DIRESP(isp)

      IMPLICIT NONE
      integer isp

C     Dielectric response of a species, using equilibrium of QLFP code
C     interpolated at local radial index INTAB of Cyrano by routine INTF0
C     isp is index of species among the ones requiring general diel. response
C     Absolute species index is ISPE

      include 'PARDIM.COPY'
      include 'COMREG.COPY'
      include 'COMSUB.COPY'
      include 'COMGEO.COPY'
      include 'COMEQU.COPY'
      include 'COMANT.COPY'
      include 'COMFIN.COPY'
      include 'COMFIC.COPY'
      include 'COMIN2.COPY'
      include 'COMROT.COPY'
      include 'COMMAG.COPY'
      include 'COMMA2.COPY'
      include 'COMPLA.COPY'
      include 'COMFOU.COPY'
      include 'COMMOD.COPY'
      include 'COMMMK.COPY'
      include 'COMPRI.COPY'
      include 'COMSWE.COPY'
      include 'COMGDR.COPY'
      include 'COMUSR2.F'
      include 'COMPHY.COPY'
 
      logical pass, cxr, toini

      save toini
      
      integer i, j, im, ix, jv, ip, iex, jev, nset0, isig, igv, igx
     ;, igv2, igx2, jev1, jev2, iex1, iex2
     ;, ispe, ii, irv, jvn, irx, ixn, icount, ic, ipov, ipox, ih
     ;, ji1, jm1, i1, im1, ib, jb
     ;, iplb, kstop1, kstop2, L, La
     ;, isrchfgt

C     act, react: im and re of iL      
c     See if they could be equivalenced to, say, VMAT in comusr2,
c     or other memory saving?
      double precision
     ;  act(2*maxpom-1,0:maxcou,-1:1)
     ;, react(2*maxpom-1,0:maxcou,-1:1)
     ;, ac(2*maxpom-1,0:maxcou,-1:1,3)
     ;, reac(2*maxpom-1,0:maxcou,-1:1,3)
      equivalence (act,ac), (react,reac)

      double precision 
     ;  v, v2, vi, xn, xna, x
      double precision 
     ;  u0, u0i, u1, u2, ru1, ru2, iu1, iu2, xtir1, khi
      double precision 
     ;  khi1, lamb, del2, del, te0, t1, t2, t3, t4, t5
     ;, cmr1, cmr2, smr1, smr12, smr2, smr22, cmr0
     ;, smr0, smr02, csc22
     ;, clmr0(0:maxcou), clmr1(0:maxcou), clmr2(0:maxcou)
     ;, er0, r00, r1, r2, r3, s0, s1, s1i, s2, s2i, omxnai, xmax
     ;, hp, hp2, hp3, hp4, hp5, hp6, omca, xdb, xdbi, sxdb, sxdbi
     ;, petita(5*maxnex), petitb(5*maxnex)
     ;, nze, N1, L0, L1, jze, J1, alp2, alt2
     ;, NL(5*maxnex,-1:maxcou+3), JL(-1:maxcou+3,5), LL(-1:maxcou+3,5)
     ;, tj1(maxcou+2), tj2(maxcou+2), tj3(maxcou+2), tj4(maxcou+2)
     ;, cmmv(5*maxnex)
     ;, kel, kelv(5*maxnex), s2k
     ;, maomi
     ;, ek(5*maxnex), ee(5*maxnex), cmm, w0, w1, w2
     ;, oo, o1, o2, o3, o4, o5, fis
     ;, bn(4,4)
     ;, ellfc, ellec, ellpic

      complex*16
     ;  JLc(-1:maxcou+1), LLc(-1:maxcou+1), cu1, cdellpic
     ;, ccmr1, csmr12
     
      equivalence (JLc, JL), (LLc, LL)

      external ellfc, ellec, ellpic, cdellpic

      data toini/.true./

c       Coefficients for NL recurrences:
        if(toini)then
          do j = 1, maxcou+2
          tj1(j) = dfloat(2*j) / (dfloat(j)-0.5d0)
          tj2(j) = (dfloat(j)+0.5d0) / (dfloat(j)-0.5d0)
          tj3(j) = dfloat(2*j) / (dfloat(j)+0.5d0)
          tj4(j) = (dfloat(j)-0.5d0) / (dfloat(j)+0.5d0)
          end do
        toini = .false.
        end if
      
c     indices: mbar, m1-m2, p

      ispe = ispgdr(isp)
      
      nset0 = (2*maxpom-1)*(maxcou+1)*3
      call dset(3*nset0, 0.d0, ac, 1)      
      call dset(3*nset0, 0.d0, reac, 1)
      call dset(nset0, 0.d0, tact, 1)      
      call dset(nset0, 0.d0, treact, 1)
      clmr0(0) = 1.d0
      clmr1(0) = 1.d0
      clmr2(0) = 1.d0

c     Normalized x at separatrix (co-passing):
      xnsep = - 2.d0 * delb(intab) / bmax(intab) 
c     Indices of left and right limits of x intervals:
c     Resp. co-passing, 'co-trapped', 'counter-trapped', counter-passing
        if(gigdr)then
c       For Gauss points loop
        igl(1) = 1
        igl(2) = isrchfgt(lxg, xngaug, 1, xnsep)
        igl(3) = isrchfgt(lxg, xngaug, 1, 0.d0)
        igl(4) = isrchfgt(lxg, xngaug, 1, - xnsep)
        idl(4) = lxg
        idl(1) = igl(2) - 1
        idl(2) = igl(3) - 1
        idl(3) = igl(4) - 1
        igv2 = ngauv
        igx2 = ngaux
        else
c       For element nodes loop
        call icopy(nxgreg, ifiax, 1, igl, 1)
        call icopy(nxgreg, ilaax, 1, idl, 1)
        igv2 = 1
        igx2 = 1
        end if

      omca = qom(ispe) * B0
      xtir1 = omca / omegag
      oo = 1.d0 / dabs(xtir1)
      te0 =  PI / (qom(ispe)*hachi(intab)*delb(intab))
      s0 = bbar(intab) / delb(intab)
      xmax = B0 / bmin(intab)
      maomi = bmax(intab) / bmin(intab)
      s1 = B0 / delb(intab)
      s1i = 1.d0 / s1
      s2 = bmin(intab) / delb(intab)
      s2i = delb(intab) / bmin(intab)
c     Normaliz. of radial response density:
      r00 = PI * (EEL*ZCH(ispe))**2
cERN ; * 4.d0 * eqta1d(intab,6) * abscis(intab)
     ; * 4.d0 * dPsidr_n(intab) * abscis(intab)
     ; / (2.d0*MH*AMASS(ispe) * B0 * hachi(intab) * omegag)

      write(6,*)'diresp: before generate NL passing'      
c     Generate the NORMALIZED NL on x mesh:
c     Passing regions:
      do ii = 1, 4, 3
        do i = igl(ii), idl(ii)
        xn = xngaug(i)
        xna = dabs(xn)
        petita(i) = dsqrt(xna)
        petitb(i) = dsqrt(maomi*(xna+xnsep))
        kel = dsqrt(1.d0-(petitb(i)/petita(i))**2)
        kelv(i) = kel
        call ellfec(kel, ek(i), ee(i))
c       Normalized NL:
        nze = 2.d0 / petita(i) * ek(i)
        if(kel.ne.0.d0)then
        cmm = - 0.5d0 * (kel + 1.d0 / kel)
        cmmv(i) = cmm
        s2k = dsqrt(2.d0 / kel)
        N1 = 2.d0 / petita(i) * (ek(i) +2.d0/(kel*kel)*(ee(i) - ek(i)))
        else
        N1 = 0.d0
c NB:   actual value is - infinity:
        cmmv(i) = 0.d0
        end if
c       backward recurrence is stable for passing: check starting index for
c       overall accuracy!
        NL(i,maxcou+3) = 0.d0
        NL(i,maxcou+2) = 1.d0
          do j = maxcou+2, 1, -1
          NL(i,j-1) = tj1(j)*cmm * NL(i,j) - tj2(j) * NL(i,j+1)
          end do
          do j = 1, maxcou+3
          NL(i,j) = NL(i,j) * (nze / NL(i,0))
          end do
        NL(i,0) = nze
        NL(i,-1) = NL(i,1)
        end do
      end do

      write(6,*)'diresp: before generate NL trapped'      
c     Trapped regions:
      do ii = 2, 3
        do i = igl(ii), idl(ii)
        xn = xngaug(i)
        xna = dabs(xn)
        petita(i) = dsqrt(xna)
c       Convention, actual b is complex for trapped:
        petitb(i) = 0.d0
        omxnai = 1.d0 / (1.d0 - xna)
        cmm = (1.d0 - xna * s0) * omxnai
        cmmv(i) = cmm
        kel = dsqrt(0.5d0 * (1.d0 - cmm))
        kelv(i) = kel
        call ellfec(kel, ek(i), ee(i))
        x = (1.d0 - dabs(xn)) * xmax
        sxdbi = dsqrt(s1 / x)
c       Normalized NL:
        NL(i,0) = 2.d0 * kel * ek(i) / petita(i)
        NL(i,1) = 2.d0 * kel * (2.d0*ek(i) - ee(i)) / petita(i)
        NL(i,-1) = NL(i,1)
          do j = 1, maxcou
          NL(i,j+1) = tj3(j)*cmm*NL(i,j) - tj4(j)*NL(i,j-1)
          end do
        end do
      end do
c      write(6,*)'diresp: before NL recursion'      
c     All: not ok, stable for trapped but unstable for passing!
c      do j = 1, maxcou
c      w2 = 1.d0 / (dfloat(j)+0.5d0)
c      w0 = dfloat(2*j) * w2
c      w1 = (dfloat(j)-0.5d0) * w2
c        do i = 1, lxg
c        NL(i,j+1) = w0 * cmmv(i) * NL(i,j) - w1 * NL(i,j-1)
c        end do
c      end do
      write(66,*)'NL'  
      write(66,*)lxg, maxcou+3
        do i = 1, lxg    
        write(66,*)xngaug(i), (NL(i,j),j=0,maxcou+1)     
        end do

      write(6,*)'DIRESP: before 2D loop'  
c     V region loop:
      do irv = 1, nvgreg
        if(gigdr)then
C       jev is an element index
        jev1 = ifielv(irv)
        jev2 = ilaelv(irv)
        else
c       jev is a node index; jv = jev
        jev1 = ifiav(irv)
        jev2 = ilaav(irv)
        end if

c     X region loop:
      do irx = 1, nxgreg

c     Element length dependent matrices (here for piecewise uniform meshes):
      call normat(vlelr(irv), xlelr(irx))
      call intnod

c       This assumes 4 x regions!:
        if(irx.le.2)then
        isig = 1
        else 
        isig = -1
        end if
      fis = dfloat(isig)
      pass = irx.eq.1 .or. irx.eq.4
        if(gigdr)then
c       iex is an element index
        iex1 = ifielx(irx)
        iex2 = ilaelx(irx)
        else
c       iex is a node index; ix = iex
        iex1 = ifiax(irx)
        iex2 = ilaax(irx)
        end if
      
c     V element loop:
c     when gigdr=.F., loop over nodes instead
      do jev = jev1, jev2
c     ===================
c     right node index (gigdr=.t.):
      jvn = irv + jev

c     X element loop:
c     when gigdr=.F., loop over nodes instead
      do iex = iex1, iex2
c     -------------------
c     right node index (gigdr=.t.):
      ixn = irx + iex
c     f0 and first derivatives at Gauss points:
      icount = 0
      if(gigdr)call fagp(jev, iex, jvn, ixn, isp)
      
      do igv = 1, igv2
c     ================
      jv = (jev - 1) * igv2 + igv
      v = vgaug(jv)
      v2 = v * v
      vi = 1.d0 / v
      er0 = r00 * v2
      r1 = er0 * v2 * PI
      
        do igx = 1, igx2
c       ~~~~~~~~~~~~~~~~
        ix = (iex - 1) * igx2 + igx
        xn = xngaug(ix)
        x = (1.d0 - dabs(xn)) * xmax
        khi1 = x / xtir1
C       x*delta/B0:
        xdb = x / s1
        xdbi = 1.d0 / xdb
        sxdb = dsqrt(xdb)
        sxdbi = 1.d0 / sxdb
        
c     Find resonances:
c loop over ip=±1 to add below!?
      ip = 1
      khi = ip * khi1

c     Fragment stored in GIGDR.f: Gaussian integration in x and v, to update if ever needed!
c     INCLUDE 'GIGDR.f'

      if(.not.gigdr)then
c     Semi analytical integration over elements
c     =========================================
c     Groupe selon coef. de f0, df0/dv, df0/dx (cf dernier indice ac et reac)

      r3 = er0 * PI * s1 * 0.25d0 * omegag / dabs(omca)

c     Local terms are indep. of k//:
      o1 = 1.d0 - x * bbar(intab) / B0
c     Using normalized NL:
      do i = 0, maxcou
      reac(1,i,0,1) = er0 * (0.5d0 * (1.d0 - 3.d0 * o1) * NL(ix,i)
     ;                       - 0.75d0 * xdb * (NL(ix,i+1) + NL(ix,i-1)))
      reac(1,i,0,3) = - 2.d0 * x * er0 * (
     ;0.5d0 * o1 * NL(ix,i) 
     ;+ 0.25d0 * xdb * (NL(ix,i+1) + NL(ix,i-1)))

      o2 = - 0.25d0 * xdb * (s0 * NL(ix,i) - 0.5d0 * (NL(ix,i+1) + NL(ix,i-1)))
      reac(1,i, 1,1) = er0 * (0.5d0 * NL(ix,i) + 3.d0 * o2)
      reac(1,i,-1,1) = reac(1,i,1,1)
      reac(1,i, 1,2) = er0 * 0.25d0 * xdb * NL(ix,i) * v / (xtir1) * s1 
      reac(1,i,-1,2) = - reac(1,i,1,2)
      o3 = er0 * 0.5d0 * xdb * NL(ix,i) * s1
      reac(1,i, 1,3) = er0 * 2.d0 * x * o2 + o3
      reac(1,i,-1,3) = reac(1,i,1,3)
      reac(1,i, 1,3) = reac(1,i,1,3) - o3 * khi1
      reac(1,i,-1,3) = reac(1,i,-1,3) + o3 * khi1
      end do
        do ic = 1, 3
          do i = 0, maxcou
            do im = 2, 2*nmoant-1
            reac(im,i, 0,ic) = reac(1,i, 0,ic)
            reac(im,i, 1,ic) = reac(1,i, 1,ic)
            reac(im,i,-1,ic) = reac(1,i,-1,ic)
            end do
          end do
        end do
            
      do im = 1, 2*nmoant-1
C     >>>>>>>>>>>>>>>>>>>>>
c     p=0:
      u0i = (allkpa(im) * v) / omegag
      u0 = 1.d0 / u0i
      cmr0 = s0 - xdbi * (1.d0 - u0 * u0)
      smr02 = 1.d0 - cmr0 * cmr0
      csc22 = 2.d0 / (1.d0 - cmr0)
C       p=0 only absorbs if isig*k//>0, |u0i|>1, and x between...
c        if(isig .eq. dint(dsign(1.d0,u0i)) .and. dabs(u0i).gt.1.d0)then
        if(u0*fis.gt.petitb(ix) .and. u0*fis.lt.petita(ix))then
        smr0 = dsqrt(smr02)
        clmr0(1) = cmr0
c      r2 = r1 * xdbi * hp

        t3 = er0 * pi * xdbi * dabs(u0)**3 / smr0
        ac(im,0,0,2) = - v * t3
        ac(im,0,0,3) =  2.d0 * x * t3
        ac(im,1,0,2) = cmr0 * ac(im,0,0,2)
        ac(im,1,0,3) = cmr0 * ac(im,0,0,3)
          do i = 2, maxcou
          clmr0(i) = 2.d0 * cmr0 * clmr0(i-1) - clmr0(i-2)
          ac(im,i,0,2) = ac(im,0,0,2) * clmr0(i)
          ac(im,i,0,3) = ac(im,0,0,3) * clmr0(i)
          end do
        end if

c     Reactive, p=0:
C     Delta L,0 term:
      reac(im,0,0,2) = reac(im,0,0,2) + er0 * fis * 0.5d0 * v * PI * u0
      reac(im,0,0,3) = reac(im,0,0,3) - er0 * fis * PI * u0 * x
C     LL and JL sequences for u0:
      LL(0,1) = 0.d0
c       Passing:
        if(pass)then
        alp2 = kelv(ix)**2 / (1.d0 - (u0/petita(ix))**2)
        JL(0,1) = u0 / petita(ix) * alp2 * ellpic(-alp2, kelv(ix))
        if(alp2.lt.1.d0)JL(0,1) = JL(0,1) + fis * 0.5*pi*alp2/dsqrt(1.d0-alp2)
c        JL(0,1) = - 4.d0 / (petita(ix) + petitb(ix)) * 
c     ;  (xdb * ek(ix) / (petita(ix) + u0) - csc22 * petita(ix) * 
c     ;   ellpic(-kelv(ix) * (petita(ix) + u0) / (petita(ix) - u0), kelv(ix)) )
        JL(1,1) = fis * PI
     ;          + cmr0 * JL(0,1) + u0 * NL(ix,0)
        LL(1,1) = 0.5d0 * (smr02 * JL(0,1) - fis * cmr0 * PI) - u0 * (
     ;  ek(ix)*(cmr0+cmmv(ix))+2.d0*ee(ix)/kelv(ix)**2) / petita(ix)
          do i = 1, maxcou-1
          LL(i+1,1) = - LL(i-1,1) + 2.d0 * cmr0 * LL(i,1) - 0.5d0 * u0 *
     ;    (NL(ix,i+1) - NL(ix,i-1))
          JL(i+1,1) = JL(i-1,1) - 4.d0 * LL(i,1)
          end do
        else
c@here
c       Trapped:
        JL(0,1) = 2.d0 * u0 * csc22 * ellpic(-csc22*kelv(ix)**2 ,kelv(ix))
        JL(1,1) = cmr0 * JL(0,1) + u0 * NL(ix,0)
        LL(1,1) = smr02 * 0.5d0 * JL(0,1) - u0 * sqrt2 * sxdbi * (2.d0*ee(ix) - ek(ix) / csc22)
          do i = 1, maxcou-1
          LL(i+1,1) = - LL(i-1,1) + 2.d0 * cmr0 * LL(i,1) - u0 *
     ;    (NL(ix,i+1) - NL(ix,i-1))
          JL(i+1,1) = JL(i-1,1) - 4.d0 * LL(i,1)
          end do
        end if
        o4 = er0 * u0 * u0
        do i = 0, maxcou
        o5 = NL(ix,i) + u0 * sxdb * JL(i,1)
        reac(im,i,0,2) = reac(im,i,0,2) + o4 * v * 0.5d0 * o5
        reac(im,i,0,3) = reac(im,i,0,3) - o4 * x * o5
        end do

c     p=±1: 2 roots to resonance equations.
      do ip = -1, 1, 2
c     ````````````````
      khi = dfloat(ip) * khi1
      lamb = 0.5 * u0i * khi
      del2 = lamb * lamb + 1.d0 - khi
      cxr = del2 .lt. 0.d0

        if(cxr)then
        del = dsqrt(- del2)
        ru1 = lamb
        ru2 = lamb
        iu1 = del
        iu2 = - del
        cu1 = dcmplx(lamb, del)
        t1 = dfloat(ip) / (qom(ispe) * delb(intab))
        t2 = omegag - dfloat(ip) * qom(ispe) * bbar(intab)
        ccmr1 = t1 * (allkpa(im) * v * cu1 - t2)
        csmr12 = 1.d0 - ccmr1 * ccmr1

c       reac(im,0,-1:1)
c       JL sequence for complex u1:
      LLc(0) = 0.d0
      LLc(0) = 0.d0
          if(pass)then
          t3 = - 4.d0 * xdb / (petita(ix) + petitb(ix))
          JLc(0) = t3 * (
     ;ek(ix) / (cu1 + fis*petita(ix))
     ;+ 2.d0 * fis * petita(ix) / (cu1**2 - petita(ix)**2)
     ;* cdellpic(-kelv(ix)*(petita(ix)+fis*cu1)/(petita(ix)-fis*cu1), kelv(ix))
     ;)
          JLc(1) = fis * pi
c finir:
          LLc(1) = 0.5 * fis * pi * ccmr1
          else
          t3 = 2.d0 / (1.d0 - ccmr1)
          JLc(0) = 2.d0 * cu1 * t3 * 
     ;    cdellpic(dcmplx(-t3*kelv(ix)**2,0.d0), kelv(ix))
          JLc(1) = 0.d0
c finir:
          LLc(1) = 0.d0
          end if
        JLc(1) = JLc(1) + ccmr1 * JLc(0) + cu1 * NL(ix,0)
        JLc(-1) = JLc(1)
        LLc(1) = LLc(1) + 0.5 * csmr12 * JLc(0)
        LLc(-1) = - LLc(1)
            do i = 2, maxcou
            LLc(i) = - LLc(i-2) + 2.d0 * ccmr1 * LLc(i-1)
            JLc(i) = JLc(i-2) + 4.d0 * LLc(i-1)
            end do
          t3 = er0 * s1 / (dfloat(ip)*xtir1 * 4.d0 * del)
          t4 = v * t3
          t5 = t3 * 2.d0 * (x - dfloat(ip) * xtir1)
            do i = 0, maxcou
            o5 = dimag((1-u1*u1)*JLc(i))
            reac(im,i,ip,2) = reac(im,i,ip,2) - t4 * o5
            reac(im,i,ip,3) = reac(im,i,ip,3) + t5 * o5
            end do

        else
        del = dsqrt(del2)
        u1 = lamb + del
        u2 = lamb - del
        t1 = dfloat(ip) / (qom(ispe) * delb(intab))
        t2 = omegag - dfloat(ip) * qom(ispe) * bbar(intab)
        cmr1 = t1 * (allkpa(im) * v * u1 - t2)
        cmr2 = t1 * (allkpa(im) * v * u2 - t2)
c       ac(im,0,-1:1)
c see: two roots can take place together
          if(dabs(cmr1) .le. 1.d0)then
          smr12 = 1.d0 - cmr1 * cmr1
          smr1 = dsqrt(smr12)
          clmr1(1) = cmr1
          t3 = r3 * (1.d0 - u1 * u1) / (del * smr1)
          t4 = v * t3
          t5 = t3 * 2.d0 * (x - ip * xtir1)
          ac(im,0,ip,2) = ac(im,0,ip,2) - t4
          ac(im,1,ip,2) = ac(im,1,ip,2) - t4 * cmr1
          ac(im,0,ip,3) = ac(im,0,ip,3) + t5
          ac(im,1,ip,3) = ac(im,1,ip,3) + t5 * cmr1
            do i = 2, maxcou
            clmr1(i) = 2.d0 * cmr1 * clmr1(i-1) - clmr1(i-2)
            ac(im,i,ip,2) = ac(im,i,ip,2) - t4 * clmr1(i)
            ac(im,i,ip,3) = ac(im,i,ip,3) + t5 * clmr1(i)
            end do
          end if
          if(dabs(cmr2).le.1.d0)then
          smr22 = 1.d0-cmr2*cmr2
          smr2 = dsqrt(smr22)
          clmr2(1) = cmr2
          t3 = r3 * (1-u2*u2) / (del*smr2)
          t4 = v * t3
          t5 = t3 * 2.d0 * (x - ip * xtir1)
          ac(im,0,ip,2) = ac(im,0,ip,2) - t4
          ac(im,1,ip,2) = ac(im,1,ip,2) - t4 * cmr2
          ac(im,0,ip,3) = ac(im,0,ip,3) + t5
          ac(im,1,ip,3) = ac(im,1,ip,3) + t5 * cmr2
            do i = 2, maxcou
            clmr2(i) = 2.d0 * cmr2 * clmr2(i-1) - clmr2(i-2)
            ac(im,i,ip,2) = ac(im,i,ip,2) - t4 * clmr2(i)
            ac(im,i,ip,3) = ac(im,i,ip,3) + t5 * clmr2(i)
            end do
          end if
c       reac(im,0,-1:1)
c       JL sequences for u1 and u2:
      LL(0,1) = 0.d0
      LL(0,2) = 0.d0
          if(pass)then
          t3 = - 4.d0 * xdb / (petita(ix) + petitb(ix))
          JL(0,1) = t3 * (
     ;ek(ix) / (u1 + fis*petita(ix))
     ;+ 2.d0 * fis * petita(ix) / (u1**2 - petita(ix)**2)
     ;* ellpic(-kelv(ix)*(petita(ix)+fis*u1)/(petita(ix)-fis*u1), kelv(ix))
     ;)
          JL(0,2) = t3 * (
     ;ek(ix) / (u2 + fis*petita(ix))
     ;+ 2.d0 * fis * petita(ix) / (u2**2 - petita(ix)**2)
     ;* ellpic(-kelv(ix)*(petita(ix)+fis*u2)/(petita(ix)-fis*u2), kelv(ix))
     ;)
          JL(1,1) = fis * pi
          JL(1,2) = fis * pi
c finir:
          LL(1,1) = 0.5 * fis * pi * cmr1
          LL(1,2) = 0.5 * fis * pi * cmr2
          else
          t3 = 2.d0 / (1.d0 - cmr1)
          JL(0,1) = 2.d0 * u1 * t3 * ellpic(-t3*kelv(ix)**2,kelv(ix))
          t4 = 2.d0 / (1.d0 - cmr2)
          JL(0,2) = 2.d0 * u2 * t4 * ellpic(-t4*kelv(ix)**2,kelv(ix))
          JL(1,1) = 0.d0
          JL(1,2) = 0.d0
c finir:
          LL(1,1) = 0.d0
          LL(1,2) = 0.d0
          end if
        JL(1,1) = JL(1,1) + cmr1 * JL(0,1) + u1 * NL(ix,0)
        JL(1,2) = JL(1,2) + cmr2 * JL(0,2) + u2 * NL(ix,0)
        JL(-1,1) = JL(1,1)
        JL(-1,2) = JL(1,2)
        LL(1,1) = LL(1,1) + 0.5 * smr12 * JL(0,1)
        LL(1,2) = LL(1,2) + 0.5 * smr22 * JL(0,2)
        LL(-1,1) = - LL(1,1)
        LL(-1,2) = - LL(1,2)
            do i = 2, maxcou
            LL(i,1) = - LL(i-2,1) + 2.d0 * cmr1 * LL(i-1,1)
            LL(i,2) = - LL(i-2,2) + 2.d0 * cmr2 * LL(i-1,2)
            JL(i,1) = JL(i-2,1) + 4.d0 * LL(i-1,1)
            JL(i,2) = JL(i-2,2) + 4.d0 * LL(i-1,2)
            end do
          t3 = er0 * s1 / (dfloat(ip)*xtir1) * 0.125 / del
          t4 = v * t3
          t5 = t3 * 2.d0 * (x - dfloat(ip) * xtir1)
            do i = 0, maxcou
            o5 = (1-u1*u1)*JL(i,1) - (1-u2*u2)*JL(i,2)
            reac(im,i,ip,2) = reac(im,i,ip,2) - t4 * o5
            reac(im,i,ip,3) = reac(im,i,ip,3) + t5 * o5
            end do
        end if
      end do
c     ``````
      end do
C     >>>>>>

c     Loops over 'top right' nodes of surrounding elements:
c     Implements analytical integration of products of basis functions times linear interpolation of remaining coefficients
c     The loops over j and i sum contributions of surrounding elements times dielectric coefficients at the current mesh node.
        do j = max(0, ifiav(irv)-jev+1), min(1, ilaav(irv)-jev)
        ipov = 1 - j
        ji1 = jev + j
        jm1 = ji1 - 1
        do i = max(0, ifiax(irx)-iex+1), min(1, ilaax(irx)-iex)
        ipox = 1 - i
        i1 = iex + i
        im1 = i1 - 1
      
        bn(1,1) = F0i(im1, jm1, 1, isp)
        bn(2,1) = F0i(im1, jm1, 2, isp)
        bn(1,2) = F0i(im1, jm1, 3, isp)
        bn(2,2) = F0i(im1, jm1, 4, isp)
        bn(1,3) = F0i(i1, jm1, 1, isp)
        bn(2,3) = F0i(i1, jm1, 2, isp)
        bn(1,4) = F0i(i1, jm1, 3, isp)
        bn(2,4) = F0i(i1, jm1, 4, isp)
        bn(3,1) = F0i(im1, ji1, 1, isp)
        bn(4,1) = F0i(im1, ji1, 2, isp)
        bn(3,2) = F0i(im1, ji1, 3, isp)
        bn(4,2) = F0i(im1, ji1, 4, isp)
        bn(3,3) = F0i(i1, ji1, 1, isp)
        bn(4,3) = F0i(i1, ji1, 2, isp)
        bn(3,4) = F0i(i1, ji1, 3, isp)
        bn(4,4) = F0i(i1, ji1, 4, isp)
c---------------------------------------------------
cPL 13/5/05: item to check: when current node is at a region boundary, may find elements with a length different from 
c the value in current region!
c---------------------------------------------------
c       Loop over v basis functions
        do ib = 1, 4
c         Loop over x basis functions
          do jb = 1, 4
            do ip = -1, 1
            do ih = 0, maxcou
            do im = 1, 2*nmoant-1

c A few explanations (PL26/5/2005):

c im: index of average poloidal mode; ih: corresponds to l index of theory; ip: cyclotron harmonic (-1,0,1)

c Last index of ac and reac: 1, 2 or 3, correspond to coefficient of f0, df0/dv, df0/dx respectively.

c bn array: has received coefficients of bicubic basis function representation of the distribution function on current element
c

c prods array:
c 1st index v basis function
c 2nd index v derivative index, 0 or 1
c 3rd index labels function of reduced variable ksi to integrate over v element: 0 for (1-ksi), 1 for ksi
c 4th index x basis function
c 5th index x derivative index, 0 or 1
c 6th index labels function of reduced variable ksi to integrate over x element: 0 for (1-ksi), 1 for ksi

c I confirm a factor * lv * lx is missing in these expressions / pull out of nested loops as much as possible!

c           Coeff. of f0, df0/dv, df0/dx:
            tact(im, ih, ip) = tact(im, ih, ip) + bn(ib,jb) * (
     ;        ac(im,ih,ip,2) * prods(ib,1,ipov,jb,0,ipox)
     ;      + ac(im,ih,ip,3) * prods(ib,0,ipov,jb,1,ipox)    )
            treact(im, ih, ip) = treact(im, i, ip) + bn(ib,jb) * (
     ;        reac(im,ih,ip,1) * prods(ib,0,ipov,jb,0,ipox)
     ;      + reac(im,ih,ip,2) * prods(ib,1,ipov,jb,0,ipox)
     ;      + reac(im,ih,ip,3) * prods(ib,0,ipov,jb,1,ipox)      )
            end do
            end do
            end do
          end do
        end do

        end do
        end do

      end if

        end do
C       ~~~~~~
      end do
C     ======
      end do     
c     ------
      end do     
c     ======
      end do    
      end do    

c     In the end: convolution with geom. coeffs...
c     Assembly of contributions: add to vmat, same assembly
c     as Maxwellian contributions will follow in ASPLAS.

      iplb = 0
        DO MAV2 = 2*MINF(IEL), 2*MSUP(IEL)
        im = mav2 - 2*MINF(IEL) + 1
          IF(MAV2 .LE. MINF(IEL)+MSUP(IEL))THEN
          KSTOP2 = MIN(KLIM, MAV2-2*MINF(IEL))
          ELSE
          KSTOP2 = MIN(KLIM, 2*MSUP(IEL)-MAV2)
          END IF
        IF(MOD(MAV2-KSTOP2,2).NE.0)KSTOP2 = KSTOP2 - 1
        KSTOP1 = - KSTOP2
 
          DO K = KSTOP1, KSTOP2, 2
          IPLB = IPLB + 1
c@ wrong array bounds!
            do L = k - ncresp, k + ncresp
            La = iabs(L)
            vmat(1,1,iplb) = vmat(1,1,iplb) + gcdr(L-k,k) * 
     ;      dcmplx(act(im,La,1), react(im,La,1))
            vmat(3,3,iplb) = vmat(3,3,iplb) + gcdr(L-k,k) * 
     ;      dcmplx(act(im,La,-1), react(im,La,-1))
            vmat(5,5,iplb) = vmat(5,5,iplb) + gcdr(L-k,k) * 
     ;      dcmplx(act(im,La,0), react(im,La,0))
            end do
          END DO
        END DO

      return
      end