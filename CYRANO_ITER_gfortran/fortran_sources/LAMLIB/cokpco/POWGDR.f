      subroutine powgdr(isp, pic10, plan, picm10)

      implicit none
      integer isp
      double precision pic10, plan, picm10

C     Absorption of a species, using equilibrium of QLFP code
C     interpolated at local radial index INTAB of Cyrano by routine INTF0
C     isp is index of species among the ones requiring general diel. response
C     Absolute species index is ISPE

      include 'pardim.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comgeo.copy'
      include 'comequ.copy'
      include 'comant.copy'
      include 'comfin.copy'
      include 'comfic.copy'
      include 'comin2.copy'
      include 'comrot.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'compla.copy'
      include 'comfou.copy'
      include 'commod.copy'
      include 'commmk.copy'
      include 'compri.copy'
      include 'comswe.copy'
      include 'comgdr.copy'
      include 'dynou2.copy'
      include 'comphy.copy'
 
      logical pass, wribb
      
      integer i, j, im, ix, jv, ip, iex, jev, nset0, isig, igv, igx
     ;, igv2, igx2, jev1, jev2, iex1, iex2
     ;, ispe, ii, irv, jvn, irx, ixn, icount, ic, ipov, ipox
     ;, j1, jm1, i1, im1, ib, jb
     ;, kstop1, kstop2, m1, m2, L, La
     ;, isrchfgt

C     act: im of iL      
      double precision
     ;  act(2*maxpom-1,-maxcou:maxcou,-1:1)
     ;, ac(2*maxpom-1,-maxcou:maxcou,-1:1,3), pac(-maxcou:maxcou,-1:1,3)
      equivalence (act,ac)
      double precision 
     ;  v, vsq, vi, xn, xna, x
      double precision 
     ;  u0, u0i, u1, u2, xtir1, khi
      double precision 
     ;  khi1, lamb, del2, del, te0, t1, t2, t3, t4, t5
     ;, cmr1, cmr2, smr1, smr12, smr2, smr22, cmr0
      double precision 
     ;  smr0, smr02, csc22
      double precision 
     ;  clmr0(0:maxcou), clmr1(0:maxcou), clmr2(0:maxcou)
     ;, er0, r00, r1, r3, s0, s1, s1i, s2, s2i, omxnai, xmax
     ;, omca, xdb, xdbi, sxdb, sxdbi
     ;, petita(5*maxnex), petitb(5*maxnex)
     ;, cmmv(5*maxnex)
     ;, kel
     ;, maomi
     ;, cmm
     ;, oo, fis
     ;, bn(4,4)
     ;, tacno(3), tacnop(-1:1,3), t

c     indices: mbar, m1-m2, p

      ispe = ispgdr(isp)
      wribb = .false.
        do i = 1, nraddr
        if(intab.eq.irafp(i))wribb = .true.
        end do
      pic10 = 0.d0
      plan = 0.d0

      nset0 = (2*maxpom-1)*(2*maxcou+1)*3
      call dset(3*nset0, 0.d0, ac, 1)      
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
     ; / (2.d0*mh*amass(ispe) * b0 * hachi(intab) * omegag)

      write(6,*)'powgdr: before passing'      
c     Passing:
      do ii = 1, 4, 3
        do i = igl(ii), idl(ii)
        xn = xngaug(i)
        xna = dabs(xn)
        petita(i) = dsqrt(xna)
        petitb(i) = dsqrt(maomi*(xna+xnsep))
        kel = (petita(i)-petitb(i)) / (petita(i)+petitb(i))
        if(kel.ne.0.d0)then
        cmm = - 0.5d0 * (kel + 1.d0 / kel)
        cmmv(i) = cmm
        else
c see:
        cmmv(i) = 0.d0
        end if
        end do
      end do

      write(6,*)'powgdr: before trapped'      
c     Trapped:
      do ii = 2, 3
        do i = igl(ii), idl(ii)
        xn = xngaug(i)
        xna = dabs(xn)
        omxnai = 1.d0 / (1.d0 - xna)
        cmm = (1.d0 - xna * s0) * omxnai
        cmmv(i) = cmm
        end do
      end do

c        if(wribb)then
c        write(unbblo,*)tacno(2), tacno(3)
c NB: write format: coeffs of df0/dv, df0/dx
c      do irv = 1, nvgreg
c      do irx = 1, nxgreg
c      do jev = jev1, jev2
c      do iex = iex1, iex2
c        end if

      write(6,*)'powgdr: before 2D loop'  
c     V region loop:
      do irv = 1, nvgreg
        if(gigdr)then
C       jev is an element index
        jev1 = ifielv(irv)
        jev2 = ilaelv(irv)
        else
C       jev is a node index; jv = jev
        jev1 = ifiav(irv)
        jev2 = ilaav(irv)
        end if

c     X region loop:
      do irx = 1, nxgreg

C     Element length dependent matrices (here for piecewise uniform meshes):
      call normat(vlelr(irv), xlelr(irx))
      call intnod

C       This assumes 4 x regions!:
        if(irx.le.2)then
        isig = 1
        else 
        isig = -1
        end if
      fis = dfloat(isig)
      pass = irx.eq.1 .or. irx.eq.4
        if(gigdr)then
C       iex is an element index
        iex1 = ifielx(irx)
        iex2 = ilaelx(irx)
        else
C       iex is a node index; ix = iex
        iex1 = ifiax(irx)
        iex2 = ilaax(irx)
        end if
      
c     V element loop:
c     when gigdr=.F., loop over nodes instead
      do jev = jev1, jev2
c     ===================
C     right node index (gigdr=.t.):
      jvn = irv + jev

c     X element loop:
c     when gigdr=.F., loop over nodes instead
      do iex = iex1, iex2
c     -------------------
C     right node index (gigdr=.t.):
      ixn = irx + iex
C     f0 and first derivatives at Gauss points:
      icount = 0
      if(gigdr)call fagp(jev, iex, jvn, ixn, isp)
      
      do igv = 1, igv2
C     ================
      jv = (jev - 1) * igv2 + igv
      v = vgaug(jv)
      vsq = v * v
      vi = 1.d0 / v
      er0 = r00 * vsq
      r1 = er0 * vsq * PI
      
        do igx = 1, igx2
C       ~~~~~~~~~~~~~~~~
        ix = (iex - 1) * igx2 + igx
        xn = xngaug(ix)
        x = (1.d0 - dabs(xn)) * xmax
        khi1 = x / xtir1
C       x*delta/B0:
        xdb = x / s1
        xdbi = 1.d0 / xdb
        sxdb = dsqrt(xdb)
        sxdbi = 1.d0 / sxdb
        
c     Code fragment for case of Gaussian integration, to update!
c     INCLUDE 'POGDRG.f'

      if(.not.gigdr)then
c     Semi analytical integration over elements
C     =========================================
c     Groupe selon coef. de f0, df0/dv, df0/dx (cf dernier indice ac et reac)

      r3 = er0 * PI * s1 * 0.25d0 / dabs(xtir1)
      call dset(3*nset0, 0.d0, ac, 1)      
      call dset(3, 0.d0, tacno, 1)      
      call dset(3*3, 0.d0, tacnop, 1)      
            
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

        t3 = er0 * PI * xdbi * dabs(u0)**3 / smr0
        pac(0,0,2) = - v * t3
        pac(0,0,3) =  2.d0 * x * t3
        pac(1,0,2) = cmr0 * pac(0,0,2)
        pac(1,0,3) = cmr0 * pac(0,0,3)
          do i = 2, maxcou
          clmr0(i) = 2.d0 * cmr0 * clmr0(i-1) - clmr0(i-2)
          pac(i,0,2) = pac(0,0,2) * clmr0(i)
          pac(i,0,3) = pac(0,0,3) * clmr0(i)
          end do
        end if

c     p=±1: 2 roots to resonance equations.
      do ip = -1, 1, 2
c     ````````````````
      khi = dfloat(ip) * khi1
      lamb = 0.5 * u0i * khi
      del2 = lamb * lamb + 1.d0 - khi

        if(del2 .ge. 0.d0)then
        del = dsqrt(del2)
        u1 = lamb + del
        u2 = lamb - del
        t1 = dfloat(ip) / (qom(ispe) * delb(intab))
        t2 = omegag - dfloat(ip) * qom(ispe) * bbar(intab)
        cmr1 = t1 * (allkpa(im) * v * u1 - t2)
        cmr2 = t1 * (allkpa(im) * v * u2 - t2)
        pac(0,ip,2) = 0.d0
        pac(1,ip,2) = 0.d0
        pac(0,ip,3) = 0.d0
        pac(1,ip,3) = 0.d0
c see if two roots can take place together:
          if(dabs(cmr1) .le. 1.d0)then
          smr12 = 1.d0 - cmr1 * cmr1
          smr1 = dsqrt(smr12)
          clmr1(1) = cmr1
          t3 = r3 * (1.d0 - u1 * u1) / (del * smr1)
          t4 = v * t3
          t5 = t3 * 2.d0 * (x - ip * xtir1)
          pac(0,ip,2) = pac(0,ip,2) - t4
          pac(1,ip,2) = pac(1,ip,2) - t4 * cmr1
          pac(0,ip,3) = pac(0,ip,3) + t5
          pac(1,ip,3) = pac(1,ip,3) + t5 * cmr1
            do i = 2, maxcou
            clmr1(i) = 2.d0 * cmr1 * clmr1(i-1) - clmr1(i-2)
            pac(i,ip,2) = pac(i,ip,2) - t4 * clmr1(i)
            pac(i,ip,3) = pac(i,ip,3) + t5 * clmr1(i)
            end do
          end if
          if(dabs(cmr2).le.1.d0)then
          smr22 = 1.d0-cmr2*cmr2
          smr2 = dsqrt(smr22)
          clmr2(1) = cmr2
          t3 = r3 * (1-u2*u2) / (del*smr2)
          t4 = v * t3
          t5 = t3 * 2.d0 * (x - ip * xtir1)
          pac(0,ip,2) = pac(0,ip,2) - t4
          pac(1,ip,2) = pac(1,ip,2) - t4 * cmr2
          pac(0,ip,3) = pac(0,ip,3) + t5
          pac(1,ip,3) = pac(1,ip,3) + t5 * cmr2
            do i = 2, maxcou
            clmr2(i) = 2.d0 * cmr2 * clmr2(i-1) - clmr2(i-2)
            pac(i,ip,2) = pac(i,ip,2) - t4 * clmr2(i)
            pac(i,ip,3) = pac(i,ip,3) + t5 * clmr2(i)
            end do
          end if
        end if
      end do
c     ``````

c     Convolution: see GENEQ case later on.
      mav2 = im + 2*minf(iel) - 1
        if(mav2 .le. minf(iel)+msup(iel))then
        kstop2 = min(klim, mav2-2*minf(iel))
        else
        kstop2 = min(klim, 2*msup(iel)-mav2)
        end if
      if(mod(mav2-kstop2,2).ne.0)kstop2 = kstop2 - 1
c        if(geneq)then
c        kstop1 = - kstop2
c        else
c       Check: FORTRAN mod function gives negative result when 1st arg < 0!!!
        kstop1 = mod(mav2, 2)
c        end if
 
c     Nodal values:
      do ip = -1, 1
c     3, 5, 1:
      ic = 5 - ip * (1 + 3 * ip)
        do k = kstop1, kstop2, 2
        m1 = (mav2 + k) / 2 + 1 - minf(ielm)
        m2 = m1 - k
          do L = k - ncresp, k + ncresp
          La = iabs(L)
            ac(im,k,ip,2) = ac(im,k,ip,2) + dreal(gcdr(L-k,k)) * pac(La,ip,2)
            ac(im,k,ip,3) = ac(im,k,ip,3) + dreal(gcdr(L-k,k)) * pac(La,ip,3)
          end do
        t = dreal(dconjg(xpmp(ic,intab,m2)) * xpmp(ic,intab,m1))
        if(k.ne.0)t = t * 2.d0
        ac(im,k,ip,2) = ac(im,k,ip,2) * t
        ac(im,k,ip,3) = ac(im,k,ip,3) * t
c        tacno(2) = tacno(2) + ac(im,k,ip,2)
c        tacno(3) = tacno(3) + ac(im,k,ip,3)
        tacnop(ip,2) = tacnop(ip,2) + ac(im,k,ip,2)
        tacnop(ip,3) = tacnop(ip,3) + ac(im,k,ip,3)
        end do
      end do

      end do
C     >>>>>>
      do ip = -1, 1
      tacno(2) = tacno(2) + tacnop(ip,2)
      tacno(3) = tacno(3) + tacnop(ip,3)
      end do
c     Building blocks for QLFP:
      if(wribb)then
      write(unbblo,*)tacno(2), tacno(3)
      end if

c     Loops over 'top right' nodes of surrounding elements:
        do j = max(0, ifiav(irv)-jev+1), min(1, ilaav(irv)-jev)
        ipov = 1 - j
        j1 = jev + j
        jm1 = j1 - 1
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
        bn(3,1) = F0i(im1, j1, 1, isp)
        bn(4,1) = F0i(im1, j1, 2, isp)
        bn(3,2) = F0i(im1, j1, 3, isp)
        bn(4,2) = F0i(im1, j1, 4, isp)
        bn(3,3) = F0i(i1, j1, 1, isp)
        bn(4,3) = F0i(i1, j1, 2, isp)
        bn(3,4) = F0i(i1, j1, 3, isp)
        bn(4,4) = F0i(i1, j1, 4, isp)
      
c       Loop over v basis functions
        do ib = 1, 4
c         Loop over x basis functions
          do jb = 1, 4
c         Coeff. of df0/dv, df0/dx:
          picm10 = picm10  + bn(ib,jb) * (
     ;      tacnop(-1,2) * prods(ib,1,ipov,jb,0,ipox)
     ;    + tacnop(-1,3) * prods(ib,0,ipov,jb,1,ipox)    )
          pic10 = pic10  + bn(ib,jb) * (
     ;      tacnop(1,2) * prods(ib,1,ipov,jb,0,ipox)
     ;    + tacnop(1,3) * prods(ib,0,ipov,jb,1,ipox)    )
          plan = plan  + bn(ib,jb) * (
     ;      tacnop(0,2) * prods(ib,1,ipov,jb,0,ipox)
     ;    + tacnop(0,3) * prods(ib,0,ipov,jb,1,ipox)    )
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

      return
      end