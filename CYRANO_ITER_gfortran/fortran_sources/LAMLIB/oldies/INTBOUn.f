      subroutine intboun(side, lblo)

      implicit none

c     Side of boundary ('L' or 'R'):
      character*1 side
      integer lblo
      
c     General constraint at boundary between 2 media.
c     This version with loops over average poloidal mode index.

c     To see: general diel response contrib. to [Erho] and FLR!
c     Currently treated as Maxwellian here!

c 12/12/03: contribution to right-hand side of linear system
c from jump conditions is now assembled here

@@@@@@@@@@@@ not done yet!!!!!!!!!!!!!!
      
      include 'pardim.copy'
      include 'comusr.f'
      include 'comequ.copy'
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
      include 'comphy.copy'

      logical trafo, rtp

      character*3  eletyp
      
      integer 
     ;  i, j, i2, index(6), icolo, ibulo, isigne, iregb, iell, ielb, ielr
     ;, idcato, jdcato, minfc, msupc, nmodc, iregne, ipath, impk
     ;, idcat, jdcat, kstop1, kstop2, kp, im, idum
     ;, ini
c     ;, nfft
     
      double precision tabaf(npfft,4)
      
      complex*16 
     ;  vl(6,6), vls(6,6), v2(0:npfft,19), ztra(6,6)
     ;, zhf(5,5), zhf2(5,6), fac, signe
     ;, work(3*npft2+2), foutra(-npft2:npft2,4)
     ;, cmat(5,6), cmati(6,5)

      external zset
      
      save ini, cmat, cmati, fac
      data ini/0/

        if(ini.eq.0)then
c       5 by 6 matrices for conversions to/from +,-,//:
        call zset(5*6, czero, cmat, 1)
        call zset(6*5, czero, cmati, 1)
        cmat(1,1) = sqrt2i
        cmat(1,3) = sqrt2i
        cmat(2,1) = dcmplx(0.d0,-sqrt2i)
        cmat(2,3) = dcmplx(0.d0, sqrt2i)
        cmat(3,2) = dcmplx(0.d0,-sqrt2i)
        cmat(3,4) = dcmplx(0.d0, sqrt2i)
        cmat(4,5) = cun
        cmat(5,6) = cun
        
        cmati(1,1) = sqrt2i
        cmati(3,1) = sqrt2i
        cmati(1,2) = dcmplx(0.d0, sqrt2i)
        cmati(3,2) = dcmplx(0.d0,-sqrt2i)
        cmati(2,3) = dcmplx(0.d0, sqrt2i)
        cmati(4,3) = dcmplx(0.d0,-sqrt2i)
        cmati(5,4) = cun
        cmati(6,5) = cun
        fac = dcmplx(0.d0, -1.d0 / (omegag * mu0 * rnorm))
        ini = 1
        end if
  
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      eletyp = styp(isubr)
 
        if(side .eq. 'L' .or. side .eq. 'l')then
        isigne = - 1
        signe = - cun
        iregb = ireg
        iell = iel
        ielb = iel
        ielr = min0(iel + 1, nele)
        idcato = lblo
        jdcato = 0
        intab = iell * (ngauss + 1) + iregb
         
        else if(side .eq. 'R' .or. side .eq. 'r')then
        isigne = 1
        signe = cun
        iregb = ireg - 1
        iell = max0(1, iel - 1)
        ielb = iel - 1
        ielr = iel
        idcato = 0
        jdcato = lblo
        intab = iell * (ngauss + 1) + iregb + 1
        end if
 
      minfc = minf(ielb)
      msupc = msup(ielb)
      nmodc = nmode(ielb)
      y = rx0(iregb)
      ig = 1
      yinv = abscni(intab)

      coljum = coldpl(ireg) .or. vacuum(ireg) .or. .not.flrops(ireg)
      iregne = ireg - isigne
      if(iregne.ge.1 .and. iregne.le.nreg)coljum = coljum .and.
     ;( coldpl(iregne) .or. vacuum(iregne) .or. .not.flrops(iregne) )

c     Disables 'hot' jump condition if some species are RK's nonmaxwellian:
      if(.not.glomax)coljum = .true.

      if(iregb.gt.0 .and. iregb.lt.nreg)then
c     internal boundaries:
        do i = 1, ncstr(iregb)
        index(i) = i
        end do
      else
c     Wall: constraints on Etheta, Ephi; also zero radial particle
c     current for hot plasma.
        if(coljum)then
        index(1) = 2
        index(2) = 4
        else
        index(1) = 1
        index(2) = 2
        index(3) = 4
        end if
      end if

      if(circ)then
c     ============
      call fouco
      call loctri(iel, .false.)
      call locti2
      call loctr2
  
c     Plasma contributions:
c     --------------------
C     Remains to see: metal walls with hot plasma!:
      if(   .not. vacuum(ireg)
     ;.and. .not.(coljum .and. rbtyp(iregb).eq.'MET')
     ;.and. (coljum .or. (.not.coldpl(ireg) .and. flrops(ireg)))
     ; )then
 
c     Switch indicating the need for poloidal ffts:
      trafo = .not. polsym
      ipath = 1
c   new:
c   see with general diel response!
      do 2 mav2 = 2*minf(ielb), 2*msup(ielb)
cPL      call polfft(v2(0,1), npfft+1, .false., idum, .true., 2, trafo, ipath, nfft
      call polfft(v2(0,1), npfft+1, .false., idum, .true., 2, trafo, ipath
     ;, onlyab_now)
        if(trafo)then
          if(mav2 .le. minf(ielb) + msup(ielb))then
          kstop2 = min(klim, mav2 - 2 * minf(ielb))
          else
          kstop2 = min(klim, 2 * msup(ielb) - mav2)
          end if
        else
        kstop2 = 0
        end if
      if(mod(mav2-kstop2, 2) .ne. 0)kstop2 = kstop2 - 1
      kstop1 = - kstop2
 
      do 2 k = kstop1, kstop2, 2
      mpk = (mav2 + k) / 2
      m = mpk - k
      impk = mpk + 1 - minf(ielb)
      jdcat = jdcato + icolo * (impk - 1)
      kp = iabs(k)
      im = m + 1 - minfc
      idcat = idcato + ncstr(iregb) * (im - 1)
 
cPL fixed bug 26/5/04          call plasmb(vl, vls, k, k, v2(0,1), 0, npfft, 2, .false.)
      call plasmb(vl, vls, k, k, v2(0,1), 0, npfft+1, 2, .false.)
cPL      call plasmb(vl, vls, k, k, v2(0,1), 0, npfft, 2, nfft, .false.)
c     constraints expressed in r,t,p coord.:
        if(coljum .or. rbtyp(iregb).eq.'MET')then
        i2 = 1
c       bc2ih 5*6, rectangular matrices; bc2ih has 2 indices.
        call mul2c3_fast(ztra, bc2ih, vl, 5, 6, 6, 6, 5, 6)
        else
        i2 = ncstr(iregb)
c       bcih 6*6*(0;...), square matrices; mind: local value stored in 3rd arg
c        = 1, but index starts at 0!
        call mul2c_fast(ztra, bcih(1,1,1), vl, 6)
c check normalization of ztra
        end if

      if(eletyp.eq.'M23')call mul2c3_fast(ztra, ztra, bc2i, 5, 6, 5, 6, 6, 6)

c !! M23 AND HOT PLASMA NOT CHECKED
c Voir these eq. 4.51 - 4.53: rows 1, 3, 5 only have nonzero contributions!
        do i = 1, i2, 2
          do j = 1, icolo
          ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j) + signe * ztra(i,j)
          end do
        end do
    2 continue
      end if

c     Vacuum terms: 
c     ------------
      do im = 1, nmode(ielb)
      m = im - 1 + minfc
      idcat = idcato + ncstr(iregb) * (im - 1)
      kstop1 = - min0(ncrot, m - minf(ielb))
      kstop2 =   min0(ncrot, msup(ielb) - m)

        do k = kstop1, kstop2
        kp = iabs(k)
        mpk = m + k
        impk = mpk + 1 - minf(ielb)
        jdcat = jdcato + icolo * (impk - 1)
        call zset(25, czero, zhf, 1)

c    The normalization below gives the following constraints:
c    [Erho + ...(plasma)] = ...                    (charge / epsilon0)
c    [Etheta] = ...                                (0)
c    [dEtheta/dy -dErho/(y dtheta)] = ...          (-i omega mu0 rnorm jtheta)
c    [Ephi] = ...                                  (0)
c    (R/R0)*[dEphi/dy -rnorm*dErho/(R dphi)] = ... (-i omega mu0 rnorm Rjphi/R0)

          if(k .eq. 0)then
          if(coljum)zhf(1,1) = cun
          zhf(2,2) = cun
          zhf(3,1) = - ci * mpk * yinv
          zhf(3,3) = cun
          zhf(4,4) = cun
          zhf(5,1) = - ci * kprn
          end if

        zhf(5,5) = afou(kp)
     
          if(eletyp.eq.'HEC')then
          call mul2c3_fast(zhf2, zhf, bc2, 5, 5, 6, 5, 5, 5)
            do i = 1, ncstr(iregb)
              do j = 1, icolo
              ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j)
     ;                               + signe * zhf2(index(i),j)
              end do
            end do
          else if(eletyp.eq.'M23')then
            do i = 1, ncstr(iregb)
              do j = 1, icolo
              ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j)
     ;                               + signe * zhf(index(i),j)
              end do
            end do
          end if
        end do

      end do

      else
c     ====
c     D shape or general equilibrium:
c     Select which components to use for jump conditions:
      rtp = styp(isubr).eq.'M23'
c     When two sides of boundary use different element types, give priority 
c     to r,t,p:
      if(iregb.ne.nreg)rtp = rtp .or. styp(isubr-isigne).eq.'M23'

c     Vacuum terms: 
c     generates CVC(K,L), K=-npfft/2,npfft/2
c     requires npfft/2 >= ncrot (checked in main)
      call foucoj(rtp)

      if(eletyp.eq.'M23')then
C     +++++++++++++++++++++++  

      do im = 1, nmode(ielb)
      m = im - 1 + minfc
      idcat = idcato + ncstr(iregb) * (im - 1)
      kstop1 = - min0(ncrot, m-minf(ielb))
      kstop2 =   min0(ncrot, msup(ielb)-m)

        do k = kstop1, kstop2
        kp = iabs(k)
        mpk = m + k
        impk = mpk + 1 - minf(ielb)
        jdcat = jdcato + icolo * (impk - 1)
        call zset(5*5, czero, zhf, 1)
     
c    The normalization below gives the following constraints:
c    [Erho + ...(plasma)] = ...                    (charge / epsilon0)
c    [Etheta] = ...                                (0)
c    [dEtheta/dy -dErho/(y dtheta)] = ...          (-i omega mu0 rnorm jtheta)
c    [Ephi] = ...                                  (0)
c    (R/R0)*[dEphi/dy -rnorm*dErho/(R dphi)] = ... (-i omega mu0 rnorm Rjphi/R0)

          if(k .eq. 0)then
          if(coljum)zhf(1,1) = cun
          zhf(2,2) = cun
          zhf(3,3) = cun
          zhf(4,4) = cun
          zhf(5,5) = cun
          end if

        zhf(3,1) = cvc(k,2) + ci * mpk * cvc(k,3)
        zhf(5,1) = cvc(k,1)
          do i = 1, ncstr(iregb)
            do j = 1, icolo
            ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j)
     ;      + signe * zhf(index(i),j)
            end do
          end do
        end do

      end do

      else if(eletyp.eq.'HEC' .and. rtp)then
c     ++++++++++++++++++++++++++++++++++++++
c     Constraints are in rtp, but unknowns are in +-// !
      do im = 1, nmode(ielb)
      m = im - 1 + minfc
      idcat = idcato + ncstr(iregb) * (im - 1)
      kstop1 = - min0(ncrot, m - minf(ielb))
      kstop2 =   min0(ncrot, msup(ielb) - m)

        do k = kstop1, kstop2
        kp = iabs(k)
        mpk = m + k
        impk = mpk + 1 - minf(ielb)
        jdcat = jdcato + icolo * (impk - 1)
        call zset(5*5, czero, zhf, 1)
     
          if(k .eq. 0)then
          if(coljum)zhf(1,1) = cun
          zhf(2,2) = cun
          zhf(4,4) = cun
          end if
c       [Hphi] eq.: Etheta'= cos * Eeta' + sin * E//':
        zhf(3,3) = cvc(k,6)
        zhf(3,5) = cvc(k,7)
c       [Htheta] eq.: Ephi'= - sin * Eeta' + cos * E//':
        zhf(5,3) = - zhf(3,5)
        zhf(5,5) = zhf(3,3)

        zhf(3,1) = cvc(k,2) + ci * mpk * cvc(k,3)
        zhf(5,1) = cvc(k,1)
c       Columns: rho,eta,// to +,-,//:
        call mul2c2_fast(zhf2, zhf, cmat, 5, 5, 6)

          do i = 1, ncstr(iregb)
            do j = 1, icolo
            ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j)
     ;      + signe * zhf2(index(i),j)
            end do
          end do
        end do

      end do

      else if(eletyp.eq.'HEC')then
C     ++++++++++++++++++++++++++++
c       jumps expressed in rho, eta, // ; check!
          do im = 1, nmode(ielb)
          m = im - 1 + minfc
          idcat = idcato + ncstr(iregb) * (im - 1)
          kstop1 = - min0(ncrot, m-minf(ielb))
          kstop2 =   min0(ncrot, msup(ielb)-m)
             
            do k = kstop1, kstop2
            kp = iabs(k)
            mpk = m + k
            impk = mpk + 1 - minf(ielb)
            jdcat = jdcato + icolo * (impk - 1)
            call zset(5*5, czero, zhf, 1)
              if(k .eq. 0)then
              if(coljum)zhf(1,1) = cun
              zhf(2,2) = cun
              zhf(3,3) = cun
              zhf(4,4) = cun
              zhf(5,5) = cun
              end if
            zhf(3,1) = cvc(k,3) + ci * mpk * cvc(k,4)
            zhf(5,1) = cvc(k,1) + ci * mpk * cvc(k,2)  
            call mul2c2_fast(zhf2, zhf, cmat, 5, 5, 6)

              do i = 1, ncstr(iregb)
                do j = 1, icolo
                ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j)
     ;          + signe * zhf2(index(i),j)
                end do
              end do
            end do
          end do

      end if
C     ++++++

c plasma: at work...
C     Remains to see: metal walls with hot plasma!:
      if(   .not. vacuum(ireg)
     ;.and. .not.(coljum .and. rbtyp(iregb).eq.'MET')
     ;.and. (coljum .or. (.not.coldpl(ireg) .and. flrops(ireg)))
     ; )then
 
      trafo = .not. polsym
      ipath = 1
c see with general diel response; now use contrib. as Maxwellian:
      do 3 mav2 = 2*minf(ielb), 2*msup(ielb)
cPL      call polfft(v2(0,1), npfft+1, .false., idum, .true., 2, trafo, ipath, nfft
      call polfft(v2(0,1), npfft+1, .false., idum, .true., 2, trafo, ipath
     ;, onlyab_now)
        if(trafo)then
          if(mav2 .le. minf(ielb) + msup(ielb))then
          kstop2 = min(klim, mav2 - 2 * minf(ielb))
          else
          kstop2 = min(klim, 2 * msup(ielb) - mav2)
          end if
        else
        kstop2 = 0
        end if
      if(mod(mav2 - kstop2, 2) .ne. 0)kstop2 = kstop2 - 1
      kstop1 = - kstop2
 
        do 3 k = kstop1, kstop2, 2
        mpk = (mav2 + k) / 2
        m = mpk - k
        impk = mpk + 1 - minf(ielb)
        jdcat = jdcato + icolo * (impk - 1)
        kp = iabs(k)
        im = m + 1 - minfc
        idcat = idcato + ncstr(iregb) * (im - 1)
         
cPL fixed bug 26/5/04          call plasmb(vl, vls, k, k, v2(0,1), 0, npfft, 2, .false.)
        call plasmb(vl, vls, k, k, v2(0,1), 0, npfft+1, 2, .false.)
cPL        call plasmb(vl, vls, k, k, v2(0,1), 0, npfft, 2, nfft, .false.)
c@
          if(eletyp.eq.'HEC')then
c         Constraints expressed in R,T,P coord, between unknowns in +-//:
c         Contribution to cond. [Drho]=0 (radial displacement vector is C0)
          ael(idcat+1,jdcat+1) = ael(idcat+1,jdcat+1) + signe * sqrt2i * vl(1,1)
          ael(idcat+1,jdcat+3) = ael(idcat+1,jdcat+3) + signe * sqrt2i * vl(3,3)
          else if(eletyp.eq.'M23')then
c         to do
          end if
    3 continue
      end if

      end if
c     ======

C     WRITE(6,*)((CSTR(I,J),J=1,6),I=1,NCSTR(IREGB))

c     Conjugate block (columns <-> Lagrange multipliers in global matrix):
      do i = jdcato + 1, jdcato + nmode(ielb) * icolo
        do j = idcato + 1, idcato + nmodc * ncstr(iregb)
        ael(i,j) = dconjg( ael(j,i) )
        end do
      end do

      return
      end
