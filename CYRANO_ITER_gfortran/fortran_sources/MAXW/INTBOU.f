      subroutine intbou(side, lblo)

      implicit none

c     Side of boundary ('L' or 'R'):
      character*1 side
      integer lblo
      
c     General constraint at boundary between 2 media.
c     This version with loops over average poloidal mode index for plasma terms.

c     To see: general diel response contrib. to [Erho] and FLR!
c     Currently treated as Maxwellian here!
      
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
      include 'cokpco.copy'

      logical trafo, rtp

      character*3  eletyp
      
      integer 
     ;  i, j, i2, index(6), icolo, ibulo, isigne, iregb, iell, ielb, ielr
     ;, idcato, jdcato, minfc, msupc, nmodc, iregne, ipath, impk
     ;, idcat, jdcat, kstop1, kstop2, kp, im, idum
     ;, ini
c     ;, nfft
     
c      double precision tabaf(npfft,4)
      
      complex*16 
     ;  vl(6,6), vls(6,6), v2(0:npfft,28), ztra(6,6)
     ;, zhf(5,5), zhf2(5,6), fac, signe
     ;, cmat(5,6), cmati(6,5)
c     ;, work(3*npft2+2), foutra(-npft2:npft2,4)

      external zset
      
      save ini, cmat, cmati, fac
      data ini/0/

      write(nofile,*)'Enter intbou'

        if(ini.eq.0)then
C       5 by 6 matrices for conversions to/from +,-,//:
cERN        call zset(5*6, czero, cmat, 1)
cERN        call zset(6*5, czero, cmati, 1)
	  cmat = czero
	  cmati = czero

c       rho, eta, eta', //, //' to +, +', -, -', //, //':
        cmat(1,1) = dcmplx(sqrt2i, 0.d0)
        cmat(1,3) = dcmplx(sqrt2i, 0.d0)
        cmat(2,1) = dcmplx(0.d0,-sqrt2i)
        cmat(2,3) = dcmplx(0.d0, sqrt2i)
        cmat(3,2) = dcmplx(0.d0,-sqrt2i)
        cmat(3,4) = dcmplx(0.d0, sqrt2i)
        cmat(4,5) = cun
        cmat(5,6) = cun
        
        cmati(1,1) = dcmplx(sqrt2i, 0.d0)
        cmati(3,1) = dcmplx(sqrt2i, 0.d0)
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
c        intab = iell * (ngauss + 1) + iregb
c@2004 this will fail when gicurl, giplas false if radial table is reduced to nodes: / see other occurences
	  intab = ifiabs(ilael(ireg)) + ngauss + 1
        else if(side .eq. 'R' .or. side .eq. 'r')then
        isigne = 1
        signe = cun
        iregb = ireg - 1
        iell = max0(1, iel - 1)
        ielb = iel - 1
        ielr = iel
        idcato = 0
        jdcato = lblo
c        intab = iell * (ngauss + 1) + iregb + 1
	  intab = ifiabs(ifiel(ireg))
        end if
 
      minfc = minf(ielb)
      msupc = msup(ielb)
      nmodc = nmode(ielb)
      y = rx0(iregb)
      ig = 1
      yinv = abscni(intab)

      coljum = coldpl(ireg) .or. vacuum(ireg) .or. .not.flrops(ireg)
      iregne = ireg - isigne
      if(iregne.ge.1 .and. iregne.le.nreg)coljum = coljum .and. (coldpl(iregne) .or. vacuum(iregne) .or. .not.flrops(iregne))

c     Disables 'hot' jump condition if some species are RK's nonmaxwellian:
      if(.not.glomax)coljum = .true.

      if(iregb.gt.0 .and. iregb.lt.nreg)then
c     Internal boundaries:
        do i = 1, ncstr(iregb)
        index(i) = i
        end do
      else
c     Wall: constraints on Etheta, Ephi; also zero radial particle current for hot plasma.
c@2004 Has this an effect on magn axis?
        if(coljum)then
        index(1) = 2
        index(2) = 4
        else
        index(1) = 1
        index(2) = 2
        index(3) = 4
        end if
      end if

      if(circ .and. (.not.cokpco .or. vacuum(ireg)))then
c     ==================================================
cPL6Dec04: uses new version of fouco (covers cokpco circular case)
      call fouco
      call loctri(iel, .false.)
      call locti2
      call loctr2
  
      if(.not. vacuum(ireg)
     ;   .and. .not.(coljum .and. rbtyp(iregb).eq.'MET') .and. (coljum .or. (.not.coldpl(ireg) .and. flrops(ireg)))
     ;  )then
c     Eventual plasma contributions:
c     -----------------------------
c     Remains to see: metal walls with hot plasma!:
 
c     Switch indicating the need for poloidal ffts:
      trafo = .not. polsym
      ipath = 1
c   see with general diel response!
        do mav2 = 2*minf(ielb), 2*msup(ielb)
        call polfft(v2(0,1), npfft+1, .false., idum, .true., 2, trafo, ipath, onlyab_now)
cPL      call polfft(v2(0,1), npfft+1, .false., idum, .true., 2, trafo, ipath, nfft, onlyab)
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
 
          do k = kstop1, kstop2, 2
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
c           bc2ih 5*6, rectangular matrices; 
cERN(underflow) call mul2c3(ztra, bc2ih, vl, 5, 6, 6, 6, 5, 6)
			  call mul2c3_fast(ztra(1:5,1:6), bc2ih, vl, 5, 6, 6, 6, 5, 6)
            else
            i2 = ncstr(iregb)
c           bcih 6*6*(0;...), square matrices; mind: local value stored in 3rd arg
c            = 1, but index starts at 0!
cERN            call mul2c(ztra, bcih(1,1,1), vl, 6)
			call mul2c_fast(ztra, bcih(1,1,1), vl, 6)
c check normalization of ztra
            end if

          if(eletyp.eq.'M23')then
cERN	        call mul2c3(ztra, ztra, bc2i, 5, 6, 5, 6, 6, 6)
	    call mul2c3_fast(ztra, ztra, bc2i, 5, 6, 5, 6, 6, 6)
          else if(eletyp.eq.'CAR')then
	    print *, 'INTBOU: CAR element not written for circular case'
	    stop
	    end if

c !! M23 AND HOT PLASMA NOT CHECKED
c Voir these eq. 4.51 - 4.53: rows 1, 3, 5 only have nonzero contributions!
            do i = 1, i2, 2
              do j = 1, icolo
              ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j) + signe * ztra(i,j)
              end do  ! j
            end do  ! i
          end do  ! k
        end do  ! mav2
      end if ! Plasma contribution

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

c    The normalization below gives the following constraints:
c    [Erho + ...(plasma)] = ...                    (charge / epsilon0)
c    [Etheta] = ...                                (0)
c    [dEtheta/dy -dErho/(y dtheta)] = ...          (-i omega mu0 rnorm jtheta)   [this is +rnorm times phi component of curl E]
c    [Ephi] = ...                                  (0)
c    (R/R0)*[dEphi/dy -rnorm*dErho/(R dphi)] = ... (-i omega mu0 rnorm Rjphi/R0) [this is -R/R0*rnorm times theta component of curl E]

cPL        call zset(25, czero, zhf, 1)
        zhf = czero
          if(k .eq. 0)then
          if(coljum)zhf(1,1) = cun
          zhf(2,2) = cun
          zhf(3,1) = - ci * mpk * yinv
          zhf(3,3) = cun
          zhf(4,4) = cun
          zhf(5,1) = - ci * kprn
          end if

        zhf(5,5) = afou(kp)

cPL6Dec04: cokpco circular case, assuming dphibar/drho and dphibar/dtheta radially continuous::
c    [dEtheta/dy -dErho/(y dtheta)-dErho/(y dphi)*dphibar/dtheta] = ...      (-i omega mu0 rnorm jtheta)
          if(cokpco)then
	    zhf(3,1) = zhf(3,1) - ci * n * yinv * acokfou(k)
          end if
               
          if(eletyp.eq.'HEC')then
               call mul2c3(zhf2, zhf, bc2, 5, 5, 6, 5, 5, 5)
cJAC	    call mul2c3_fast(zhf2, zhf, bc2, 5, 5, 6, 5, 5, 5)
            do i = 1, ncstr(iregb)
              do j = 1, icolo
              ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j) + signe * zhf2(index(i),j)
              end do
            end do
          else if(eletyp.eq.'M23')then
            do i = 1, ncstr(iregb)
              do j = 1, icolo
              ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j) + signe * zhf(index(i),j)
              end do
            end do
          else if(eletyp.eq.'CAR')then
	    print *, 'INTBOU: CAR element not written for circular case'
          stop
          end if  ! eletyp
        end do  ! k

      end do  ! im

      else
c     ====
c     D shape or general equilibrium, or const. k// coordinates (noncircular or circular with plasma):
c     Select which components to use for jump conditions:
      rtp = styp(isubr).eq.'M23' .or. styp(isubr).eq.'CAR'
c     When two sides of boundary use different element types, give priority to r,t,p:
      if(iregb.ne.nreg)rtp = rtp .or. styp(isubr-isigne).eq.'M23' 
     ;                           .or. styp(isubr-isigne).eq.'CAR'
c
c     ccccccccccccccccc
c     c Vacuum terms: c
c     ccccccccccccccccc
c  
c     generates CVC(K,L), K=-npfft/2,npfft/2
c     requires npfft/2 >= ncrot (checked in main)
      call foucoj(rtp)

      if(eletyp.eq.'M23')then
c     +++++++++++++++++++++++  

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
        end do  ! k

      end do  ! im

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
c 2004  zhf: rows are in r,t,p, columns in rho, eta, //     
          if(k .eq. 0)then
          if(coljum)zhf(1,1) = cun
cPL 5/4/2004: bug: columns of zhf must be in rho, eta, //:
c          zhf(2,2) = cun
c          zhf(4,4) = cun
          end if
cPL 5/4/2004: bug: columns of zhf in rho, eta, //:
c       [Etheta] eq.: Etheta =   cos * Eeta + sin * E//:
        zhf(2,2) = cvc(k,6)
        zhf(2,4) = cvc(k,7)
c       [Ephi] eq.:   Ephi   = - sin * Eeta + cos * E//:
        zhf(4,2) = - zhf(2,4)
        zhf(4,4) = zhf(2,2)
cPL ...already better

c       [Hphi] eq.: Etheta'= cos * Eeta' + sin * E//':
c NB terms including angle radial derivative are missing
cERN  cvc(k,8) = d/dr(cos THETA) and cvc(k,9) = d/dr(sin THETA) 
c     ADDED in FOUCOJ. Need to update below!!

        zhf(3,3) = cvc(k,6)
        zhf(3,5) = cvc(k,7)
c       [Htheta] eq.: Ephi'= - sin * Eeta' + cos * E//':
c NB terms including angle radial derivative are missing
        zhf(5,3) = - zhf(3,5)
        zhf(5,5) = zhf(3,3)

        zhf(3,1) = cvc(k,2) + ci * mpk * cvc(k,3)
        zhf(5,1) = cvc(k,1)

c        print*, zhf(:,:)
c        print*, signe

c       Convert columns of zhf from rho,eta,// to +,-,// in zhf2:
        call mul2c2(zhf2, zhf, cmat, 5, 5, 6)
cJAC		call mul2c2_fast(zhf2(1:5,1:6), zhf, cmat, 5, 5, 6)
c        print*, zhf2(1:5,:1)
c        print*, zhf2(1:5,:2)
c        print*, zhf2(1:5,:3)
c        print*, zhf2(1:5,:4)
c        print*, zhf2(1:5,:5)
c        print*, zhf2(1:5,:6)
c        print*, idcat,jdcat,icolo,ncstr(iregb),nmode(iel)
c        print*,index(1:6)
        
          do i = 1, ncstr(iregb)
            do j = 1, icolo
            ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j) + signe * zhf2(index(i),j)
            end do
c            if(idcat+i.eq.312) print*, ael(306+6,1:10)
          end do




        end do  ! k

      end do  ! im


      else if(eletyp.eq.'CAR')then
c     ++++++++++++++++++++++++++++
c     Constraints are in rtp, but unknowns are in R,Y,phi !
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
        call zset(5*6, czero, zhf2, 1)
c       zhf2: rows are in r,t,p, columns in R, Y, phi
c       Coefficients cvc(6,:), cvc(7,:) from FOUCOJ are for cos and sin of angle alpha
c       between radial e1 and eR basis vectors.
          if(k .eq. 0)then
c         [Ephi] eq.:
          zhf2(4,5) = cun
c         [Htheta] eq.: Ephi'= Ephi':
          zhf2(5,6) = cun
          end if
c       [Erho]   eq.: Erho =     cos * ER + sin * EY:
        zhf2(1,1) = cvc(k,6)
        zhf2(1,3) = cvc(k,7)
c       [Etheta] eq.: Etheta = - sin * ER + cos * EY:
        zhf2(2,1) = - cvc(k,7)
        zhf2(2,3) =   cvc(k,6)

c       [Hphi] eq.: Etheta'= - sin * ER' - cos * EY':
c NB terms including angle radial derivative are missing
cERN  cvc(k,8) = d/dr(cos THETA) and cvc(k,9) = d/dr(sin THETA) 
c     ADDED in FOUCOJ. Need to update below!!
        zhf2(3,2) = - cvc(k,7)
        zhf2(3,4) = - cvc(k,6)

c@ see
c        zhf(3,1) = cvc(k,2) + ci * mpk * cvc(k,3)
c        zhf(5,1) = cvc(k,1)

          do i = 1, ncstr(iregb)
            do j = 1, icolo
            ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j)
     ;      + signe * zhf2(index(i),j)
            end do
          end do
        end do  ! k

      end do  ! im

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
cERN            call mul2c2(zhf2, zhf, cmat, 5, 5, 6)
		call mul2c2_fast(zhf2(1:5,1:6), zhf, cmat, 5, 5, 6)

              do i = 1, ncstr(iregb)
                do j = 1, icolo
                ael(idcat+i,jdcat+j) = ael(idcat+i,jdcat+j)
     ;          + signe * zhf2(index(i),j)
                end do
              end do
            end do  ! k
          end do  ! im

      end if
C     ++++++
c--------------------------------------------------------------------------
c
c     ccccccccccccccccc
c     c Plasma terms: c
c     ccccccccccccccccc
c  
C     Remains to see: metal walls with hot plasma!:
      if(   .not. vacuum(ireg)
     ;.and. .not.(coljum .and. rbtyp(iregb).eq.'MET')
     ;.and. (coljum .or. (.not.coldpl(ireg) .and. flrops(ireg)))
     ;  )then
 
      trafo = .not. polsym
c      write(603,*),'trafo=',trafo
      ipath = 1
c see with general diel response; now use contrib. as Maxwellian:
      if(.not. cokpco)then
c     ~~~~~~~~~~~~~~~~~~~~
        do mav2 = 2*minf(ielb), 2*msup(ielb)
c	  NB: polfft uses current value of mav2, passed through commod
c                                                 Boundary terms switch:
c                                                            V

cERN: Problems with minf/msup values (they are 0 sometimes!!!!)
cERN: By printing them to screen, the values are OK (WHY ????)
cERN  write(603,*)'intbou:', minf(ielb),msup(ielb) 

        call polfft(v2(0,1), npfft+1, .false., idum, .true., 2, trafo, ipath
     ;  , onlyab_now)

c            write(603,*) mav2, v2(1:10,1)
c            print*,v2(1:10,1),mav2

cPL     ;  , nfft, onlyab)
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
 
          do k = kstop1, kstop2, 2
          mpk = (mav2 + k) / 2
          m = mpk - k
          impk = mpk + 1 - minf(ielb)
          jdcat = jdcato + icolo * (impk - 1)
          kp = iabs(k)
          im = m + 1 - minfc
          idcat = idcato + ncstr(iregb) * (im - 1) 
c                                              Boundary terms switch:
c                                                         V
cPL fixed bug 26/5/04          call plasmb(vl, vls, k, k, v2(0,1), 0, npfft, 2, .false.)
          call plasmb(vl, vls, k, k, v2(0,1), 0, npfft+1, 2, .false.)
cPL         call plasmb(vl, vls, k, k, v2(0,1), 0, npfft, 2, nfft, .false.)
c            write(603,*) mav2, k, vl(1,1)
            if(eletyp.eq.'HEC' .and. rtp)then
c           Constraints expressed in R,T,P coord, between unknowns in +-//:
c           Contribution to cond. [Drho]=0 (radial displacement vector is C0)
            ael(idcat+1,jdcat+1) = ael(idcat+1,jdcat+1) + signe*sqrt2i * vl(1,1)
            ael(idcat+1,jdcat+3) = ael(idcat+1,jdcat+3) + signe*sqrt2i * vl(3,3)
c            else if(eletyp.eq.'M23')then
cc           M23 in plasma to do
c            else
cc           plasma and not rtp to do
            end if

	    end do  ! k
        end do    ! mav2
c        print*, signe, sqrt2i,vl(1,1)
c      print*, ael(306+6,1:20)
c      stop

      else
c     ~~~~
c     Const. k// coordinates
cJAC	call m12cokpco(.false.,idum)
        do mav2 = 2*minf(ielb), 2*msup(ielb)
           if(mav2 .le. minf(ielb) + msup(ielb))then
           kstop2 = min(klim, mav2 - 2 * minf(ielb))
           else
           kstop2 = min(klim, 2 * msup(ielb) - mav2)
           end if
        if(mod(mav2 - kstop2, 2) .ne. 0)kstop2 = kstop2 - 1
        kstop1 = - kstop2
 
         do k = kstop1, kstop2, 2
         mpk = (mav2 + k) / 2
         m = mpk - k
         impk = mpk + 1 - minf(ielb)
         jdcat = jdcato + icolo * (impk - 1)
         kp = iabs(k)
         im = m + 1 - minfc
         idcat = idcato + ncstr(iregb) * (im - 1)
cERN: next lines replace the call to plasmb of usual coordinates         
cE	     vl(1,1) = m12left(mav2-2*minf(ielb)+1,k+klim+1)   ! p = +1
cE	     vl(3,3) = m12right (mav2-2*minf(ielb)+1,k+klim+1) ! p = -1
cE	     vl(5,5) = m12landau(mav2-2*minf(ielb)+1,k+klim+1) ! p = 0 


c
           if(eletyp.eq.'HEC')then
c          Constraints expressed in R,T,P coord, between unknowns in +-//:
c          Contribution to cond. [Drho]=0 (radial displacement vector is C0)
           ael(idcat+1,jdcat+1) = ael(idcat+1,jdcat+1) + signe*sqrt2i * vl(1,1)
           ael(idcat+1,jdcat+3) = ael(idcat+1,jdcat+3) + signe*sqrt2i * vl(3,3)
           else if(eletyp.eq.'M23')then
c          to do
           end if
         end do  ! k
        end do  ! mav2

      end if
c     ~~~~~~
      end if

        end if
c       ======
c--------------------------------------------------------------------------

c     Conjugate block (columns <-> Lagrange multipliers in global matrix):
      do i = jdcato + 1, jdcato + nmode(ielb) * icolo
        do j = idcato + 1, idcato + nmodc * ncstr(iregb)
        ael(i,j) = dconjg( ael(j,i) )
        end do
      end do

c      write(6,*)'Constraints as assembled in ael: ', nmodc*ncstr(iregb)
c	  do i = idcato + 1, idcato + nmodc * ncstr(iregb)
c	  write(nofile,*)i
c       write(6,1000)(ael(i,j),j=jdcato + 1, jdcato + nmode(ielb) * icolo)
c	  end do
c 1000 format(6(1h(, g11.3, 1h,, g11.3, 1h)))

      write(nofile,*)'Exit intbou'
      return
      end
