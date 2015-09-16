      subroutine ascurl

      implicit none

c     Assembly of curl.curl and displacement current terms for current element.

c     No account of symmetry in curl.term
c     07/09/89: worked out further optimization in small matrix products

      include 'pardim.copy'

c     Special use of COMSOL: redefined larger!!
c     Assumes v, totop not used in that common elsewhere!!
c     NB: continuation line is not allowed after cf77 include '...' !
c     Careful here:
c     File COMUSR2 must be updated every time file COMUSR changes!
c      COMPLEX*16 V(6,6,2*MAXPOM**2), TOTOC(12,12,BLLEN)
      include 'comusr2.f'
c     ;, v, totop

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
      include 'comphy.copy'

      character*3 eletyp

      logical sing, sing1, inl, inr, jnl, jnr

      integer 
     ;  in(6), jn(6), index(12), jndex(12)
     ;, ipop1(3), jpop1(3), ipop2(3), jpop2(3), npopi, npopj
     ;, i, j
     ;, nmodc, indext, m1, m2, ni, nj, k1, k2, ip, jp, kp
     ;, ideca, jdeca, id, jd, nset0
     ;, ivab, ivabs, ivabs0, i1, i2, j1, j2, jndext
     ;, idllo, icolo, ibulo, icpblo
     ;, ielet
     ;, ns0s, nsingb, nnsib, insosb(maxpom**2), ib
     ;, nsum, i1st

      integer
     ;  ir, jr, lr, lc, itr, itc
      
      double precision 
     ;  cabs1, rle
     ;, second
     ;, rea

      complex*16 zdum

      complex*16 
     ;  axl(3), axc(3), tem3, t, fac
     ;, rotax(6,6,maxpom*maxpom), rotaxs(6,6,maxpom*maxpom)
c     ;, to(12,12), toto(12*12)
c     ;, tos(12,12), totos(12*12)
c     ;, big(12,12), bigs(12,12), big1(12*12), big1s(12*12)

c      equivalence (toto,to)
c     ;, (totos,tos), (big1,big), (big1s,bigs)

      external second, nsum, i1st

      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
  
c     FAC: common factor to all matrix elements
C     Remember field spectrum has dimensions:
C       - V/m in toroidal geometry
C       - V   in the cylinder limit
      fac = glofac * twopi ** 2
      if(.not.cyl)fac = fac * ra

      eletyp = styp(isubr)
      ielet = istyp(isubr)
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo

      rle = fl(iel)
      nmodc = max0(nmode(iel-1), nmode(iel))
 
C     Index list for assembly:
      indext = nmodc * icpblo
        do  i = 1, icolo
        index(i) = i
        jndex(i) = i
        index(i+icpblo) = i + indext
        jndex(i+icpblo) = i + indext
        end do
 
      if(iel.eq.1)then
c     ================
c     first element at magnetic axis or crown inner boudary
c     !nmode(0) = nmode(1) assumed!

      sing1 = .not. crown 

c       nnsib: number of nonsing. blocks
        if(ncrot .ge. nmodc-1)then
        nnsib = nmodc ** 2
        else
        nnsib = (ncrot + 1) * (2 * nmodc - ncrot) - nmodc
        end if
        nset0 = nnsib
c     ns0s = nmodc
        if(circ)then
c       sing terms occur only when k=0 (poloidal mode with itself)
        ns0s = nmodc
        else
        ns0s = nset0
        end if
      ivabs0 = nset0

      if(sing1)then
c     nset0 to include singular blocks:
      nset0 = nset0 + ns0s 
      end if

      if(nset0.gt.bllen)write(nofile,*)'error ascurl elt.1: nset0=',nset0
     ;, ', but bllen=',bllen
     
cPL30/5/04      call zset(144*nset0, czero, totop, 1)
      call zset(144*nset0, czero, totoc, 1)
      m1 = minf(iel)
      m2 = msup(iel)

c     For first element, compute at magnetic axis too with ig=0.
c     This is used in treatment of singular terms.
      do ig = 0, ngauss
c     ~~~~~~~~~~~~~~~~~
      intab = ig + 1
      y = abscno(intab)
      yinv = abscni(intab)

        if(circ)then
        call fouco
        else
        call foucog(1)
        end if

      ivab = 0
      nsingb = 0
      ivabs = ivabs0

      do mr = 1, nmodc
c     ----------------
c     Mode indices: here m for test function, mpk=m+k for electric field.
c     With this convention, the + sign is required in poloidal FFT exponent (FOUCO and FOUCOG).
      m = mr - 1 + m1
      k1 = max0(-ncrot, m1-m)
      k2 = min0(ncrot, m2-m)
 
        do k = k1, k2
        kp = iabs(k)
        mpk = m + k
        ivab = ivab + 1
          if(circ)then
          sing = sing1 .and. k.eq.0
          else
          sing = sing1
          end if
          if(sing)then
c         increment number of singular blocks:
          nsingb = nsingb + 1
          ivabs = ivabs + 1
c         store index of associated nonsingular block:
          insosb(nsingb) = ivab
          end if
        
        if(ig.eq.0)then
        if(sing)call vacmbg(rotax(1,1,nsingb), rotaxs(1,1,nsingb), sing)
        else
        call vacmbg(vmat(1,1,ivab), vmat(1,1,ivabs), sing)
c  2003 what did this mean:???
c       In circular concentric case, singular terms are indep. of y,
c       hence next block does not contribute:
          if(sing .and. .not.circ)then
cERN            do j = 1, 6
cERN              do i = 1, 6
cERN              t = vmat(i,j,ivabs) - rotaxs(i,j,nsingb)
cERN              vmat(i,j,ivab) = vmat(i,j,ivab) + t * yinv
cERN              end do
cERN            end do
              vmat(1:6,1:6,ivab) = vmat(1:6,1:6,ivab) + (vmat(1:6,1:6,ivabs) - rotaxs(1:6,1:6,nsingb)) * yinv
          end if
        end if
        
        end do

      end do
c     ------
 
c     Nonsingular terms:
c     Left and right multipliy vmat by element basis functions, store in totop
        if(ig.ne.0)then
        tem3 = rle * wga(ig)
        call mu3tf2bl_fast(totoc, tem3, ldif(1,1,ig), vmat, ivab, idllo)
cPL30/5/04        call mu3tf2bl(totop, tem3, ldif(1,1,ig), vmat, ivab, idllo)
        end if

      end do
C     ~~~~~~

c      write(nofile,*)' '
c      write(nofile,*)'genera, element 1: raw totop block #1:'
c      write(nofile,3000)((totop(i,j,1),j=1,idllo),i=1,idllo)
c 3000 format(1h ,12(1x,g11.3))

c To install in the future, replaces LDIF multiplication...
C NEW: WOULD ALLOW DISCARDING LDIF MATRIX IF DONE FOR OTHER ELTS. AND PLASMA:
C IS THIS FASTER?
c change big...
c
c      call zset(12*12, czero, big1, 1)
c      do ir = 1, nficom
c      itr = ibafty(ir,ielet)
c      do i1 = 0, 1
c      i2 = 2*ir-1 + i1
c
c       do jr = 1, nficom
c       itc = ibafty(jr,ielet)
c       do j1 = 0, 1
c       j2 = 2*jr-1 + j1
cc          t = totop(i2,j2,)
c
cc              if(t.ne.czero)then
c                 do lr = 1, nbafel(ir,ielet)
c                 i = jbfloc(lr,ir)
c                   do lc = 1, nbafel(jr,ielet)
c                   j = jbfloc(lc,jr)
c
c                   do ib = 1, @
c                   big(i,j) = big(i,j) + totop(i2,j2,@)
c     ;              * bafgn(lr,i1,itr,ig) * bafgn(lc,j1,itc,ig)
c     ;              * rle**(ideriv(lr,itr)-i1+ideriv(lc,itc)-j1)
c                   end do
c
c                   end do
c               end do
cc              end if
c
c       end do
c       end do
c
c      end do
c      end do
c      
c      call zaxpy(12*12, tem3, big1, 1, toto, 1)

c     Now deal with singular terms:

      if(sing1)then
c     -------------
c     Loop over field components, row-wise:
      do ir = 1, nficom
      itr = ibafty(ir,ielet)
c     Field and its derivative:
      do i1 = 0, 1
      i2 = 2 * ir - 1 + i1

c       Loop over field components, column-wise:
        do jr = 1, nficom
        itc = ibafty(jr,ielet)
c         Field and its derivative:
          do j1 = 0, 1
          j2 = 2 * jr - 1 + j1
        
c           Loop over basis functions, row-wise:
            do lr = 1, nbafel(ir,ielet)
            i = jbfloc(lr,ir)
c             Loop over basis functions, column-wise:
              do lc = 1, nbafel(jr,ielet)
              j = jbfloc(lc,jr)

                if(i1.eq.0 .and. j1.eq.0 .and. lr.eq.1 .and. lc.eq.1)then
c               The test selects the only terms with singularities:
c               non-derivative dofs of the first basis functions.
c               (space was already set to zero in totop singular blocks)
C               The sings. are removed by the boundary conds; 
C               we must only retain terms with a nonvanishing
C               power of the reduced coordinate (i.e. prior to division by y).
C
C     Warning if code transformed to solve another equation:
C     Here, an assumption: derivative of first basis function must vanish
C     at left element node; if not, corresp. element of ROTAXS must vanish.
C     This is the case for curl.curl operator provided Lagrange basis functions
C     (type 2) are only used with component Erho. No problem with standard
C     element types 'HEC'(#1) and 'M23'(#2).
C     ROTAXS only has nonzero entries in row/col 1, 3, 5.
C     @May 2004 New element type 'CAR'(#3): check!
C
                  do ig = 1, ngauss
c                 yinv = abscni(ig+1)
c                 rea = bafgn(lr,i1,itr,ig) * bafgn(lc,j1,itc,ig)
c     ;                - bafgn(lr,i1,itr,0) * bafgn(lc,j1,itc,0)
cc or:
c       rea =
c     ;(bafgn(lr,i1,itr,ig) - bafgn(lr,i1,itr,0)) * bafgn(lc,j1,itc,ig) +
c     ; bafgn(lr,i1,itr,0) * (bafgn(lc,j1,itc,ig) - bafgn(lc,j1,itc,0))
c       rea = rea * wga(ig) * yinv *
c     ;                  rle**(1+ideriv(lr,itr)-i1+ideriv(lc,itc)-j1)
cC Power of RLE: 1 for integration, ideriv for each 'derivative' basis func.,
cC -i1, -j1 for field derivatives.

c     To avoid division by y:
                  rea = bafsin(lr,i1,itr,ig) * bafgn(lc,j1,itc,ig) +
     ;                  bafgn(lr,i1,itr,0) * bafsin(lc,j1,itc,ig)
                  rea = rea * wga(ig) *
     ;                  rle**(ideriv(lr,itr)-i1+ideriv(lc,itc)-j1)
c Net power of RLE: 1 for integration, ideriv for each 'derivative' basis func.,
c -i1, -j1 for field derivatives, -1 for normalis. of BAFSIN.

                      do ib = 1, nsingb
cPL30/5/04                      totop(i,j,insosb(ib)) = totop(i,j,insosb(ib))
                      totoc(i,j,insosb(ib)) = totoc(i,j,insosb(ib)) + rea * rotaxs(i2,j2,ib)
                      end do
                  end do

                else
C               add other terms normally, in TO:
                  do ig = 1, ngauss
                  yinv = abscni(ig+1)
c                 rea = wga(ig) * yinv
c     ;                 * bafgn(lr,i1,itr,ig) * bafgn(lc,j1,itc,ig)
c     ;                  * rle**(1+ideriv(lr,itr)-i1+ideriv(lc,itc)-j1)
cC Power of RLE: 1 for integration, ideriv for each 'derivative' basis func.,
cC -i1, -j1 for field derivatives.
c     To avoid division by y:
                    if(leadeg(lr,itr)-i1.gt.0)then
                    rea = bafsin(lr,i1,itr,ig) * bafgn(lc,j1,itc,ig)
                    else if(leadeg(lc,itc)-j1.gt.0)then
                    rea = bafgn(lr,i1,itr,ig) * bafsin(lc,j1,itc,ig)
                    else
                    rea = bafgn(lr,i1,itr,ig) * bafgn(lc,j1,itc,ig) * yinv * rle
                    end if
                  rea = rea * wga(ig)
     ;                  * rle**(ideriv(lr,itr)-i1+ideriv(lc,itc)-j1)
C Power of RLE: 1 for integration, ideriv for each 'derivative' basis func.,
C -i1, -j1 for field derivatives, -1 for normalis. of BAFSIN.

                      do ib = 1, nsingb
cPL30/5/04                      totop(i,j,insosb(ib)) = totop(i,j,insosb(ib))
                      totoc(i,j,insosb(ib)) = totoc(i,j,insosb(ib)) + rea * rotaxs(i2,j2,ib)
                      end do
                  end do

                end if

            end do
          end do
        end do
        end do

      end do
      end do

      end if  ! sing1 case
c     ------ 

cPL30/5/04      call zscal(12*12*nset0, fac, totop, 1)
      call zscal(12*12*nset0, fac, totoc, 1)
c      write(nofile,*)' '
c      write(nofile,*)'genera, element 1: totop bl.#1 after treatment of sing:'
c      write(nofile,3000)((totop(i,j,1),j=1,idllo),i=1,idllo)

c# progressively replaced by above:
      indext = nmodc * icpblo
        do i = 1, icolo
        index(i+icpblo) = i + indext
        jndex(i+icpblo) = i + indext
        end do
      m1 = minf(iel)
      m2 = msup(iel)
      ivab = 0
      ivabs = ivabs0
      
                do 12 mr = 1, nmode(iel)
c               ------------------------
      m = mr - 1 + m1
      if(.not.bcinb1)call bouas1('R', ni, npopi, axl, ipop1, ipop2, in)

      if(eletyp.eq.'M23')
     ;index(icolo+1) = (m - m1) * ibulo + 1 + (msup(iel-1) - m + 1) * icolo

      k1 = max0(-ncrot, m1-m)
      k2 = min0(ncrot, m2-m)
      ideca = (mr - 1) * icolo
 
      do 120 k = k1, k2
c     -----------------
      kp = iabs(k)
      mpk = m + k
      ivab = ivab + 1
      jdeca = ideca + k * icolo

        if(circ)then
        sing = sing1 .and. k.eq.0
        else
        sing = sing1
        end if
      if(sing)ivabs = ivabs + 1
 
      if(eletyp.eq.'M23')jndex(icolo+1) = (mpk - m1) * ibulo + 1 + (msup(iel-1) - mpk + 1) * icolo
 
      if(.not.bcinb1)then
c     +++++++++++++++++++
c     Traditional option for axis boundary conditions:
      call bouas1('C', nj, npopj, axc, jpop1, jpop2, jn)
        DO IP = 1, NPOPI
          do j = 1, idllo
cPL30/5/04   totop(ipop1(ip),j,ivab) = totop(ipop1(ip),j,ivab) + totop(ipop2(ip),j,ivab) * axl(ip)
          totoc(ipop1(ip),j,ivab) = totoc(ipop1(ip),j,ivab) + totoc(ipop2(ip),j,ivab) * axl(ip)
          end do
        end do
        do jp = 1, npopj
          do i = 1, idllo
cPL30/5/04   totop(i,jpop1(jp),ivab) = totop(i,jpop1(jp),ivab) + totop(i,jpop2(jp),ivab) * axc(jp)
          totoc(i,jpop1(jp),ivab) = totoc(i,jpop1(jp),ivab) + totoc(i,jpop2(jp),ivab) * axc(jp)
          end do
        end do
 
      do j = 1, nj
      jd = jdeca + jndex(jn(j))
        do i = 1, ni
        id = ideca + index(in(i))
cPL30/5/04        ael(id,jd) = ael(id,jd) + totop(in(i),jn(j),ivab)
        ael(id,jd) = ael(id,jd) + totoc(in(i),jn(j),ivab)
        end do
        do i = icpblo + 1, idllo
        id = ideca + index(i)
cPL30/5/04        ael(id,jd) = ael(id,jd) + totop(i,jn(j),ivab)
        ael(id,jd) = ael(id,jd) + totoc(i,jn(j),ivab)
        end do
      end do
 
      do j = icpblo+1, idllo
      jd = jdeca + jndex(j)
        do i = 1, ni
        id = ideca + index(in(i))
cPL30/5/04        ael(id,jd) = ael(id,jd) + totop(in(i),j,ivab)
        ael(id,jd) = ael(id,jd) + totoc(in(i),j,ivab)
        end do
        do i = icpblo+1, idllo
        id = ideca + index(i)
cPL30/5/04        ael(id,jd) = ael(id,jd) + totop(i,j,ivab)
        ael(id,jd) = ael(id,jd) + totoc(i,j,ivab)
        end do
      end do

      else
c     ++++
c     New option at axis for general geometry; essential (+ nat. if required)
c     conds. enforced by call to AXISBC in GENERA2.
c     
      do j = 1, idllo
      jd = jdeca + jndex(j)
        do i = 1, idllo
        id = ideca + index(i)
cPL30/5/04        ael(id,jd) = ael(id,jd) + totop(i,j,ivab)
        ael(id,jd) = ael(id,jd) + totoc(i,j,ivab)
        end do
      end do
        
      end if
c     ++++++

  120 continue
   12 continue
c     --------

c      write(nofile,*)' '
c      write(nofile,*)'genera, element 1: totop bl.#1 after BOUAS:'
c      write(nofile,3000)((totop(i,j,1),j=1,idllo),i=1,idllo)
 
c     End of code for first element.

      else
c     ====
c     Other elements
C
c     Now decided no special treatment of 1/y terms, so that
      sing = .false.
 
        if(nmode(iel-1).gt.nmode(iel))then
        ielm = iel - 1
        else
        ielm = iel
        end if

      nmodc = nmode(ielm)
      if(nmodc.gt.0)then
c     ================== 
        if(ncrot.ge.nmodc-1)then
        nset0 = nmodc ** 2
        else
        nset0 = (ncrot + 1) * (2 * nmodc - ncrot) - nmodc
        end if

      if(circ)then
c     sing terms could occur only when k=0
      ns0s = nmodc
      else
      ns0s = nset0
      end if
      ivabs0 = nset0
c     nset0 includes singular blocks:
      if(sing)nset0 = nset0 + ns0s

cPL30/5/04      call zset(144*nset0, czero, totop, 1)
      call zset(144*nset0, czero, totoc, 1)
      m1 = minf(ielm)
      m2 = msup(ielm)
 
      do 703 ig = 1, ngauss
c     ---------------------
      intab = ig + (iel - 1) * (ngauss + 1) + ireg
      y = abscno(intab)
      yinv = abscni(intab)

        if(circ)then
        call fouco
        else
        call foucog(1)
        end if

      ivab = 0
      ivabs = ivabs0

      do mr = 1, nmodc
c     ----------------
c     Mode indices: here m for test function, mpk=m+k for electric field.
c     With this convention, the + sign is required in poloidal FFT exponent.
      m = mr - 1 + m1
      k1 = max0(-ncrot, m1-m)
      k2 = min0(ncrot, m2-m)
 
        do k = k1, k2
        kp = iabs(k)
        mpk = m + k
        ivab = ivab + 1
c       sing = k.eq.0
c       if(sing)ivabs = ivabs + 1
        call vacmbg(vmat(1,1,ivab), vmat(1,1,ivabs), sing)
 
c 	  if(abscis(intab)>0.1) then
c		 call OUT_VMAT(m,k,ivab)
c	  end if
 
        end do

      end do
c     ------
c	call close_VMAT

cERN	ccccccccccccccc OUTPUT for DEBUGGING ccccccccccccccccccccccc
c	( should be commented out in usual runs)


cERN	  DEBUGGING
c	  if(abscis(intab)>0.4) then
c	     write(11111,"(f15.5)"), dreal(vmat(1,1,:))
c	     write(22222,"(f15.5)"), dreal(vmat(3,3,:))
c	     write(33333,"(f15.5)"), dreal(vmat(5,5,:))
c	     stop
c	  end if


c	  if(abscis(intab)>0.1) then
c	     stop
c	  end if


c	open (UNIT = 777, FILE = 'vmat_stand_re.dat', 
c     ;      STATUS = "REPLACE", ACTION = "WRITE")
c		do j = 1, 1000	
c		   write(777,'(6f16.8)')(dreal(vmat(k,1,j)), k = 1, 6)
c		end do
c	close(777)
c	open (UNIT = 777, FILE = 'vmat_stand_im.dat', 
c    ;      STATUS = "REPLACE", ACTION = "WRITE")
c		do j = 1, 1000	
c		   write(777,'(6f16.8)')(dimag(vmat(k,1,j)), k = 1, 6)
c		end do
c	close(777)
c	Results are equivalent for dshape = FALSE
c      print *,  vmat(1,1:6,1)
c	print *
c	print *, 'Stop at ASCURL'
c	stop
cERN  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      tem3 = rle * wga(ig) * fac
cPL30/5/04      call mu3tf2bl(totop, tem3, ldif(1,1,ig), vmat, ivab, idllo)
      call mu3tf2bl_fast(totoc, tem3, ldif(1,1,ig), vmat, ivab, idllo)
cERN      call mu3tf2bl(totoc, tem3, ldif(1,1,ig), vmat, ivabs, idllo)
 703  continue
c     --------

c      ivabs1 = nset0 - ns0s
c      ivabs = ivabs1
      ivab = 0
      ivabs = ivabs0

      sing = .false.
      if(sing)then
      do m = m1, m2
      k1 = max0(-ncrot, m1-m)
      k2 = min0(ncrot, m2-m)
        do k = k1, k2
        ivab = ivab + 1
        kp = iabs(k)
c       sing = k.eq.0
          if(sing)then
          ivabs = ivabs + 1
	    totoc(1:12,1:12,ivab) = totoc(1:12,1:12,ivab) + totoc(1:12,1:12,ivabs) * yinv
          end if
        end do
      end do
      end if
 
      ivab = 0
                do 222 mr = 1, nmodc
c               --------------------
      m = mr - 1 + m1
c     index list for assembly:
      inl = nmode(iel-1).gt.0 .and. m.ge.minf(iel-1) .and.
     ;      m.le.msup(iel-1)
      inr = nmode(iel).gt.0 .and. m.ge.minf(iel) .and.
     ;      m.le.msup(iel)
      indext = nmodc*ibulo
      if(inl) indext = indext + (msup(iel-1)-m+1) * icolo
      if(inr) indext = indext + (m-minf(iel)) * icolo

        do i = 1, icolo
        index(i+icpblo) = i + indext
        end do
        
        if(eletyp.eq.'M23')then
        index(icolo+1) = (m - m1) * ibulo + 1
        if(inl)index(icolo+1) = index(icolo+1) + (msup(iel-1) - m + 1) * icolo
        end if

        if(inl)then
        i1 = 1
        else
        i1 = icolo + 1
        end if
        if(inr)then
        i2 = idllo
        else
        i2 = icpblo
        end if
      k1 = max0(-ncrot, m1-m)
      k2 = min0(ncrot, m2-m)
 
      do 222 k = k1, k2
c     -----------------
      kp = iabs(k)
      mpk = m + k
      ivab = ivab + 1

c     index list for assembly:
      jnl = nmode(iel-1).gt.0 .and. mpk.ge.minf(iel-1) .and. mpk.le.msup(iel-1)
      jnr = nmode(iel).gt.0 .and. mpk.ge.minf(iel) .and. mpk.le.msup(iel)
      jndext = nmodc*ibulo
      if(jnl)jndext = jndext + (msup(iel-1)-mpk+1)*icolo
      if(jnr)jndext = jndext + (mpk-minf(iel))*icolo

        do i = 1, icolo
        jndex(i+icpblo) = i + jndext
        end do
      
        if(eletyp.eq.'M23')then
        jndex(icolo+1) = (mpk - m1) * ibulo + 1
        if(jnl)jndex(icolo+1) = jndex(icolo+1) + (msup(iel-1) - mpk + 1) * icolo
        end if
 
        if(inl)then
        ideca = (m - minf(iel-1)) * icolo
        else
        ideca = nmode(iel-1) * icolo
        end if

        if(jnl)then
        j1 = 1
        jdeca = (mpk - minf(iel-1)) * icolo
        else
        j1 = icolo + 1
        jdeca = nmode(iel-1) * icolo
        end if

        if(jnr)then
        j2 = idllo
        else
        j2 = icpblo
        end if
 
        do j = j1, j2
        jd = jdeca + jndex(j)
          do i = i1, i2
          id = ideca + index(i)
cPL30/5/04          ael(id,jd) = ael(id,jd) + totop(i,j,ivab)
          ael(id,jd) = ael(id,jd) + totoc(i,j,ivab)
cERN	-----------------------------------------------
c 	  if(abscis(intab)>0.1) then
c	     write(1111,"(i10, i10, f15.6, f15.6)"), m, k, 
c     ;                  dreal(ael(id,jd)), dimag(ael(id,jd))
cc		print *, m, k
c	  end if
c	-----------------------------------------------
          end do
        end do
 
 222  continue
c     --------
cERN	----------------------------------------------- 
c	  if(abscis(intab)>0.1) then
c	     stop
c	  end if
c	-----------------------------------------------


      end if
c     ======
 
 
      end if
c     ======
 
 
      return

      end
