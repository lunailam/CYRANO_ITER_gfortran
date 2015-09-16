      subroutine outrf2

      implicit none

c     Builds output tables, to be plotted elsewhere.
c     New, shortened version (1998) of OUTPUT, free of graphics software.
 
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
      include 'comphy.copy'
      include 'coequi.copy' ! ERN       
      logical writr
 
      character*3  eletyp
      character*9  charfi(2), charco(6), charpa(3)
      character*27 charto
	character(50) filename	!ERN
      character*1  charfi2(2) !ERN
      character*3  charco2(6) !ERN
      character*2	 charpa2(3) !ERN

      integer 
     ;  i, j, aqp(3), aux
     ;, iplom1, ipol, ii, ii2
     ;, is1, is2 
     ;, ideca
     ;, idllo, icolo, ibulo, icpblo, im
     ;, ie1, nmodc
     ;, ima, iblk, ireq, j1, j2
     ;, istatu, ie2, nintma
     ;, minfc, msupc, nblk
     ;, ipir, ico, ifi, ihi, ipict, itab
     ;
     ;, nsum, i1st, ceilq
     ;, idamin, idmin, idmax
c	 ;, incrho, nrho

      double precision 
     ;  rle
     ;, tabmin, tabmax
     ;, cabs2, rveno
     ;, norfac
     ;, trapeze, gauss
 
      complex*16
     ;  c1, c2
     ;, zdotc, zdotu
     ;, tra1
 
      external 
     ;  nsum, i1st, ceilq, trapeze, gauss, rveno
     ;, zdotc, zdotu, zset, dset, dsum, mucrvz
     ;, idamin, idmin, idmax
 
      double precision 
     ;  second

      external second
c
      data charfi/'Electric ','Magnetic '/
     ;    ,charco/'radial ','poloidal ','toroidal ',
     ;            'L.H.(+) ','R.H.(-) ','parallel '/
     ;    ,charpa/'real','imag','modulus'/
cERN
      data charfi2/'E','H'/
     ;    ,charco2/'rad','pol','tor',
     ;             '_LH','_RH','par'/
     ;    ,charpa2/'re','im','mo'/
 
      cabs2(tra1) = dreal(tra1 * dconjg(tra1))

	 
c**********************************************************************
  
      if(plot2d)then
c     ~~~~~~~~~~~~~~
      write(603,*)'enter outrf2'
      iplom1 = nploth + 1
c     Caution: this must come after parpow has been filled in outpow!

cPL24/01/2005 @ should then add test if(use_outpow)... and put use_outpow in a common

C     Field normalization factor for TRANSP input. Assumes power RFPOW into 1 toroidal mode.
C     sqrt(2) comes from Dirk's convention on E+, E_.
      norfac = 1.d0
c assumes flr?
c      if((plopow.or.transp).and. parpow(nspec+1,2).gt.0.d0)norfac = dsqrt(2.d0 * rfpow / parpow(nspec+1,2))
c      if(transp)write(notrsp,1588)istp(1),iplom1,1
      

cERN	20/12/05
cccccccccccccccccccccccccccccccc Open files for QLFP ccccccccccccccccccccccccccccccccccccccc

c	call date(dia)
c	call time(hora)

c	write(4321,*) 'Saving 2D field profiles for BATCH ...'
c	write(4321,*) 'Folder:', paplda
c	write(4321,*)'Time:', hora

		open (UNIT = 1770, FILE = paplda // 'QLFP_Eplus_2D_re.dat', 
     ;            STATUS = "REPLACE", ACTION = "WRITE")	
c		write(4321,*)'QLFP_Eplus_2D_re.dat (', hora , ')'	
		write(1770,*), 'E+(real): lines=rho(0..rp), cols=theta(0..2pi)'
	      write(1770,*)'Total power = '
              write(1770,*) TOTALPOWER_poynt
		write(1770,*), ' Nrad  ' , ' Npol  ' 
		write(1770,"(i7,i7)"), ipla, iplom1

		open (UNIT = 1771, FILE = paplda // 'QLFP_Eplus_2D_im.dat', 
     ;            STATUS = "REPLACE", ACTION = "WRITE")	
c		write(4321,*)FILE
		write(1771,*), 'E+(imag): lines=rho(0..rp), cols=theta(0..2pi)'
c	      write(1771,*)'Total power = ', POWERFAC*TOTALPOWER_poynt
	      write(1771,*)'Total power = '
              write(1771,*) TOTALPOWER_poynt
		write(1771,*), ' Nrad  ' , ' Npol  '
		write(1771,"(i7,i7)"), ipla, iplom1

		open (UNIT = 1772, FILE = paplda // 'QLFP_Eminus_2D_re.dat', 
     ;            STATUS = "REPLACE", ACTION = "WRITE")	
c		write(4321,*)FILE
		write(1772,*), 'E-(real): lines=rho(0..rp), cols=theta(0..2pi)'
c	      write(1772,*)'Total power = ',  POWERFAC*TOTALPOWER_poynt
	      write(1772,*)'Total power = '
              write(1772,*) TOTALPOWER_poynt
		write(1772,*), ' Nrad  ' , ' Npol  ' 
		write(1772,"(i7,i7)"), ipla, iplom1

		open (UNIT = 1773, FILE = paplda // 'QLFP_Eminus_2D_im.dat', 
     ;            STATUS = "REPLACE", ACTION = "WRITE")	
c		write(4321,*)FILE
		write(1773,*), 'E-(imag): lines=rho(0..rp), cols=theta(0..2pi)'
c	      write(1773,*)'Total power = ',  POWERFAC*TOTALPOWER_poynt
	      write(1773,*)'Total power = '
              write(1773,*) TOTALPOWER_poynt
		write(1773,*), ' Nrad  ' , ' Npol  '
		write(1773,"(i7,i7)"), ipla, iplom1

		open (UNIT = 1774, FILE = paplda // 'QLFP_Epar_2D_re.dat', 
     ;            STATUS = "REPLACE", ACTION = "WRITE")	
c		write(4321,*)FILE
		write(1774,*), 'E//(real): lines=rho(0..rp), cols=theta(0..2pi)'
c	      write(1774,*)'Total power = ',  POWERFAC*TOTALPOWER_poynt
	      write(1774,*)'Total power = '
              write(1774,*) POWERFAC*TOTALPOWER_poynt
		write(1774,*), ' Nrad  ' , ' Npol  '
		write(1774,"(i7,i7)"), ipla, iplom1

		open (UNIT = 1775, FILE = paplda // 'QLFP_Epar_2D_im.dat', 
     ;            STATUS = "REPLACE", ACTION = "WRITE")	
c		write(4321,*)FILE
		write(1775,*), 'E//(imag): lines=rho(0..rp), cols=theta(0..2pi)'
c	      write(1775,*)'Total power = ',  POWERFAC*TOTALPOWER_poynt
	      write(1775,*)'Total power = '
              write(1775,*) TOTALPOWER_poynt
		write(1775,*), ' Nrad  ' , ' Npol  ' 
		write(1775,"(i7,i7)"), ipla, iplom1

c           New: (R,Z)-grid for QLFP

		open (UNIT = 1776, FILE = paplda // 'QLFP_Rmesh.dat', 
     ;            STATUS = "REPLACE", ACTION = "WRITE")	
c		write(4321,*)FILE
		write(1776,*), 'R-mesh: lines=rho(0..rp), cols=theta(0..2pi)'
		write(1776,*), ' Nrad  ', ipla
            write(1776,2222)(abscis(j), j = 1, ipla)   ! rho
		write(1776,*), ' Npol  ', iplom1
            write(1776,2222)(polplo(j), j = 1, iplom1) ! theta
		  do i = 1, ipla	
		     write(1776,2222)(eqt(i, j, 1), j = 1, iplom1)
		  end do
            close(1776)

		open (UNIT = 1777, FILE = paplda // 'QLFP_Zmesh.dat', 
     ;            STATUS = "REPLACE", ACTION = "WRITE")	
c		write(4321,*)FILE
		write(1777,*) 'R-mesh: lines=rho(0..rp), cols=theta(0..2pi)'
		write(1777,*) ' Nrad  ', ipla
            write(1777,2222)(abscis(j), j = 1, ipla)   ! rho
		write(1777,*) ' Npol  ', iplom1
            write(1777,2222)(polplo(j), j = 1, iplom1) ! theta
		  do i = 1, ipla	
		     write(1777,2222)(eqt(i, j, 2), j = 1, iplom1)
		  end do
           close(1777)

cccccccccccccccccccccccccccccccccccccccccccccc


            print*, 'POWERFAC = ', POWERFAC
            print*, 'Total RF power =',  POWERFAC*TOTALPOWER_poynt


c      include 'ope2df.copy'
 
c     plots in poloidal plane required by ploty:
c     -------------------------------------------
      write(603,*)'enter do 18'
      do 18 ihi = 1, 3 * nficom
        if(ploty(1,ihi).or.ploty(2,ihi).or.ploty(3,ihi) .or. (transp .and. (ihi.eq.7 .or. ihi.eq.8)))then
        call dset( (maxthp+1)*nabplo*3, 0.d0, tab(1,1,1), 1)
          if(ihi.gt.ndof)then
c         1,3,5:
          ii = 2 * (ihi - ndof) - 1
          else if(ihi.gt.nficom)then
c         1,2,3:
          ii = ihi - nficom
          else
c         1,2,3:
          ii = ihi
          end if
 
          if(ihi.ge.4 .and. ihi.le.6)then
          ifi = 2
          else
          ifi = 1
          end if
          if(ihi.le.3)then
          ico = ihi
          else
          ico = ihi - 3
          end if
 
        iblk = 1
        nblk = ceilq(nwn*nabplo*nmode(imax), iobll)

          do 309 imoto = 1, ntotor
c     Read modal components:
      if(.not.monoto)then
      nblk = ceilq(nwn*5*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqread(aqp,xrtp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*3*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqread(aqp,hrtp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*ndof*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqread(aqp,xpmp,iblk,nblk,ireq,1,istatu)
      iblk = iblk + nblk
      nblk = ceilq(nwn*nficom*nabplo*nmode(imax),iobll)
      ireq = ireq + 1
      call aqread(aqp,xpmp2,iblk,nblk,ireq,0,istatu)
      iblk = iblk + nblk
      call aqwait(aqp,istatu)
      end if
 
      i = 0
      do 17 ireg = 1, nreg
c     ===================      writr = transp .and. ireg.eq.1
      is1 = nsum(ireg-1, ns) + 1
      is2 = is1 + ns(ireg) - 1
      do 17 isubr = is1, is2
c     =====================      eletyp = styp(isubr)
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo
      ie1 = i1st(isubr) + 1
      ie2 = ie1 + iele(isubr) - 1
 
      do 17 iel = ie1, ie2
c     --------------------
c     write(6,*)'in do 17 iel= ',iel
 
      rle = fl(iel)
      if(outgau)then
      j1 = 1
      if(iel.eq.ifiel(ireg))j1 = 0
      j2 = ngauss + 1
c      j2 = ngauss
c      if(iel.eq.ilael(ireg))j2 = ngauss+1
      else
      nintma = ninter
      if (iel .eq. ilael(ireg)) nintma = ninter + 1
      j1 = 0
      j2 = nintma
      end if
      ielm = iel
      if(nmode(iel-1) .gt. nmode(iel))ielm = iel - 1
      nmodc = nmode(ielm)
      minfc = minf(ielm)
      msupc = msup(ielm)
      
                 do 17 j = j1, j2
c                ----------------
      i = i + 1
 
      call zset((ndof-1)*nploth, czero, eloc, 1)
      call zset(ndof*nploth, czero, epmloc, 1)
      call zset((nficom)*nploth, czero, hloc, 1)
      
        do im = 1, nmodc
        m = im - 1 + minfc
        ima = m + 1 - minf(imax)
        call zcopy(iplom1, tabexp(ima,1), maxpom, copisi, 1)
 
          if(ihi.le.ndof)then
          ii2 = 2 ** (ii - 1)
            if(ihi .le. nficom)then
            c1 = xrtp(ii2,i,im)
            call zaxpy(nploth, c1, copisi, 1, eloc(ii2,1), ndof-1)
            else
            c2 = hrtp(ii,i,im)
            call zaxpy(nploth, c2, copisi, 1, hloc(ii,1), nficom)
            end if
          else
          c1 = xpmp(ii,i,im)
          call zaxpy(nploth, c1, copisi, 1, epmloc(ii,1), ndof)
          end if
c following: eloc and hloc must be defined!!! eloc def has changed.
cc    Compute local Poynting vector:
c     call ccrosp(poyloc, eloc(1,ipol), hloc(1,ipol))
c     poyloc(1) = 0.5*poyloc(1)
c     poyloc(2) = 0.5*poyloc(2)
c     poyloc(3) = 0.5*poyloc(3)
c     if(.not.glovac)then
cc    Poynting vector real part poloidal projection:
c     e1 = dreal(poyloc(1))
c     e2 = dreal(poyloc(2))
c     else
cc    Poynting vector imag.part poloidal projection:
c     e1 = dimag(poyloc(1))
c     e2 = dimag(poyloc(2))
c     end if
c 203 continue
c 201 continue
c     end if 
        end do
 
        if(ihi.le.ndof)then
          if(ihi.le.nficom)then
	    if(cokpco)eloc(ii2,1:nploth) = eloc(ii2,1:nploth) * einphb(i,1:nploth)
            do ipol = 1, nploth
            tab(ipol,i,1) = tab(ipol,i,1) + dreal(eloc(ii2,ipol))
            tab(ipol,i,2) = tab(ipol,i,2) + dimag(eloc(ii2,ipol))
            tab(ipol,i,3) = tab(ipol,i,3) + cdabs(eloc(ii2,ipol))
            end do
          else
	    if(cokpco)hloc(ii,1:nploth) = hloc(ii,1:nploth) * einphb(i,1:nploth)
            do ipol = 1, nploth
            tab(ipol,i,1) = tab(ipol,i,1) + dreal(hloc(ii,ipol))
            tab(ipol,i,2) = tab(ipol,i,2) + dimag(hloc(ii,ipol))
            tab(ipol,i,3) = tab(ipol,i,3) + cdabs(hloc(ii,ipol))
            end do
          end if
        else
	    if(cokpco)epmloc(ii,1:nploth) = epmloc(ii,1:nploth) * einphb(i,1:nploth)
          do ipol = 1, nploth
          tab(ipol,i,1) = tab(ipol,i,1) + dreal(epmloc(ii,ipol))
          tab(ipol,i,2) = tab(ipol,i,2) + dimag(epmloc(ii,ipol))
          tab(ipol,i,3) = tab(ipol,i,3) + cdabs(epmloc(ii,ipol))
          end do
        end if
  17  continue

 309  continue
 
        do itab = 1, 3
c       For polar coords.: repeat data of first angle
        call dcopy(ist11, tab(1,1,itab), maxthp+1, tab(iplom1,1,itab), maxthp+1)
        end do
 

cERN  13/01/06 cccccccccccccccccc QLFP fields OUTPUT ccccccccccccccccccccccccccc

c           2D fields: rho = 1..Npla (lines), theta = 1..2pi (cols)
c           ---------
c          

c           (Factor 1/sqrt(2) still MISSING !!!!!!! -> DONE by DIRK!!!!!)

c           1) E+
               if(ihi.eq.7)then
	            do i = 1, ipla
	               write(1770,2222)(tab(j,i,1), j = 1, iplom1)
                  end do
	            print *, 'Re E+ for QLFP ... OK!'
	            do i = 1, ipla
	               write(1771,2222)(tab(j,i,2), j = 1, iplom1)
                  end do
	            print *, 'Im E+ for QLFP ... OK!'
	         end if

c           2) E-
               if(ihi.eq.8)then
	            do i = 1, ipla
	               write(1772,2222)(tab(j,i,1), j = 1, iplom1)
                  end do
	            print *, 'Re E- for QLFP ... OK!'
	            do i = 1, ipla
	               write(1773,2222)(tab(j,i,2), j = 1, iplom1)
                  end do
	            print *, 'Im E- for QLFP ... OK!'
	         end if

c           3) E//
               if(ihi.eq.9)then
	            do i = 1, ipla
	               write(1774,2222)(tab(j,i,1), j = 1, iplom1)
                  end do
	            print *, 'Re E// for QLFP ... OK!'
	            do i = 1, ipla
	               write(1775,2222)(tab(j,i,2), j = 1, iplom1)
                  end do
	            print *, 'Im E// for QLFP ... OK!'
	         end if

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cERN       Commented on Friday 20 Oct 2006 (before JET)

c          Uncommented on 26/03/2008
   
cERN	ccccc Writing 2D fields to separate files cccccccccccccc
        do ipict = 1, 3
          if(ploty(ipict,ihi))then
          ipir = ipict
          ideca = (maxthp+1) * (nabplo*(ipir-1))
          tabmin = tab2(ideca + idmin((maxthp+1)*ist11, tab(1,1,ipir), 1))
          tabmax = tab2(ideca + idmax((maxthp+1)*ist11, tab(1,1,ipir), 1))
            if(tabmin.ne.tabmax  .and. (abs(tabmin).gt.1.d-60 .or. abs(tabmax).gt.1.d-60))then
		iglplo = iglplo + 1
		charto = charfi(ifi) // charco(ico) // charpa(ipict)
		filename = charfi2(ifi) // charco2(ico) // "_2D_"// charpa2(ipict) // ".dat"
		print *, filename
			 
		write(nofile,1000) iglplo, '2d ' // charto
		ploleg(iglplo) = '2d ' // charto
		open (UNIT = 777, FILE = paplda // filename, STATUS = "REPLACE", ACTION = "WRITE")	
		write(777,*), charto
		write(777,*), ' rad ' , ' pol ', ' pla ' 
		write(777,"(i5,i5,i5)"), ist11, iplom1, ipla
	        do i = 1, ist11
	        aux = intaplot(i)
              write(777,2222)(tab(j,aux,ipir), j = 1, iplom1)
              end do
	      close(777)
            end if ! (tabmin.ne.tabmax)
          end if ! (ploty(ipict,ihi))
        end do ! (ipict = 1, 3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        end if ! (?????)
 
  18  continue  ! (ihi = 1, 3*nficom)

      end if ! (plot2d)
c     ~~~~~~

	close(1770)
	close(1771)
	close(1772)
	close(1773)
	close(1774)
	close(1775)

      if(.not.monoto)call aqclose(aqp,istatu)
      
      return
 1000 format(1h ,'Plot #',i4,'  ',a50)
 1488 format((' ',4(1x,1p,d11.3)))    
 1588 format(' ',3(1x,i4))          
 2222 format(1000(g14.6))  
      end

