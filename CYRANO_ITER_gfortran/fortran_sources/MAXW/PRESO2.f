      subroutine presol(FLAG, OPT)

      implicit none

c     This routine prepares some index arrays for Cyrano runs.
 
      include 'pardim.copy'

      include 'comgeo.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comber.copy'
      include 'comant.copy'
      include 'comfin.copy'
      include 'comfic.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comin2.copy'
      include 'commod.copy'
      include 'comequ.copy'
      include 'compla.copy'
      include 'comphy.copy'
      include 'complp.copy'

	 logical, intent(in) :: FLAG, OPT		!ERN	
         !! integer, intent(in) :: FLAG, OPT
c	FLAG = .false. : normal use, do not re-compute number of pol. modes
c	FLAG = .true.  : re-compute number of pol. modes according to surface perimeter
c					 and the value pollength given in Cyrano.f 	
      character*3 eletyp

      integer 
     ;  ireg, is1, is2, ifirse, ilaste, isubr, idllo, icolo, ibulo
     ;, icpblo, iel, nmomax, nmodc, nmomlo, imaxlo, ich
     ;, ie1, ie2, nonec, ithoma, ifim1, totvec, totbub
     ;, i1st, nsum, aux, msup_aux, minf_aux
	 
	integer :: ms_1, mi_1, ms_2, mi_2, iel_1, iel_2, fo1, fo2
	real*8 :: alpha_1, alpha_2, x

c      double precision rle
      
      external i1st, nsum
	
cERN	21/03/05: Needed if OPT = 2 (Assume 3 regions, antenna at 2-3 interface)    
	rhoant = rx0m(2)

	k = 0
      nonec = 0
      gloty2 = .true.
      do 4 ireg = 1, nreg
c     ===================
      is1 = nsum(ireg-1, ns) + 1
      ifiel(ireg) = i1st(is1) + 1
      ifirse = ifiel(ireg)
      is2 = is1 + ns(ireg) - 1
      ilael(ireg) = i1st(is2) + iele(is2)
      ilaste = ilael(ireg)

      do 4 isubr = is1, is2
c     =====================
      eletyp = styp(isubr)
      gloty2 = gloty2 .and. eletyp.eq.'M23'
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo-icolo

      do 4 iel = i1st(isubr)+1, i1st(isubr)+iele(isubr)
c     -------------------------------------------------
c      rle = fl(iel)
      if(.not.autome)then
       if(monomo)then
      	minf(iel) = moant(1)
      	msup(iel) = moant(1)
      	nmode(iel) = 1
       else
      	minf(iel) = modva1
      	msup(iel) = modva2
cERN  21/03/05 ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c	OPT = 1: Compute modes according to surface perimeter and 'pollength'

	   	if(FLAG)then !!!! .and. OPT .eq. 1)then
		   if(ireg.eq.1)then
			  aux = min(iel+1, nele)
			  minf(iel) = max(modva1,-nint(surf_perim(ifiabs(aux))/pollength) )
			  minf(iel) = min(minf(iel),-2)
			  msup(iel) = min(modva2,+nint(surf_perim(ifiabs(aux))/pollength) )	
			  msup(iel) = max(msup(iel),+2)
		   end if
		end if	! OPT = 1

c	OPT = 2: Compute modes according to exponential function.
c			 Min/Max at antenna: modva1/modva2
c			 Min/Max at axis: modax1/modax2 (Cyrano line 256)

c		if(FLAG .and. OPT .eq. 2)then
c		   if(ireg.eq.1 .or. ireg.eq.2)then
c			  aux = min(iel+1, nele)
c			  x = abscis(ifiabs(aux))
c			  fo1 = modax1
c			  alpha_1 = 1.0d0/rhoant*log(dfloat(modva1)/dfloat(fo1)) 
c			  fo2 = modax2
c			  alpha_2 = 1.0d0/rhoant*log(dfloat(modva2)/dfloat(fo2)) 
c			  msup(iel) = nint(fo2*dexp(alpha_2*x))
c			  minf(iel) = nint(fo1*dexp(alpha_1*x))
c			  ms_2 =  msup(iel)	 
c			  mi_2 =  minf(iel)	
c		   elseif(ireg .eq. 3)then        
c		      msup(iel) = ms_2
c		      minf(iel) = mi_2 
c		   end if
c		end if	! OPT = 2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		nmode(iel) = msup(iel) + 1 - minf(iel)
       end if	! (not monomo)
      end if

      nonec = nonec + nmode(iel)**2 * (idllo**2 - icolo**2)
 
   4  continue
c     --------

      minf(1) = minf(2)
	msup(1) = msup(2)
      nmode(1) = msup(1) - minf(1) + 1
      minf(0) = minf(1)
      msup(0) = msup(1)
      nmode(0) = msup(0) - minf(0) + 1
      ithoma = 0
      nk0(1) = 0
      nkm(1) = 0
      imax = 1
      nmomax = 0
      nlca = 0
      mfilin = 0

cERNcccccccccccccccccccccccccccccccccccccc
c	do iel = 0,nele
c		print *, iel, minf(iel), msup(iel)
c	end do
ccccccccccccccccccccccccccccccccccccccccccc

      call bafaga

      ngelim = 0
              
c     write(nofile,*)'     iel          nk0        nk1   shift'

      do 1 ireg = 1, nreg
c     ===================
      is1 = nsum(ireg-1, ns) + 1
      ifirse = ifiel(ireg)
      ifim1 = max0 (1, ifirse-1)
      is2 = is1 + ns(ireg) - 1
      ilaste = ilael(ireg)
      nmodc = nmode(ifirse-1)
c     Status of axis or inner wall constraints:
c     Solver needs nonsingular diagonal block to start. Inner boundary
c     conditions are imposed either by special assembly of first element
c     (bcinb1=.f.), or by writing a block of boundary conditions and performing
c     Gaussian elimination of the resulting Lagrange multipliers before the
c     solver starts (bcinb1=.t.). From the solver point of view,
c     there is in any case no extra block for the inner boundary conditions.
      	if(ireg.gt.1 .and. nmodc.gt.0)then
c         Internal boundary condition equation block, between regions:
c 9/12/03 NB not valid at magnetic axis, because there are different numbers of
c         constraints for different poloidal modes!
      	ngelim = ngelim + ncstr(ireg-1) * nmodc
      	ithoma = ithoma + 1
      	lblock(ithoma) = ncstr(ireg-1) * nmodc
c      	write(6,*)ithoma, lblock(ithoma)
      	end if
      gloshi(ireg) = ngelim

      do 3 isubr = is1, is2
c     =====================
      eletyp = styp(isubr)
      idllo = iddl(isubr)
      icolo = iconn(isubr)
      ibulo = ibub(isubr)
      icpblo = idllo - icolo
      nmomlo = 0
      imaxlo = 0

      do 2 iel = i1st(isubr)+1, i1st(isubr)+iele(isubr)
c     -------------------------------------------------
      if( iel.eq.ifirse )then
      	if( ireg .gt. 1 )then
      	nk0(iel) = nk1(iel-1) + nmode(iel-1)
      	nkm(iel) = nkm(iel-1) + nmode(iel-1)
      	end if
      ngelim = ngelim + nmode(iel-1) * icolo
      else
      ich = nmode(iel-1)
      nk0(iel) = nk0(iel-1) + ich
      nkm(iel) = nkm(iel-1) + ich
      end if

      if(nmode(iel).gt.nmomlo)then
      imaxlo = iel
      nmomlo = nmode(iel)
      end if

      ie1 = max0(iel-1, ifirse)
      ie2 = max0(iel-2, ifirse)
      nk1(iel) = nk0(iel) + nmode(iel-1)
      nmodc = max0(nmode(iel), nmode(iel-1))
        if(nmodc.gt.0)then
        ngelim = ngelim + ibulo * nmodc
     ;                  + icolo * nmode(iel)
        ithoma = ithoma + 1
        lblock(ithoma) = nmode(iel-1)*icolo
     ;                 + max0(nmode(iel-1),nmode(iel))*ibulo
c        write(6,*)ithoma, lblock(ithoma)
        end if
c     write(nofile,100) iel, nk0(iel), nk1(iel)
 
   2  continue
      if(nmomlo .gt. nmomax)imax = imaxlo
      nmomax = max0(nmomax, nmomlo)
      nlca = max0(nlca, iddl(isubr)*nmomlo-1)
c     This fill-in for finite element problem:
      mfilin = max0(mfilin, nmomlo*(iconn(isubr)+ibub(isubr)))
   3  continue

      iel = ilael(ireg)
      isubr = is2
c       Last block in region: 'right part' of its last element
        if(nmode(iel) .gt. 0)then
        ithoma = ithoma + 1
        lblock(ithoma) = nmode(iel) * iconn(isubr)
c        write(6,*)ithoma,lblock(ithoma)
        end if
   1  continue
c     --------
      iel = nele+1
      ie1 = iel - 1
      ie2 = max0( iel-2, ifiel(nreg) )
      	if(nmode(nele).gt.0)then
c       Wall boundary condition block:
      	ngelim = ngelim + nmode(nele) * ncstr(nreg)
      	ithoma = ithoma + 1
      	lblock(ithoma) = ncstr(nreg) * nmode(nele)
c     	write(6,*)ithoma, lblock(ithoma)
      	end if
c     Number of blocks for Thomas (= block tridiagonal) algorithm:
      nthoma = ithoma
      gloshi(nreg+1) = ngelim
      nuca = nlca
      write(nofile,101) imax, nmode(imax)
      write(nofile,*) 'There are approx. ',nonec,
     ;' complex numbers in finite element blocks.'
      totvec = nk1(nele) + nmode(nele)
      totbub = nkm(nele) + nmode(nele)
      write(nofile,*)totvec,' vector blocks; vector length<= '
     ;, ndof*totvec
      write(nofile,*)'Number of unknowns =', ngelim
      write(nofile,*)'Number of Thomas blocks =', nthoma
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
cERN	Write poloidal mode couplings to file 'minf_msup.dat'
	if(FLAG) then
	  open (UNIT = 777, FILE = paplda // 'minf_msup.dat', 
     ;        STATUS = "REPLACE", ACTION = "WRITE")
		    write(777, 102)
			do iel = 0, nele
			   aux = min(iel+1, nele)
			   write(777,"(i4,f12.4,i8,i8,f12.4)"), iel, abscis(ifiabs(aux))
     ;					,minf(iel), msup(iel), surf_perim(ifiabs(aux))
		    end do
	  close(777)	
	end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
  100 format(1h ,10(i10,1x),2(d10.2,1x))
  101 format(1h ,'Element #', i4, ' has max.number of modes (', i4, ')')
  102 format(1h ,'  #      rho(m)     minf    msup     perim(m)' )  
      end
