      subroutine coordp
c
c     Writes data to plot coordinate lines
c     (to be called after TABLES)
c
      implicit none
      
      include 'pardim.copy'
      include 'comfic.copy'
      include 'comgeo.copy'
      include 'commag.copy'
      include 'comma2.copy'
      include 'comreg.copy'
      include 'comsub.copy'
      include 'comfin.copy'
      include 'comswe.copy'
      include 'comin2.copy'
      include 'comrot.copy'
      include 'compla.copy'
      include 'comfou.copy'
      include 'commod.copy'
      include 'comant.copy'
      include 'complo.copy'
      include 'complp.copy'
      include 'comwrr.copy'
      include 'comphy.copy'

      integer i, j, l, iha, nsum
      external nsum

      if(wrplta)then
c     --------------

      include 'open2d.copy'
          
      if(updsym)then
      npp = npft2 + 1
      else 
      npp = npfft + 1
      end if

c      write(nofile,*)'R, Y coordinates:'
      write(69,*)nabsci
      write(69,*)npp
        do i = 1, nabsci
          do j = 1, npp
          write(69,2000)eqt(i,j,1), eqt(i,j,2)
          end do
        end do

c      write(nofile,*)'B0, sin Theta, sin/rho:'
      write(68,*)nabsci
      write(68,*)npp
        do i = 1, nabsci
          do j = 1, npp
          write(68,2000)bmotab(i,j), eqt(i,j,15), eqt(i,j,12)
          end do
        end do

c      write(nofile,*)'sin Theta, sin/rho:'
c      do i = 1, nabsci
c      do j = 1, npp
c      write(nofile,2000)eqt(i,j,15),eqt(i,j,12)
c      end do
c      end do

      if(cokpco)then
c      write(nofile,*)'thbar,phibar,mu:'
      write(67,*)nabsci
      write(67,*)npp
        do i = 1, nabsci
          do j = 1, npp
c         write(67,2000)ckt(i,j,3), ckt(i,j,4), ckt(i,j,5)
          write(67,2001)(ckt(i,j,l),l=1,7)
          end do
        end do
	  end if
	  
      end if
c     ------
	  
      if(wrifou)then
c     --------------
      include 'openfo.copy'

	  do i = 1, nabsci, nabsci-1
	  intab = i
	  y = abscno(i)
	    if(i.eq.1)then
	    isubr = 1
	    else
	    isubr = nsum(nreg, ns)
	    end if
	  call foucog(1)
c	  write(nofile,*)'Fourier coeffs. at abscissa #',i
c	  write(nofile,*)'================================='
c	    do j = 1, nfouwr
c	    j2 = ifouwr(j)
c	    write(nofile,*)'Coeff. #',j2
      write(66,*)npfft+1
      write(66,*)nfouwr
      write(66,*)(ifouwr(j),j=1,nfouwr)
	    do iha = -npft2, npft2
	    write(66,1000)dfloat(iha)
     ;, (dreal(cvc(iha,ifouwr(j))), dimag(cvc(iha,ifouwr(j))),j=1,nfouwr)
	    end do
c	    end do
	  end do
      end if
c     ------
	  
      return

 1000 format(1h ,8(g14.6, 2x))
 2000 format(1h ,3(2x,g14.6))
 2001 format(1h ,7(2x,g14.6))

      end