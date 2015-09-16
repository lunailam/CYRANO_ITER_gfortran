c
      subroutine fbtini(iuout
     ;                 ,iumat ,newfm, labelm, lenlm
     ;                 ,iusol ,newfs, labels, lenls
     ;                 ,iuema ,newfe, labele, lenle
     ;                 ,doc, stomat)

      implicit none
c
c     Initializes frontal block-tridiagonal solver
c
c     Author: Philippe Lamalle
c     1/10/92 Version
c     24/10/2000 Modified to allow storage of element matrices
c
      include 'pardim.copy'
      include 'comzzz2.copy'
      include 'comusr.f'
c
      integer iimin, iimax
     ;      , iuout, iumat, iusol, iuema, newfm, newfs, newfe, doc
     ;      , lenlm, lenls, lenle
      double precision labelm(*), labels(*), labele(*)
      logical stomat
C-----------------------------------------------------------------------
C   IUOUT: output file number (messages, results)(input).
C   IUMAT: Fortran unit number for storage of global matrix
C          factorization (input).
C   IUSOL: Fortran unit number for storage of solutions (input).
C   IUEMA: Fortran unit number for storage of element matrices (input).
C   NEWFM: nonzero: use a new file to store matrix factorizations.
C          zero:    use a file already containing at least one previous
C                   factorization (input).
C   LABELM: real array, at user disposal to label the file of
C           factorizations (input when NEWFM is nonzero, else output).
C   LENLM:  length of LABELM (input when NEWFM is nonzero, else output).
C   NEWFS: nonzero: use a new file to store solutions.
C          zero:    use a file already containing at least one previous
C                   solution (input).
C   LABELS: real array, at user disposal to label the file of
C           solutions (input when NEWFS is nonzero, else output).
C   LENLS:  length of LABELS (input when NEWFS is nonzero, else output).
C   NEWFE: nonzero: use a new file to store element matrices.
C          zero:    use a file already containing at least one previous
C                   set of matrices (input).
C   LABELE: real array, at user disposal to label the file of element
C           matrices (input when NEWFE is nonzero, else output).
C   LENLE:  length of LABELE (input when NEWFE is nonzero, else output).
C
C   DOC: (integer) when nonzero, documentation on execution will be
C        written to Fortran unit IUOUT (input).
C   STOMAT: logical. true: store element matrices (uses additional disk storage)
C           NB: this is provided for the user's convenience. 
C           Such storage is not required by the solver.
C                    false: don't store element matrices (and save disk space)
C-----------------------------------------------------------------------
      integer nbpr, i, j, nblk, iprob, ceilq
      external iimin, iimax, ceilq
      
      nbpr = ceilq(2*(maxpro+2),iobll)

        if(stomat)then
        nuni = 3
        else
        nuni = 2
        end if
c
c     Output file number:
      nofile=iuout
C     Fortran unit number for storage of global matrix factorization:
      unit(1)=iumat
C     Fortran unit number for storage of solutions:
      unit(2)=iusol
C     Fortran unit number for storage of element matrices:
      unit(3)=iuema
c
      newfil(1) = newfm .ne. 0
      newfil(2) = newfs .ne. 0
      newfil(3) = newfe .ne. 0
      soldoc = doc .ne. 0
      if(soldoc)write(nofile,100)
      if(unit(1) .eq. unit(2))then
      write(nofile,*)'Error in routine FBTINI: matrix and solution'
     ;,' may not be stored on the same FORTRAN unit'
      stop
      end if
      if(stomat .and. (unit(2) .eq. unit(3)
     ;            .or. unit(1) .eq. unit(3)))then
      write(nofile,*)
     ;'Error in routine FBTINI: element matrices, factorization and solution'
     ;,' may not be stored on the same FORTRAN unit'
      stop
      end if

      do 1 i = 1, nuni
c     Open I/O units for matrix and solutions:
      call asop(i)
      iwri(i) = 0
      irea(i) = 0
      nblrd(i) = 0
      nblwr(i) = 0
      do 10 j = 1, lwk1
      wk1zzz(j,i) = 0.
  10  continue

c     wk1zzz to receive file heading information
c     iprdef 

      if(newfil(i))then
c     *****************
      nprobl(i) = 0
      nblk = 0
      if(i.eq.1 .and. lenlm.gt.0)then
      iblk(i) = 1 + nbpr
      wk1zzz(3,i) = lenlm
      call aswr(i, lenlm, labelm, 0, nblk)
      call aswa(i)
      end if
      if(i.eq.2 .and. lenls.gt.0)then
      iblk(i) = 1 + nbpr
      wk1zzz(3,i) = lenls
      call aswr(i, lenls, labels, 0, nblk)
      call aswa(i)
      end if
      if(i.eq.3 .and. lenle.gt.0)then
      iblk(i) = 1 + nbpr
      wk1zzz(3,i) = lenle
      call aswr(i, lenle, labele, 0, nblk)
      call aswa(i)
      end if
      iprdef(1,i) = 1 + nbpr + nblk
c
      else
c     ****
c     Fortran unit UNIT(I) already contains at least 1 factorization.
      iblk(i) = 1
      call asrd(i, iobll, wk1zzz(1,i), 0, nblk)
      iblk(i) = iblk(i) + nblk
      call aswa(i)
      nprobl(i) = wk1zzz(1,i)
      nblk = wk1zzz(2,i)-1
          if(nprobl(i).gt.maxpro)then
          write(nofile,*)'Warning: Fortran unit ',UNIT(I)
     ;,' contains ',NPROBL(I),' problems,'
     ;,' whereas this version of solver allows a maximum of ',MAXPRO
      write(nofile,*)
     ;'Only the first ',MAXPRO,' problems will be considered'
      write(nofile,*)
     ;'To cure this, recompile solver with parameter MAXPRO = ',MAXPRO
      nprobl(i) = maxpro
      nblk = nbpr
          end if
      if(nblk.gt.0)then
      call asrd(i, nblk*iobll, wk1zzz(iobll+1,i), 0, nblk)
      iblk(i) = iblk(i) + nblk
      call aswa(i)
      end if

      do 50 iprob = 1, nprobl(i)
c     I-O block index of problem heading:
      iprdef(iprob,i) = wk1zzz(2*iprob+2,i)
c     I-O block index of problem data:
      iprdat(iprob,i) = wk1zzz(2*iprob+3,i)
  50  continue
      iprdef(nprobl(i)+1,i) = wk1zzz(2*nprobl(i)+4,i)

      if(wk1zzz(3,i).gt.0)then
      iblk(i) = 1 + nbpr
          if(i.eq.1)then
          lenlm = wk1zzz(3,i)
          call asrd(i, lenlm, labelm, 0, nblk)
          else if(i.eq.2)then
          lenls = wk1zzz(3,i)
          call asrd(i, lenls, labels, 0, nblk)
          else if(i.eq.3)then
          lenle = wk1zzz(3,i)
          call asrd(i, lenle, labele, 0, nblk)
          end if
      call aswa(i)
      end if
      end if
c     ******
   1  continue

c
      return
c
 100  format(1h1,'Frontal block-tridiagonal solver'/1H ,
     ;           '********************************'/)
      end
c
      subroutine proini(nblo, lblo, nrhss, new, ipro)

      implicit none
c
c     Initializes one problem
c
c     Author: Philippe Lamalle
c     1/10/92 Version
c     Modified 24/10/2000 to allow storage of element matrices

      include 'pardim.copy'
      include 'comzzz2.copy'
      include 'comusr.f'
c
      integer nblo, lblo(nblo), nrhss, new, ipro
C-----------------------------------------------------------------------
C   NBLO: row- or column-wise number of blocks in global matrix (input
C         when NEW nonzero).
C   LBLO: array of block sizes (input when NEW nonzero).
C   NRHSS: number of right-hand-sides (input).
C   NEW: nonzero: a new factorization has to be performed
C        zero:    a previously computed factorization has to be used
C        (input)
C   IPRO: index of old factorization to be used (input only for NEW=0)
C-----------------------------------------------------------------------
      integer 
     ;  i, j, nblk, ierror, nwr
     ;, isum, imin, imax, iimin, iimax, ceilq
     ;, int
      external iimin, iimax, ceilq
      intrinsic int
c
      ierror = 0
      if(new.eq.0 .and. nrhss.eq.0)then
      write(nofile,*)'From PROINI: no r.h.s.; factorization already'
     ;//' exists: end of problem'
      return
      end if
      newmat(1) = new .ne. 0
      newmat(2) = nrhss .gt. 0
c     Case element matrices to be stored: 
c     assume same status as factorization file:
      if(nuni.eq.3)NEWMAT(3) = NEWMAT(1)
      do 1 i = 1, nuni
      do 80 j = 1, lwk2
      wk2zzz(j,i) = 0.
  80  continue
      if(newmat(i))then
c     *****************
c     New problem
      iprobl(i) = nprobl(i) + 1
      nthoma = nblo
      wk2zzz(1,i) = nthoma
      wk2zzz(2,i) = ceilq(nthoma+5,iobll)
      wk2zzz(3,i) = modele
      do 60 ithoma = 1, nthoma
      lblock(ithoma) = lblo(ithoma)
      wk2zzz(ithoma+5,i) = lblo(ithoma)
  60  continue
      else
c     ****
C     I=1: Matrix factorization is already on unit UNIT(I).
C     I=2: No r.h.s. provided.
C     I=3: Element matrices are already on unit UNIT(I).
          if(modfac.ne.1)then
      write(nofile,*)'error: to work with a previously stored matrix,'
      write(nofile,*)'the value of parameter modfac must be set to 1'
      stop
          end if
          if(i.eq.1)then
      if(ipro.le.0 .or. ipro.gt.nprobl(i))then
      write(nofile,*)'error in proini: the old factorization # ipro '
     ;,'must be between 1 and',nprobl(i)
      write(nofile,*)'but ipro=',ipro,' was given'
      ierror = 1
      end if
      iprobl(i) = ipro
      iblk(i) = iprdef(iprobl(i),i)
      call asrd(i, iobll, wk2zzz(1,i), 0, nblk)
      iblk(i) = iblk(i) + nblk
      call aswa(i)
          if(int(wk2zzz(3,i)).ne.modele)then
          write(nofile,110)unit(i),int(wk2zzz(3,i)),modele
          ierror = 1
          end if
      nthoma = wk2zzz(1,i)
      tstore = wk2zzz(4,i)
          if(wk2zzz(2,i).gt.1)then
          nblk = wk2zzz(2,i) - 1
          call asrd(i, nblk*iobll, wk2zzz(iobll+1,i), 0,
     ;              nblk)
          iblk(i) = iblk(i) + nblk
          call aswa(i)
          end if
      do 50 ithoma = 1, nthoma
      lblock(ithoma) = wk2zzz(ithoma+5,i)
  50  continue
          end if
      end if
c     ******
   1  continue

c     NB: this is nrhs of comzzz2:
      nrhs = nrhss
      if(nthoma.le.0)then
      write(nofile,200)nthoma
      ierror = 1
      else if(nthoma.gt.maxnbl)then
      write(nofile,300)nthoma,nthoma
      ierror = 1
      end if
      if(nrhs.lt.0)then
      write(nofile,400)nrhs
      ierror = 1
      else if(nrhs.gt.maxrhs)then
      write(nofile,500)nrhs,nrhs
      ierror = 1
      else if(nrhs.eq.0 .and. modfac.eq.0)then
      write(nofile,*)'NRHSS=0 is inconsistent with parameter MODFAC=0'
      ierror = 1
      end if
      if(ierror.ne.0)stop

c     Maximal and minimal size block indices:
c      write(6,*)'calling iimax: nthoma,lblock=',nthoma,lblock
      imax = iimax(nthoma,lblock,1)
c      write(6,*)'calling iimin: nthoma,lblock=',nthoma,lblock
      imin = iimin(nthoma,lblock,1)

      if(lblock(imin).le.0)then
      write(nofile,*)'block #',imin,' was given a negative size'
      ierror = 1
      end if
      if(lblock(imax).gt.maxbll)then
      write(nofile,600)imax,lblock(imax)
      ierror = 1
      end if
      if(ierror.ne.0)stop

c     Number of unknowns:
      ngelim = isum(nthoma,lblock,1)

C     Number of stored matrix and r.h.s. entries:
      stored = modfac*lblock(1)**2 + ngelim * nrhs
C     Number of disk blocks needed to store matrix and pivots:
      tstore = modfac*ceilq(lblock(1)*(nwn*lblock(1)+1),iobll)
C     Number of disk blocks needed to store element matrices and rhs:
      storel = 0
C     Number of disk blocks needed to store r.h.s.:
      blorhs = ceilq(lblock(1)*nrhs*nwn,iobll)
      do 10 i = 2, nthoma
      stored = stored + lblock(i)*(modfac*lblock(i)
     ;                             +(modfac+1)*lblock(i-1))
      tstore = tstore + modfac*(ceilq(lblock(i)*(nwn*lblock(i)+1),iobll)
     ;                + ceilq(lblock(i-1)*lblock(i)*nwn,iobll) )
     ;+ ceilq( lblock(i-1)*(lblock(i)+(1-modfac)*nrhs)*nwn, iobll)
      storel = storel + ceilq(nwn*(lblock(i)+lblock(i-1))*
     ;(lblock(i)+lblock(i-1)+nrhs),iobll)
      blorhs = blorhs + ceilq(lblock(i)*nrhs*nwn,iobll)
  10  continue
      if(modfac.eq.0)tstore = tstore
     ;+ ceilq(lblock(nthoma)*nrhs*nwn, iobll)

      if(soldoc)then
      write(nofile,*)' '
          if(nrhs.gt.0)then
      write(nofile,*)'Factorization # ',IPROBL(1),'; solution group # '
     ;,iprobl(2),':'
      write(nofile,*)'*********************************************************'
          else
      write(nofile,*)'Factorization # ',IPROBL(1)
      write(nofile,*)'*************************'
          end if
      write(nofile,*)' '
          if(newmat(1))then
      write(nofile,*)'Number of blocks on diagonal:',NTHOMA
      write(nofile,*)'Maximum block size:          ',LBLOCK(IMAX)
      write(nofile,*)'Index of maximum size block: ',IMAX
      write(nofile,*)'Number of unknowns:          ',NGELIM
      if(nwn.eq.2)then
      write(nofile,*)
     ;'Matrix and r.h.s. blocks contain a total of ',FLOAT(STORED)
     ;,' complex numbers'
      else if(nwn.eq.1)then
      write(nofile,*)
     ;'Matrix and r.h.s. blocks contain a total of ',FLOAT(STORED)
     ;,' real numbers'
      end if

      write(nofile,*)'List of block sizes:'
      write(nofile,*)'--------------------'
      nwr = ceilq(nthoma, 6)
        do 20 i = 1, nwr
        write(nofile,800)(j+i-1,lblock(j+i-1),j=1,nthoma+1-i,nwr)
  20    continue
      write(nofile,1000) tstore, iobll, iobll*1.d-6*tstore
     ;, iobll*wordlzz*1.d-6*tstore
 1000 format(1h , 'Global matrix and pivots storage requires   ', i7 
     ;, ' blocks of', i4, ' words, that is ', g10.3, ' Megaword, or '
     ;, g10.3, ' Megabyte.')        
        if(nuni.eq.3)then
        write(nofile,1001) storel, iobll, iobll*1.d-6*storel
     ; , iobll*wordlzz*1.d-6*storel
 1001 format(1h , 'Element matrices and r.h.s. storage requires', i7 
     ;, ' blocks of', i4, ' words, that is ', g10.3, ' Megaword, or '
     ;, g10.3, ' Megabyte.')        
        end if
          end if

      if(nrhs .gt. 0)
     ;write(nofile,1002) blorhs, iobll, iobll*1.d-6*blorhs
     ;, iobll*wordlzz*1.d-6*blorhs
 1002 format(1h , 'Solution storage requires                   ', i7 
     ;, ' blocks of', i4, ' words, that is ', g10.3, ' Megaword, or '
     ;, g10.3, ' Megabyte.')        
      end if

c     Initialize common COMUSR.F to zero:
c	call dset(luszzz, 0.d0, zzzzzz, 1)
      do 30 i = 1, luszzz
      zzzzzz(i) = 0.0d0
c      ctom1(i) = 0.0d0
  30  continue

      nblk = ceilq(nthoma+5, iobll)
      if(newmat(1))iprdat(iprobl(1),1) = iprdef(iprobl(1),1) + nblk
      iblk(1) = iprdat(iprobl(1),1)
      if(newmat(2))then
      iprdat(iprobl(2),2) = iprdef(iprobl(2),2) + nblk
      iblk(2) = iprdat(iprobl(2),2)
      end if
      if(newmat(3))then
      iprdat(iprobl(3),3) = iprdef(iprobl(3),3) + nblk
      iblk(3) = iprdat(iprobl(3),3)
      end if

c     Information on problem size, will be written at successful
c     completion:
      wk2zzz(1,1) = nthoma
      wk2zzz(2,1) = nblk
      wk2zzz(3,1) = modele
      wk2zzz(4,1) = tstore
      wk2zzz(1,2) = nthoma
      wk2zzz(2,2) = nblk
      wk2zzz(3,2) = modele
      wk2zzz(4,2) = blorhs
      wk2zzz(5,2) = nrhs
      if(nuni.eq.3)then
      wk2zzz(1,3) = nthoma
      wk2zzz(2,3) = nblk
      wk2zzz(3,3) = modele
      wk2zzz(4,3) = blorhs
      wk2zzz(5,3) = nrhs
      end if
      do 40 i = 1, nuni
      do 40 ithoma = 1, nthoma
      wk2zzz(ithoma+5,i) = lblock(ithoma)
  40  continue

c     Initialize block counter:
      ithoma = 0

      return

 110  format(1h ,'The old factorization on Fortran unit',I4,
     ;' was obtained with solver storage mode MODELE=',I4,'.'
     ;/1H ,'To use that factorization, recompile solver with parameter',
     ;' MODELE=',I4)
 200  format(1h ,'The number of blocks must be positive, but NTHOMA=',
     ;i4,' was given')
 300  format(1h ,'The number of blocks,',I4,
     ;', is larger than the allowed maximum.'/1H ,
     ;' Recompile solver setting MAXNBL=',I4)
 400  format(1h ,
     ;'The number of right hand sides must be positive, but NRHS=',
     ;i4,' was given')
 500  format(1h ,'The number of right hand sides,',I4,
     ;', is larger than the allowed maximum.'/1H ,
     ;' Recompile solver setting MAXRHS=',I4)
 600  format(1h ,'The size of block #',I4,
     ;' is larger than the allowed maximum.'/1H ,
     ;' Recompile solver setting MAXBLL=',I4)
 800  format(1h ,6(i4,', ',i4,' ; '))

      end
c
      subroutine fbstep(scain, sc1, sc2)

      implicit none
c
      include 'pardim.copy'
      include 'comzzz2.copy'
      include 'comusr.f'
c
      integer scain
      double precision    sc1(maxbll), sc2(maxbll)
c
c     Author: Philippe Lamalle
c     1/10/92 Version
c     Modified 24/10/2000 to allow storage of element matrices
c
C  -  When a new matrix is undergoing factorization,
C     this routine performs one step of the frontal block-tridiagonal
C     algorithm. At the end of the last step, it also performs back
C     substitution on the right hand sides provided by the user during
C     matrix assembly.
C     The solver performs scaled row pivoting:
C     SCAIN: integer (input).
C     SCAIN = 0: inverse scaling factors are the rows maximum norms
C     (output).
C     SCAIN = 1: row inverse scaling factors are input by the user
C                in SC1.
C     SCAIN not 0 or 1: row and column inverse scaling factors are
C     input by the user in SC1 and SC2, respectively.
C
C     SC1: real vector of length LBLOCK(ithoma) containing the rows
C     inverse scaling factors used in partial pivoting
C     (output when SCAIN=0, input otherwise).
C
C     SC2: real vector of length LBLOCK(ithoma) containing column
C          inverse scaling factors when SCAIN differs from 0 and 1.
C          dummy argument when SCAIN = 0 or 1.
C          SC2 influences the evaluation of the block condition number.
C
C  -  When a previous factorization is used, FBSTEP performs one step
C     of the forward substitution on the new right hand side block.
C     In this situation, the arguments of FBSTEP are all dummy.
C     At the end of the last step, back-substitution is performed and
C     the solution vectors are stored on disk.
C
      integer i, j, k, j2, iii, lbl1, lbl2, lbl3, lbl4, ising, jshi
     ;, ncw, iqueue, nblk, nwr 
     ;, info
     ;, izamax, ceilq 
     ;, int
      double precision  cabs1, rcond(maxnbl)
      complex*16 czero, zdum
c
      intrinsic int
      external ceilq, izamax
c
      data czero/(0.0d0,0.0d0)/
      cabs1(zdum) = dabs(dreal(zdum))+dabs(dimag(zdum))


      ithoma = ithoma + 1
      lbl2 = lblock(ithoma)
      if(ithoma.lt.nthoma)then
      lbl3 = lblock(ithoma+1)
      else
      lbl3 = 0
      end if
      lbl4 = lbl3 + nrhs
      if(ithoma.gt.1)then
      lbl1 = lblock(ithoma-1)
      else
      lbl1 = 0
      end if
      jshi = modele*lbl1 + lbl2 + lbl3

              if(newmat(1))then
c             *****************
c     Factorize a new matrix

      if(nuni.eq.3)then
c     This is the case where the user requires saving element matrices 
c     (solver init. FBTINI called with stomat=.true.)
c     Must store current element matrix and associated r.h.s., then add
c     overlapping block of former element and rhs that were stored in CTOM.
c     After that, restore contents of CTOM and proceed normally.
c
c     Pack and write element matrix and r.h.s.:
c     This is split into 3 instructions to ease retrieval of blocks 
c     to be assembled.

c     Case MODELE=0:
      call packma(lbl2+lbl3, lbl2+lbl3, ael, ldael, ael, 1)
      ncw = (lbl2 + lbl3) ** 2
        if(ncw .gt. 0)then
        call aswr(3, nwn*ncw, ael, 1, nblk)
        iblk(3) = iblk(3) + nblk
        end if
      call packma(lbl2+lbl3, nrhs, bel, ldbel, bel, 1)
      ncw = (lbl2 + lbl3) * nrhs
        if(ncw .gt. 0)then
        call aswr(3, nwn*ncw, bel, 0, nblk)
        iblk(3) = iblk(3) + nblk
        end if
c     wait for completion, then unpack element matrix and r.h.s.:
      call aswa(3)
      call packma(lbl2+lbl3, lbl2+lbl3, ael, ldael, ael, 2)
      call packma(lbl2+lbl3, nrhs, bel, ldbel, bel, 2)
      end if

      if(ithoma.gt.1)then
c
      if(nuni.eq.3)then
c       Overlapping parts of former element and rhs have been stored in CTOM at
c       previous call. Add them to current element matrix and r.h.s.:
        do j = 1, lbl2
          do i = 1, lbl2
          ael(i,j) = ael(i,j) + ctom(i,j)
          end do
        end do
        do j = 1, nrhs
          do i = 1, lbl2
          bel(i,j) = bel(i,j) + ctom(i,j+lbl2)
          end do
        end do
c     Now, restore CTOM to correct value by reading on factorization unit:
c     Read block CTOM and start asynchronous I-O request:
c
      if(lbl2+nrhs .gt. 0)then
        if(modfac .eq. 0)then
        ncw = lbl1 * (lbl2 + nrhs)
c       Skip back & read:
        iblk(1) = iblk(1) - ceilq(nwn*ncw,iobll)
        call asrd(1, nwn*ncw, ctom, 1, nblk)
        iblk(1) = iblk(1) + nblk
        else
          if(lbl2 .gt. 0)then
          ncw = lbl1 * lbl2
c         Skip back & read:
          iblk(1) = iblk(1) - ceilq(nwn*ncw,iobll)
          call asrd(1, nwn*ncw, ctom, 1, nblk)
          iblk(1) = iblk(1) + nblk
          end if
          if(nrhs .gt. 0)then
          ncw = lbl1 * nrhs
c         Skip back & read:
          iblk(2) = iblk(2) - ceilq(nwn*ncw,iobll)
          call asrd(2, nwn*ncw, ctom(lbl1*lbl2+1,1), 1, nblk)
          iblk(2) = iblk(2) + nblk
          end if
        end if
      end if
      
      end if

c     Wait for completion of previous writes on factorization and r.h.s. units:
      call aswa(1)
      if(modfac.eq.1 .and. nrhs.gt.0)call aswa(2)

c     Unpack previously written BTOM, CTOM:
      if(modfac.eq.1 .and. modele.eq.0)
     ;call packma(lbl2, lbl1, btom, ldbt, btom, 2)
      call packma(lbl1, lbl2+nrhs, ctom, ldct, ctom, 2)
      end if

      if(scain .eq. 0)then
c     Compute row max.norms using contents of AEL and BTOM:
      do 10 i = 1, lbl2
      sc1(i) = cabs1( ael(i,izamax(jshi, ael(i,1), ldael) ) )
  10  continue
          if( modele.eq.0 .and. ithoma .gt. 1 )then
          do 11 i = 1, lbl2
          sc1(i) = dmax1( sc1(i), cabs1(btom(i,
     ;                    izamax(lbl1, btom(i,1), ldbt) ) ) )
  11      continue
          end if
      end if

c     Check no row is zero:
      ising = 0
      do 12 i = 1, lbl2
        if(sc1(i) .eq. 0.0d0) then
        ising = i
	  go to 100
        end if
  12  continue
 100  continue
        if(ising.ne.0)then
        write(nofile,*)
     ;  'ROW #', ISING, ' IN BLOCK #', ITHOMA,' IS ZERO; ROUTINE FBSTEP STOPPED'
        stop
        end if
c
c     One step of f.b.t.algorithm:
c
      if( ithoma .gt. 1 )then
c
c     Compute new ALPHA:
c
          if(modele .eq. 0)then
      do 13 j = 1, lbl2
      do 13 k = 1, lbl1
      do 13 i = 1, lbl2
      ael(i,j) = ael(i,j) - btom(i,k) * ctom(k,j)
  13  continue
      do 14 j = 1, nrhs
      do 14 k = 1, lbl1
      do 14 i = 1, lbl2
      bel(i,j) = bel(i,j) - btom(i,k) * ctom(k,j+lbl2)
  14  continue
          else
      do 15 j = 1, lbl2
      do 15 k = 1, lbl1
      do 15 i = 1, lbl2
      ael(i,j+lbl1) = ael(i,j+lbl1) - ael(i,k) * ctom(k,j)
  15  continue
      do 16 j = 1, nrhs
      do 16 k = 1, lbl1
      do 16 i = 1, lbl2
      bel(i,j) = bel(i,j) - ael(i,k) * ctom(k,j+lbl2)
  16  continue
          end if
      end if

c     Factorize alpha and store in atom:

      iii=1
      if(scain.ne.0 .or. scain.ne.1)iii=2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call cgeco2( ael(1,1+modele*lbl1), ldael, lbl2,
     ;             ipvt, iii, sc1, sc2, rcond(ithoma), info, wkczzz )
c     Faster: to avoid condition number evaluation, replace two previous
C     lines by four following ones:
c     call cgefa2( ael(1,1+modele*lbl1), ldael, lbl2,
c    ;             ipvt, sc1, 1, info)
c     rcond(ithoma) = 1.
c     if(info .ne. 0)rcond(ithoma) = 0.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      if(info .ne. 0 .or. 1.+rcond(ithoma) .eq. 1.) then
c 6/8/98: don't exit on rcond = 0.:
      if(info .ne. 0) then
      write(nofile,*)
     ;'BLOCK ALPHA #',ITHOMA,' IS SINGULAR; SOLVER STOPPED; INFO =', INFO
      stop
      end if
      if(1.+rcond(ithoma) .eq. 1.) then
      write(nofile,*)
     ;'BLOCK ALPHA #',ITHOMA,': warning: RCOND = 0'
      end if

      if(modfac .eq. 1)then
c     Pack and write block ALPHA and pivoting sequence:
      call packma(lbl2, lbl2,  ael(1,1+modele*lbl1), ldael, atom, 1)
      ncw = lbl2**2
      j = nwn*ncw
      do 17 i = 1, lbl2
      atom1(i+j) = ipvt(i)
  17  continue
      iqueue = 1
      if(ithoma.ge.nthoma .and. nrhs.eq.0)iqueue = 0
      call aswr(1, j+lbl2, atom, iqueue, nblk)
      iblk(1) = iblk(1) + nblk
      end if
c
c     Solve ATOM * GAMMA = CTOM, store in CTOM:
c
      if(lbl3 .gt. 0)then
      do 18 j = 1, lbl3
      j2 = j + lbl2 + modele*lbl1
      do 18 i = 1, lbl2
      ctom(i,j) = ael(i,j2)
  18  continue
      end if
      if(nrhs .gt. 0)then
      do 19 j = 1, nrhs
      j2 = j + lbl3
      do 19 i = 1, lbl2
      ctom(i,j2) = bel(i,j)
  19  continue
      end if
      if(lbl4 .gt. 0)
     ;call cgesl2( ael(1,1+modele*lbl1), ldael, lbl2, ipvt,
     ;             ctom, ldct, lbl4)

      if( ithoma.lt.nthoma .or. modele.eq.1 )then
c
          if(modfac.eq.1)then
c         Pack B in btom:
          call packma(lbl3*(1-modele)+lbl2*modele, lbl2*(1-modele)
     ;          +lbl1*modele, ael(lbl2*(1-modele)+1,1), ldael, btom, 1)
          else if(modele.eq.0)then
          do 20 j = 1, lbl2
          do 20 i = 1, lbl3
          btom(i,j) = ael(lbl2+i,j)
  20      continue
          end if
      end if
c
c     Write block BTOM:
c
      if(modfac .eq. 1)then
          if(modele.eq.0)then
          ncw = lbl3 * lbl2
          else
          ncw = lbl2 * lbl1
          end if
          if(ncw .gt. 0)then
          call aswr(1, nwn*ncw, btom, 1, nblk)
          iblk(1) = iblk(1) + nblk
          end if
      end if
c
c     Write block CTOM and start asynchronous I-O request:
c
      if(lbl4 .gt. 0)then
      call packma(lbl2, lbl4, ctom, ldct, ctom, 1)
          if(modfac .eq. 0)then
          ncw = lbl2 * lbl4
          call aswr(1, nwn*ncw, ctom, 0, nblk)
          iblk(1) = iblk(1) + nblk
          else
              if(lbl3 .gt. 0)then
              ncw = lbl2 * lbl3
              call aswr(1, nwn*ncw, ctom, 0, nblk)
              iblk(1) = iblk(1) + nblk
              end if
              if(nrhs .gt. 0)then
              ncw = lbl2 * nrhs
              call aswr(2, nwn*ncw, ctom(lbl2*lbl3+1,1), 0, nblk)
              iblk(2) = iblk(2) + nblk
              end if
          end if
      end if

      if(nuni.eq.3)then
c     Use CTOM as temporary storage for lower right block of AEL and lower BEL:
c     First wait for write completion:
      call aswa(1)
      if(modfac.eq.1 .and. nrhs.gt.0)call aswa(2)
c     Copy:
        do j = 1, lbl3
          do i = 1, lbl3
          ctom(i,j) = ael(i+lbl2,j+lbl2)
          end do
        end do
        do j = 1, nrhs
          do i = 1, lbl3
          ctom(i,j+lbl3) = bel(i+lbl2,j)
          end do
        end do
      end if

c 24/10/2000      IF(MODELE.EQ.0 .AND. ITHOMA.LT.NTHOMA)THEN
      if(nuni.ne.3 .and. modele.eq.0 .and. ithoma.lt.nthoma)then
c
c     Shift contributions for next blocks and clean up:
c
      do 21 j = 1, lbl3
      j2 = j + lbl2
      do 21 i = 1, lbl3
      ael(i,j) = ael(i+lbl2,j2)
  21  continue
      do 22 j = 1, nrhs
      do 23 i = 1, lbl3
      bel(i,j) = bel(i+lbl2,j)
  23  continue
      do 22 i = lbl3+1, lbl3+lbl2
      bel(i,j) = czero
  22  continue
      do 24 j = 1, lbl3
      do 24 i = lbl3+1, lbl3+lbl2
      ael(i,j) = czero
  24  continue
c	This produces error in debug mode with array bound checking:
c      do 25 i = 1, lbl2*ldael
c      ael(i,lbl3+1) = czero
c  25  continue
      do 25 j = 1, lbl2
      do 25 i = 1, ldael
      ael(i,lbl3+j) = czero
  25  continue
c
      else
c
c     Only clean up:
c
c	This produces error in debug mode with array bound checking:
c      do 26 i = 1, jshi*ldael
c      ael(i,1) = czero
c  26  continue
c      do 27 i = 1, nrhs*ldbel
c      bel(i,1) = czero
c  27  continue
	do 26 j = 1, jshi
      do 26 i = 1, ldael
      ael(i,j) = czero
  26  continue
      do 27 j = 1, nrhs
      do 27 i = 1, ldbel
      bel(i,j) = czero
  27  continue
c
      end if
c
c
      if(ithoma.lt.nthoma)return
c
      if(soldoc)then
      write(nofile,*)' '
      write(nofile,*)'List of block condition numbers:'
      write(nofile,*)'--------------------------------'
      nwr = ceilq(nthoma, 4)
      do 28 i = 1, nwr
      write(nofile,200)(j+i-1,1./rcond(j+i-1),j=1,nthoma+1-i,nwr)
  28  continue
 200  format(1h ,4(i4,', ',g10.2,' ; '))
      end if

c     Wait for completion of last write:
      call aswa(1)

c
c     Successful factorization: update indices
c
      if(modfac .eq. 1)then
      nprobl(1) = nprobl(1) + 1
      iblk(1) = iprdef(nprobl(1),1)
      call aswr(1, nthoma+5, wk2zzz, 1, nblk)
      wk1zzz(1,1) = nprobl(1)
      wk1zzz(2,1) = ceilq(2*nprobl(1)+4,iobll)
      wk1zzz(2*nprobl(1)+3,1) = iprdat(nprobl(1),1)
      iprdef(nprobl(1)+1,1) = iprdat(nprobl(1),1) + tstore
      wk1zzz(2*nprobl(1)+4,1) = iprdef(nprobl(1)+1,1)
      iblk(1) = 1
      call aswr(1, 2*nprobl(1)+4, wk1zzz, 0, nblk)
      call aswa(1)
      if(nrhs.gt.0)call aswa(2)
      end if
c
      if(nuni .eq. 3)then
      nprobl(3) = nprobl(3) + 1
      iblk(3) = iprdef(nprobl(3),3)
      call aswr(3, nthoma+5, wk2zzz(1,3), 1, nblk)
      wk1zzz(1,3) = nprobl(3)
      wk1zzz(2,3) = ceilq(2*(nprobl(3)+2),iobll)
      wk1zzz(2*nprobl(3)+3,3) = iprdat(nprobl(3),3)
      iprdef(nprobl(3)+1,3) = iprdat(nprobl(3),3) + storel
      wk1zzz(2*nprobl(3)+4,3) = iprdef(nprobl(3)+1,3)
      iblk(3) = 1
      call aswr(3, 2*nprobl(3)+4, wk1zzz(1,3), 0, nblk)
      call aswa(3)
      end if
c
      if(nrhs.le.0)return
c
              else
c             ****
c     Use a previously computed factorization
c
      if( ithoma .gt. 1 )then
c
c     Read block B:
c
      ncw = lbl2 * lbl1
      call asrd(1, nwn*ncw, btom, 0, nblk)
      call aswa(1)
      call packma(lbl2, lbl1, btom, ldbt, btom, 2)
      if(modele .eq. 0)iblk(1) = iblk(1) + nblk
c
c     Modify r.h.s.:
c
      call aswa(2)
      call packma(lbl1, nrhs, ctom, ldct, ctom, 2)
      do 29 j = 1, nrhs
      do 29 k = 1, lbl1
      do 29 i = 1, lbl2
      bel(i,j) = bel(i,j) - btom(i,k) * ctom(k,j)
  29  continue

          if(modele .eq. 0)then
c
c         Skip block C:
c
          ncw = lbl1 * lbl2
          nblk = ceilq(nwn*ncw,iobll)
          iblk(1) = iblk(1) + nblk
          else
c
c         Skip backward on A:
c
          iblk(1) = iblk(1) - ceilq((nwn*lbl2+1)*lbl2,iobll)
          end if
      end if
c
c     Read block ALPHA and pivoting sequence:
c
      ncw = lbl2**2
      j = nwn * ncw
      call asrd(1, j+lbl2, atom1, 0, nblk)
      call aswa(1)
      do 30 i = 1, lbl2
      ipvt(i) = int( atom1(i+j) )
  30  continue
      call packma(lbl2, lbl2, atom, ldat, atom, 2)
      iblk(1) = iblk(1) + nblk

      if(modele .eq. 1)then
c
c     Skip B, Gamma, Alpha:
c
      if(ithoma.gt.1)iblk(1) = iblk(1) + ceilq(nwn*lbl2*lbl1,iobll)
      iblk(1) = iblk(1) + ceilq(nwn*lbl2*lbl3,iobll)
      if(ithoma.lt.nthoma)iblk(1) = iblk(1)
     ;+ ceilq((nwn*lbl3+1)*lbl3,iobll)
      end if
c
c     Solve lower triangular system, store its solution:
c
      call cgesl2(atom, ldat, lbl2, ipvt, bel, ldbel, nrhs)
      call packma(lbl2, nrhs, bel, ldbel, ctom, 1)
      call aswr(2, nwn*lbl2*nrhs, ctom, 0, nblk)
      iblk(2) = iblk(2) + nblk
c
c     Initialize r.h.s. for next user block assembly:
c
      if(modele.eq.0 .and. ithoma.lt.nthoma)then
      do 31 j = 1, nrhs
      do 32 i = 1, lbl3
      bel(i,j) = bel(i+lbl2,j)
  32  continue
      do 31 i = lbl3+1, lbl3+lbl2
      bel(i,j) = czero
  31  continue
      else
c      do 33 i = 1, nrhs*ldbel
c      bel(i,1) = czero
c  33  continue
      do 33 j = 1, nrhs
      do 33 i = 1, ldbel
      bel(i,j) = czero
  33  continue
      end if
      if(ithoma .lt. nthoma)return
      call aswa(2)
c
              end if
c             ******
c
c     Now solve upper triangular system, if any:
c
      iblk(1) = iprdat(iprobl(1),1) + tstore
      iblk(2) = iprdat(iprobl(2),2) + blorhs
     ;- ceilq(nwn*nrhs*lbl2,iobll)
      call packma( lbl2, nrhs, bel, ldbel, ctom, 2 )
c
      ncw = lbl2*nrhs
      call aswr(2, nwn*ncw, ctom, 0, nblk)
      call aswa(2)
      if(modfac.eq.0)then
c     Set disk block number on start of Yn:
      iblk(1) = iblk(1) - nblk
      end if
      iblk(2) = iblk(2) - ceilq(nwn*lbl1*nrhs,iobll)
c
      if(modfac .eq. 1)then
c     Set index on last block ALPHA:
      ncw = lbl2**2
      nblk = ceilq(nwn*ncw+lbl2,iobll)
      iblk(1) = iblk(1) - nblk
          if(modele .eq. 1)then
c         also skip block b:
          ncw = lbl2*lbl1
          nblk = ceilq(nwn*ncw,iobll)
          iblk(1) = iblk(1) - nblk
          end if
      end if
c
c     Backward substitution loop:
c
      do 35 ithoma = nthoma-1, 1, -1
      lbl1 = 0
      if(ithoma .gt. 1)lbl1 = lblock(ithoma-1)
      lbl2 = lblock(ithoma)
      lbl3 = lblock(ithoma+1)
      lbl4 = lbl3 + nrhs
c
c     Set index at start of block gamma:
c
      if(modfac .eq. 0)then
      ncw = lbl2*lbl4
      else
      ncw = lbl2*lbl3
      end if
      nblk = ceilq(nwn*ncw,iobll)
      iblk(1) = iblk(1) - nblk
c
c     Read and unpack block GAMMA:
c
      if(ithoma.lt.nthoma-1) call aswa(2)
      call asrd(1, nwn*ncw, ctom, 0, nblk)
      call aswa(1)
          if(modfac.eq.1)then
          call asrd(2, nwn*lbl2*nrhs, ctom(1+lbl2*lbl3,1),
     ;              0, nblk)
          call aswa(2)
          end if
      call packma( lbl2, lbl4, ctom, ldct, ctom, 2 )

      if(modfac .eq. 1)then
c
c     Skip blocks B and ALPHA:
c
          if(modele .eq. 0)then
          ncw = lbl2 * lbl3
          else
          ncw = lbl1 * lbl2
          end if
      iblk(1) = iblk(1) - ceilq(nwn*ncw,iobll)
     ;                  - ceilq((nwn*lbl2+1)*lbl2,iobll)
      end if
c
c     Modify r.h.s.:
c
      do 36 k = 1, lbl3
      do 36 j = 1, nrhs
      do 36 i = 1, lbl2
      ctom(i,j+lbl3) = ctom(i,j+lbl3) - ctom(i,k) * bel(k,j)
  36  continue

      do 37 j = 1, nrhs
      do 37 i = 1, lbl2
      bel(i,j) = ctom(i,j+lbl3)
  37  continue
c
      ncw = lbl2*nrhs
      call packma(lbl2, nrhs, ctom(1,lbl3+1), ldct, ctom, 1)
      if(ithoma.eq.1)then
      iqueue = 1
      else
      iqueue = 0
      end if
      call aswr(2, nwn*ncw, ctom, iqueue, nblk)
      iblk(2) = iblk(2) - ceilq(nwn*lbl1*nrhs,iobll)
c
  35  continue
c
c     Successful solution: update indices
c
      nprobl(2) = nprobl(2) + 1
      iblk(2) = iprdef(nprobl(2),2)
      call aswr(2, nthoma+5, wk2zzz(1,2), 1, nblk)
      wk1zzz(1,2) = nprobl(2)
      wk1zzz(2,2) = ceilq(2*(nprobl(2)+2),iobll)
      wk1zzz(2*nprobl(2)+3,2) = iprdat(nprobl(2),2)
      iprdef(nprobl(2)+1,2) = iprdat(nprobl(2),2) + blorhs
      wk1zzz(2*nprobl(2)+4,2) = iprdef(nprobl(2)+1,2)
      iblk(2) = 1
      call aswr(2, 2*nprobl(2)+4, wk1zzz(1,2), 0, nblk)
      call aswa(2)
c
c
c            if(soldoc)then
c            end if

      ithoma = 0
      iblk(1) = iprdat(iprobl(1),1)
      iblk(2) = iprdat(iprobl(2),2)

      return
      end
c
      subroutine fbtclo(char)

      implicit none
c
c     Closes f.b.t. solver:
C-----------------------------------------------------------------------
C   CHAR: character*3 (input)
C   When CHAR='MAT', only closes the global matrix I/O unit.
C   When CHAR='ELM', only closes the element matrix I/O unit
C   (case solver called with stomat=.true.).
C   Otherwise, closes all I/O units.
C-----------------------------------------------------------------------
      include 'pardim.copy'
      include 'comzzz2.copy'
      include 'comusr.f'
c
      character*3 char
c
c     Closes global matrix I/O unit:
      if(char.eq.'MAT' .or. char.eq.'mat')then
      call ascl(1)
      else if(char.eq.'ELM' .or. char.eq.'elm')then
      if(nuni.eq.3)call ascl(3)
      else
      call ascl(1)
      call ascl(2)
      if(nuni.eq.3)call ascl(3)

          if(soldoc)then
      write(nofile,*)' '
      write(nofile,*)'*****************************************'
              if(modfac.eq.1)then
                  if(nprobl(1).eq.1)then
                  write(nofile,*)'Fortran unit ',UNIT(1),
     ;            ' contains one factorization'
                  else
                  write(nofile,*)'Fortran unit ',UNIT(1),
     ;            ' contains ',NPROBL(1),' factorizations'
                  end if
              end if
                  if(nprobl(2).eq.1)then
                  write(nofile,*)'Fortran unit ',UNIT(2),
     ;            ' contains one group of solution vectors'
                  else
                  write(nofile,*)'Fortran unit ',UNIT(2),
     ;            ' contains ',NPROBL(2),' groups of solution vectors'
                  end if
                    if(nuni.eq.3)then
                  if(nprobl(3).eq.1)then
                  write(nofile,*)'Fortran unit ',UNIT(3),
     ;            ' contains one set of element matrices and element r.h.s.'
                  else
                  write(nofile,*)'Fortran unit ',UNIT(3),' contains ',NPROBL(3)
     ;            ,' groups of element matrices and element r.h.s.'
                  end if
                    end if
      write(nofile,*)'Total number of writes: ',IWRI(1)+IWRI(2)
      write(nofile,*)'Total number of reads:  ',IREA(1)+IREA(2)
      write(nofile,*)'Total number of words written: ',
     ;iobll*(nblwr(1)+nblwr(2))
      write(nofile,*)'Total number of words read:    ',
     ;iobll*(nblrd(1)+nblrd(2))
        if(nuni.eq.3)then
      write(nofile,*)' '
      write(nofile,*)'Finite element matrices and r.h.s. stored to disk:'
      write(nofile,*)'Total number of writes: ',IWRI(3)
      write(nofile,*)'Total number of reads:  ',IREA(3)
      write(nofile,*)'Total number of words written: ',IOBLL*NBLWR(3)
      write(nofile,*)'Total number of words read:    ',IOBLL*NBLRD(3)
        end if
      write(nofile,*)' '
      write(nofile,*)'Exit frontal block-tridiagonal solver'
      write(nofile,*)'*************************************'
          end if
      end if

      return
      end
c
      subroutine rdsol(ipro,iblo,x,ldx,lbl,nrhss)

      implicit none
c
c     Reads block # IBLO of solution to problem # IPRO on unit UNIT(2)
C     Stores it into complex matrix X
C     LDX is the leading dimension of X, as declared in calling program
C     LBL is the number of rows output in X
C     NRHSS is the number of columns output in X
C
      include 'pardim.copy'
      include 'comzzz2.copy'
c
      integer i, ceilq, ipro, iblo, ldx, lbl, lbl2, ncw, nrhss, nblk
      complex*16 x(ldx,*)
      external ceilq
c
      if(ipro.le.0 .or. ipro.gt.nprobl(2))then
      write(nofile,*)'Error in RDSOL: argument 1 should be between 1'
     ;,' and ',NPROBL(2)
      stop
      end if
c
      iblk(2) = iprdef(ipro,2)
      call asrd(2, iobll, wk2zzz(1,2), 0, nblk)
      iblk(2) = iblk(2) + nblk
      call aswa(2)
      nthoma = wk2zzz(1,2)
c
      if(iblo.le.0 .or. iblo.gt.nthoma)then
      write(nofile,*)'Error in RDSOL: argument 2 should be between 1'
     ;,' and ',NTHOMA
      stop
      end if
c
          if(wk2zzz(2,2).gt.1)then
          nblk = wk2zzz(2,2) - 1
          call asrd(2, nblk*iobll, wk2zzz(iobll+1,2), 0,
     ;              nblk)
          iblk(2) = iblk(2) + nblk
          call aswa(2)
          end if
      nrhs = wk2zzz(5,2)
      do 50 ithoma = 1, nthoma
      lblock(ithoma) = wk2zzz(ithoma+5,2)
  50  continue
c
      iblk(2) = iprdat(ipro,2)
      if(iblo .gt. 1)then
      do 1 i = 1, iblo-1
      iblk(2) = iblk(2) + ceilq(nwn*nrhs*lblock(i),iobll)
  1   continue
      end if
      lbl2 = lblock(iblo)
      ncw = lbl2*nrhs
      call asrd(2, nwn*ncw, x, 0, nblk)
      call aswa(2)
      call packma(lbl2,nrhs,x,ldx,x,2)
      lbl = lbl2
      nrhss = nrhs
      return
      end
c
      subroutine wtsol(ipro,iblo,x,ldx,lbl,nrhss)

      implicit none

c     3/12/92 EXTENSION FOR CYRANO
C
C     Writes complex matrix X
C     into block # IBLO of solution to problem # IPRO on unit UNIT(2)
C     LDX is the leading dimension of X, as declared in calling program
c     lbl,nrhss dummy, x assumed the correct size!
c
      include 'pardim.copy'
      include 'comzzz2.copy'
c
      integer i, ceilq, ipro, iblo, ldx, lbl, lbl2, ncw, nrhss, nblk
      complex*16 x(ldx,*)
      external ceilq
c
      if(ipro.le.0 .or. ipro.gt.nprobl(2))then
      write(nofile,*)'Error in WTSOL: argument 1 should be between 1'
     ;,' and ',nprobl(2)
      stop
      end if
c
      iblk(2) = iprdef(ipro,2)
      call asrd(2, iobll, wk2zzz(1,2), 0, nblk)
      iblk(2) = iblk(2) + nblk
      call aswa(2)
      nthoma = wk2zzz(1,2)
c
      if(iblo.le.0 .or. iblo.gt.nthoma)then
      write(nofile,*)'Error in WTSOL: argument 2 should be between 1'
     ;,' and ',nthoma
      stop
      end if
c
          if(wk2zzz(2,2).gt.1)then
          nblk = wk2zzz(2,2) - 1
          call asrd(2, nblk*iobll, wk2zzz(iobll+1,2), 0,
     ;              nblk)
          iblk(2) = iblk(2) + nblk
          call aswa(2)
          end if
      nrhs = wk2zzz(5,2)
      do 50 ithoma = 1, nthoma
      lblock(ithoma) = wk2zzz(ithoma+5,2)
  50  continue
c
      iblk(2) = iprdat(ipro,2)
      if(iblo .gt. 1)then
      do 1 i = 1, iblo-1
      iblk(2) = iblk(2) + ceilq(nwn*nrhs*lblock(i),iobll)
  1   continue
      end if
      lbl2 = lblock(iblo)
      ncw = lbl2*nrhs
      call packma(lbl2,nrhs,x,ldx,x,1)
      call aswr(2, nwn*ncw, x, 0, nblk)
      call aswa(2)
c     lbl = lbl2
c     nrhss = nrhs
      return
      end
c
      subroutine rdelt(ipro, iel, lblg, lbld, nrhss)

      implicit none
c
c     Only use when solver has been requested to store element matrices
C     (stomat=.true. in call to FBTINI)
C     Reads 'generalized' element block # IEL and associated r.h.s. for 
C     problem # IPRO on unit UNIT(3)
C     Stores them into solver complex matrices AEL and BEL, respectively.
C     LBLG and LBLD are the block sizes associated with left and right part of
C     the element. NRHSS is the number of right hand sides.
C     Thus, on output AEL contains a square matrix of order LBLG+LBLD, and
C                     BEL a rectangular (LBLG+LBLD) by NRHSS matrix.
C     N.B.: 'generalized' element means a pair of successive block rows in the
C     global matrix. This includes boundary blocks.
C     Thus numbering is not the same as for finite elements on the mesh.

      include 'pardim.copy'
      include 'comzzz2.copy'
      include 'comusr.f'
c
      integer i, ceilq, ipro, iel, lblg, lbld, lbl2, lbl3, lbl, ncw, nrhss, nblk
      external ceilq
c
      if(nuni .ne. 3)then
      write(nofile,*)
     ;'Error calling RDELT: solver not initialized to read element matrices'
      stop
      end if
      if(ipro.le.0 .or. ipro.gt.nprobl(3))then
      write(nofile,*)'Error in RDELT: argument 1 should be between 1'
     ;,' and ',nprobl(3)
      stop
      end if
c
      iblk(3) = iprdef(ipro,3)
      call asrd(3, iobll, wk2zzz(1,3), 0, nblk)
      iblk(3) = iblk(3) + nblk
      call aswa(3)
      nthoma = wk2zzz(1,3)
c
      if(iel.le.0 .or. iel.ge.nthoma)then
      write(nofile,*)'Error in RDELT: argument 2 should be between 1'
     ;,' and ',NTHOMA-1
      stop
      end if
c
        if(wk2zzz(2,3).gt.1)then
        nblk = wk2zzz(2,3) - 1
        call asrd(3, nblk*iobll, wk2zzz(iobll+1,3), 0, nblk)
        iblk(3) = iblk(3) + nblk
        call aswa(3)
        end if
      nrhs = wk2zzz(5,3)
      do 50 ithoma = 1, nthoma
      lblock(ithoma) = wk2zzz(ithoma+5,3)
  50  continue
c
      iblk(3) = iprdat(ipro,3)

      if(iel .gt. 1)then
      do 1 i = 1, iel-1
      lbl = lblock(i) + lblock(i+1)
      iblk(3) = iblk(3) + ceilq(nwn*lbl**2,iobll)
     ;                  + ceilq(nwn*nrhs*lbl,iobll)
  1   continue
      end if
      lbl2 = lblock(iel)
      lbl3 = lblock(iel+1)
      lbl = lbl2 + lbl3
      ncw = lbl ** 2
      call asrd(3, nwn*ncw, ael, 1, nblk)
      iblk(3) = iblk(3) + nblk
      ncw = lbl * nrhs
      call asrd(3, nwn*ncw, bel, 0, nblk)
      iblk(3) = iblk(3) + nblk
      call aswa(3)
      call packma(lbl, lbl, ael, ldael, ael, 2)
      call packma(lbl, nrhs, bel, ldbel, bel, 2)

      lblg = lbl2
      lbld = lbl3
      nrhss = nrhs
      return
      end
c
      subroutine cgeco2(a,lda,n,ipvt,scin,sc1,sc2,rcond,info,z)

      implicit none

      integer lda, n, ipvt(*), scin, info
      complex*16 a(lda,*), z(*)
      double precision rcond, sc1(*), sc2(*)
c
C     CGECO2 FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C     SCALED ROW PIVOTING IS PERFORMED.
C
C     IF  RCOND  IS NOT NEEDED, CGEFA2 IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW CGECO2 BY CGESL2.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        SCIN    INTEGER
C                A SWITCH FOR THE INVERSE SCALING FACTORS SC1, SC2:
C                SCIN=0: SC1 IS WORKSPACE, SC2 IS A DUMMY ARGUMENT.
C                        CGECO2 COMPUTES THE MAXIMUM NORMS OF THE ROWS
C                        OF A AND STORES THEM IN SC1.
C                SCIN=1: SC1 IS INPUT, SC2 IS A DUMMY ARGUMENT.
C                OTHERWISE: SC1 AND SC2 ARE INPUT .
C
C        SC1     REAL(N)
C                WORKSPACE WHEN SCIN=0.
C                INPUT WHEN SCIN IS NONZERO. SC1 MUST THEN
C                CONTAIN INVERSE SCALING FACTORS FOR THE ROWS OF A.
C                THESE FACTORS WILL BE USED IN THE
C                SELECTION OF PIVOTS AND THE CONDITION NUMBER ESTIMATION
C
C        SC2     REAL(N)
C                DUMMY ARGUMENT WHEN SCIN=0 OR SCIN=1.
C                INPUT OTHERWISE. SC2 MUST THEN CONTAIN INVERSE SCALING
C                FACTORS FOR THE COLUMNS OF A.
C                THESE FACTORS WILL BE USED IN THE CONDITION NUMBER
C                ESTIMATION.
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       COMPLEX(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     MODIFIED VERSION OF LINPACK CGECO. THIS VERSION DATED 1/10/92 .
C     PURPOSE: SCALED ROW PIVOTING. THE SCALING FACTORS MAY BE INPUT
C              OR OUTPUT.
C
C     SUBROUTINES AND FUNCTIONS
C
C     CGEFA2
C     BLAS ZAXPY,ZDOTC,ZDSCAL,DZASUM
C     FORTRAN DABS,DIMAG,DMAX1,DCMPLX,DCONJG,DREAL
C
C     INTERNAL VARIABLES
C
      COMPLEX*16 ZDOTC, EK, T, WK, WKM, AKKCI, ZERO, ONE
      DOUBLE PRECISION ANORM, S, SM, YNORM, RW1, RW2, SCA
      INTEGER I, J, K, KB, KP1, L, izamax
      EXTERNAL izamax
C
      COMPLEX*16 ZDUM, ZDUM1, ZDUM2, CSIGN1
      DOUBLE PRECISION CABS1, DZASUM
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C     For real version, replace preceding line by following:
C     CSIGN1(ZDUM1,ZDUM2) = DSIGN(ZDUM1,ZDUM2)
      DATA ZERO/(0.0D0,0.0D0)/,ONE/(1.0D0,0.0D0)/
C
C     COMPUTE ROW MAXIMUM NORMS WHEN SCIN=0:
C
      IF(SCIN.EQ.0)THEN
      DO 1 I = 1, N
      J = izamax(N, A(I,1), LDA)
      SC1(I) = CABS1( A(I,J) )
c       Add 6/8/98:
        if(sc1(i) .eq. 0.d0)then
        write(6,*)'CGECO2: row #',i,' is zero; solver stopped'
        stop
        end if
  1   CONTINUE
      END IF
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 31 I = 1, N
      Z(I) = DCMPLX(1.D0/SC1(I),0.0D0)
  31  CONTINUE
      DO 30 J = 1, N
      SCA = 0.
      DO 200 I = 1, N
      SCA = SCA + CABS1(A(I,J))*DREAL(Z(I))
  200 CONTINUE
         IF(SCIN.EQ.0 .OR. SCIN.EQ.1)THEN
         ANORM = DMAX1(ANORM,SCA)
         ELSE
         ANORM = DMAX1(ANORM,SCA/SC2(J))
         END IF
   30 CONTINUE
C
C     FIRST, FACTOR A WITH CGEFA2 :
C
      CALL CGEFA2(A,LDA,N,IPVT,SC1,1,INFO)
c Add 6/8/98:
      if(INFO .ne. 0)write(6,*)'CGEFA2 called by CGECO2 returned INFO =', info
c
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
C     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE CTRANS(U)*W = E
C
      EK = ONE
      DO 20 I = 1, N
      Z(I) = ZERO
  20  CONTINUE

      DO 100 K = 1, N
         RW1 = CABS1(A(K,K))
         IF (CABS1(Z(K)) .NE. 0.0D0)THEN
      IF(SCIN.EQ.0 .OR. SCIN.EQ.1)THEN
         EK = CSIGN1(EK,-Z(K))
      ELSE
         EK = CSIGN1(EK,-Z(K)) * SC2(K)
      END IF
         END IF
         RW2 = CABS1(EK-Z(K))
         IF (RW2 .GT. RW1)THEN
            S = RW1/RW2
            CALL ZDSCAL(N,S,Z,1)
            EK = DCMPLX(S,0.0D0)*EK
         END IF
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(A(K,K)) .EQ. 0.0D0) THEN
            WK = ONE
            WKM = ONE
         ELSE
            AKKCI = ONE/DCONJG(A(K,K))
            WK = WK*AKKCI
            WKM = WKM*AKKCI
         END IF
         KP1 = K + 1
         IF (KP1 .LE. N) THEN
            DO 60 J = KP1, N
               SM = SM + CABS1(Z(J)+WKM*DCONJG(A(K,J)))
               Z(J) = Z(J) + WK*DCONJG(A(K,J))
               S = S + CABS1(Z(J))
   60       CONTINUE
            IF (S .LT. SM) THEN
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*DCONJG(A(K,J))
   70          CONTINUE
            END IF
         END IF
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DZASUM(N,Z,1)
      CALL ZDSCAL(N,S,Z,1)
C
C     SOLVE CTRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + ZDOTC(N-K,A(K+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .GT. 1.0D0) THEN
            S = 1.0D0/CABS1(Z(K))
            CALL ZDSCAL(N,S,Z,1)
         END IF
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
C     SCALING CONSISTENT WITH SCALED PIVOTING:
      DO 210 I = 1, N
      Z(I) = Z(I)*SC1(I)**2
  210 CONTINUE
      S = 1.0D0/DZASUM(N,Z,1)
      CALL ZDSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL ZAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .GT. 1.0D0) THEN
            S = 1.0D0/CABS1(Z(K))
            CALL ZDSCAL(N,S,Z,1)
            YNORM = S*YNORM
         END IF
  140 CONTINUE
      S = 1.0D0/DZASUM(N,Z,1)
      CALL ZDSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         RW1 = CABS1(Z(K))
         RW2 = CABS1(A(K,K))
         IF (RW1 .GT. RW2) THEN
            S = RW2/RW1
            CALL ZDSCAL(N,S,Z,1)
            YNORM = S*YNORM
         END IF
         IF (RW2 .EQ. 0.0D0) THEN
         Z(K) = ONE
         ELSE
         Z(K) = Z(K)/A(K,K)
         END IF
         T = -Z(K)
         CALL ZAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
      IF(SCIN.NE.0 .AND. SCIN.NE.1)THEN
      DO 170 I = 1, N
      Z(I) = SC2(I) * Z(I)
  170 CONTINUE
      END IF
C     MAKE ZNORM = 1.0
      S = 1.0D0/DZASUM(N,Z,1)
      CALL ZDSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .EQ. 0.0D0) THEN
      RCOND = 0.0D0
      ELSE
      RCOND = YNORM/ANORM
      END IF
      RETURN
      END
C
      SUBROUTINE CGEFA2(A,LDA,N,IPVT,SCALE,SCIN,INFO)

      IMPLICIT NONE

      INTEGER LDA,N,IPVT(*),INFO ,SCIN
      DOUBLE PRECISION    SCALE(*)
      COMPLEX*16 A(LDA,*)
C
C     CGEFA2 FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION.
C
C     CGEFA2 IS USUALLY CALLED BY CGECO2, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        SCIN    INTEGER
C                A SWITCH FOR THE ROW INVERSE SCALING FACTOR SCALE:
C                SCIN=0: SCALE IS OUTPUT
C                SCIN NONZERO: SCALE IS AN INPUT .
C
C        SCALE   REAL(N)
C                INPUT WHEN SCIN IS NONZERO. SCALE MUST THEN CONTAIN
C                THE INVERSE SCALING FACTOR OF EACH EQUATION TO BE USED
C                IN SEARCH OF PIVOTS.
C
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        SCALE   REAL(N)
C                WHEN SCIN=0, THE MAXIMUM NORM OF EACH ROW OF A
C                USING THE CABS1 FUNCTION FOR ABSOLUTE VALUES .
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT CGESL OR CGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN CGECO2 FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     MODIFIED VERSION OF LINPACK CGEFA. THIS VERSION DATED 1/10/92 .
C     PURPOSE: SCALED ROW PIVOTING. THE SCALING FACTORS MAY BE INPUT
C              OR OUTPUT.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS ZAXPY,CSCAL,izamax
C     FORTRAN ABS,AIMAG,REAL
C
C     INTERNAL VARIABLES
C
      INTEGER izamax, I, J, K, KP1, L, NM1
      DOUBLE PRECISION THR, THR2
      COMPLEX*16 T, ONE
      EXTERNAL izamax
C
      COMPLEX*16 ZDUM
      DOUBLE PRECISION CABS1
      DATA ONE/(1.0D0,0.0D0)/
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
C
C     ROW MAXIMUM NORMS ARE COMPUTED WHEN SCIN=0:
C
      IF(SCIN.EQ.0)THEN
      DO 1 I = 1, N
      J = IZAMAX(N, A(I,1), LDA)
      SCALE(I) = CABS1( A(I,J) )
  1   CONTINUE
      END IF
C
C     GAUSSIAN ELIMINATION WITH SCALED PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
      L = K
      THR = CABS1(A(K,K))
      DO 2 I = KP1, N
      THR2 = CABS1(A(I,K))
      IF(THR2*SCALE(L).GT.THR*SCALE(I))THEN
      L = I
      THR = THR2
      END IF
  2   CONTINUE
      IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (THR .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .NE. K) THEN
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
            END IF
C
C           COMPUTE MULTIPLIERS
C
            T = -ONE/A(K,K)
            CALL zscal(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .NE. K) THEN
                  A(L,J) = A(K,J)
                  A(K,J) = T
               END IF
               CALL ZAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF(CABS1(A(N,N)) .EQ. 0.0D0)INFO = N
      RETURN
      END
C
      SUBROUTINE CGESL2(A,LDA,N,IPVT,B,LDB,NRHS)

      IMPLICIT NONE

C
      INTEGER LDA, LDB, N, NRHS, IPVT(*)
      COMPLEX*16 A(LDA,*), B(LDB,*)
C
C     CGESL2 SOLVES THE COMPLEX SYSTEM
C     A * X = B
C     WHERE X AND B ARE N * NRHS MATRICES,
C     USING THE FACTORIZATION OF A PREVIOUSLY STORED IN A BY CGECO2
C     OR CGEFA2.
C
C     ON ENTRY
C
C        A       COMPLEX(LDA, N)
C                THE OUTPUT FROM CGECO2 OR CGEFA2.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM CGECO2 OR CGEFA2.
C
C        B       COMPLEX(N,NRHS)
C                THE RIGHT HAND SIDE MATRIX
C
C        LDB     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY B .
C
C        NRHS    INTEGER
C                THE NUMBER OF RIGHT HAND SIDE VECTORS IN B .
C                THE RIGHT HAND SIDES ARE STORED IN COLUMNS OF B.
C
C     ON RETURN
C
C        B       THE NRHS SOLUTION VECTORS X, REPLACING THE
C                CORRESPONDING INPUT RIGHT HAND SIDES.
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF CGECO2 HAS SET RCOND .GT. 0.0
C        OR CGEFA2 HAS SET INFO .EQ. 0 .
C
C
C     MODIFIED VERSION OF LINPACK CGESL. THIS VERSION DATED 1/10/92 .
C     PURPOSE: VECTORIZATION OVER LARGE NUMBER OF RIGHT HAND SIDE
C              VECTORS (NRHS).
C
C     INTERNAL VARIABLES
C
      INTEGER I, J, K, L, NM1
      COMPLEX*16 ZERO, ONE, TEMP
      DATA ZERO/(0.0D0,0.0D0)/, ONE/(1.0D0,0.0D0)/
C
      NM1 = N - 1
C
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .GT. 0) THEN
         DO 10 K = 1, NM1
            L = IPVT(K)
            IF(L .NE. K)THEN
            DO 20 J = 1, NRHS
            TEMP = B(K,J)
            B(K,J) = B(L,J)
            B(L,J) = TEMP
   20       CONTINUE
            END IF
            DO 10 I = K+1, N
            TEMP = A(I,K)
            IF(TEMP.NE.ZERO)THEN
            DO 30 J = 1, NRHS
            B(I,J) = B(I,J) + TEMP * B(K,J)
   30       CONTINUE
            END IF
   10    CONTINUE
         END IF
C
C        NOW SOLVE  U*X = Y
C
         DO 40 K = N, 1, -1
            TEMP  = ONE/A(K,K)
            DO 50 J = 1, NRHS
            B(K,J) = TEMP * B(K,J)
   50       CONTINUE
         DO 40 I = 1, K-1
         TEMP = A(I,K)
         IF(TEMP.NE.ZERO)THEN
         DO 60 J = 1, NRHS
         B(I,J) = B(I,J) - TEMP * B(K,J)
   60    CONTINUE
         END IF
   40    CONTINUE
C
          RETURN
          END
C

      SUBROUTINE PACKMA(NRO, NCOL, A, LDA, VEC, IPATH)

      IMPLICIT NONE

C     PACKS COMPLEX MATRIX A INTO COMPLEX VECTOR VEC   (IPATH=1)
C  OR UNPACKS COMPLEX VECTOR VEC INTO COMPLEX MATRIX A (IPATH=2)
C     A AND VEC CAN SHARE THE SAME STORAGE LOCATION WITH A(1,1)=VEC(1)

      INTEGER NRO, NCOL, LDA, IPATH
      COMPLEX*16 A(LDA,NCOL), VEC(NRO*NCOL)
      
      INTEGER I, J, IVEC

      IF(IPATH .EQ. 1)THEN
      IVEC = 0
      DO 1 J = 1, NCOL
      DO 1 I = 1, NRO
      IVEC = IVEC + 1
      VEC(IVEC) = A(I,J)
   1  CONTINUE

      ELSE IF(IPATH.EQ.2)THEN
      IVEC = NRO * NCOL
      DO 2 J = NCOL, 1, -1
      DO 2 I = NRO, 1, -1
      A(I,J) = VEC(IVEC)
      IVEC = IVEC - 1
   2  CONTINUE
      END IF

      RETURN
      END
C
      INTEGER FUNCTION CEILQ(I, J)

      IMPLICIT NONE

      INTEGER I, J, MOD
      INTRINSIC  MOD
      
      CEILQ = I / J
      IF(MOD(I,J).NE.0)CEILQ = CEILQ + 1
      RETURN
      END
c
c     Calls to input-output routines; case of Cray asynchronous queued
c     I/O:
c
      SUBROUTINE ASOP(IUNI)

      IMPLICIT NONE
C
C     Opens unit UNIT(IUNI), IUNI= 1, 2 or 3.
C
c     Author: Philippe Lamalle
c     1/10/92 Version
c     Modified 24/10/2000 to allow storage of element matrices

      INCLUDE 'PARDIM.COPY'
      INCLUDE 'COMZZZ2.COPY'
C
      INTEGER IUNI, ISTATU
C
      IF(IUNI.GE.1 .AND. IUNI.LE.3)THEN
C     Next line is for non-Cray version: check no effect on Cray!
      AQP(1,IUNI) = UNIT(IUNI)
C     Cray AQIO routine:
      CALL AQOPEN(AQP(1,IUNI),NAQP,UNIT(IUNI),ISTATU)
C     In case of anomalous exit:
      IF(ISTATU.NE.0)
     ;WRITE(NOFILE,*)'AQOPEN: UNIT #',UNIT(IUNI),' ISTATU=',ISTATU
      ELSE
      WRITE(NOFILE,*)'No such unit id.', IUNI
      END IF
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE ASCL(IUNI)

      IMPLICIT NONE

C
C     Closes unit UNIT(IUNI), IUNI = 1, 2 or 3
C
c     Author: Philippe Lamalle
c     1/10/92 Version
c     Modified 24/10/2000 to allow storage of element matrices

      INCLUDE 'PARDIM.COPY'
      INCLUDE 'COMZZZ2.COPY'
C
      INTEGER IUNI, ISTATU
C
      IF(IUNI.GE.1 .AND. IUNI.LE.3)THEN
C     Cray AQIO routine:
      CALL AQCLOSE(AQP(1,IUNI),ISTATU)
      ELSE
      WRITE(NOFILE,*)'No such unit id.', IUNI
      END IF
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE ASRD(IUNI, NWORD, VEC, IQUEUE, NBLK)

      IMPLICIT NONE

C
C     Reads NWORD words on unit UNIT(IUNI), starting at current record
C     IBLK(IUNI), and stores them into real vector VEC.
C
C     IQUEUE=0: asks to process asynchronous request on UNIT(IUNI)
C     IQUEUE=1: queues the request without initiating input.
C
C     Subroutine output: NBLK, the number of blocks of length IOBLL read
C
C     N.B.: as blocked input is used, VEC must be large enough to
C     accomodate NBLK blocks of IOBLL words, that is the smallest number
C     of I/O blocks containing NWORD words.
C
c     Author: Philippe Lamalle
c     1/10/92 Version
c     Modified 24/10/2000 to allow storage of element matrices

      INCLUDE 'PARDIM.COPY'
      INCLUDE 'COMZZZ2.COPY'
C
      INTEGER IUNI, NWORD, IQUEUE, CEILQ, NBLK, ISTATU
      DOUBLE PRECISION VEC(*)
C
      IF(IUNI.GE.1 .AND. IUNI.LE.3)THEN
      NBLK = CEILQ(NWORD,IOBLL)
      IREA(IUNI) = IREA(IUNI) + 1
C     Cray AQIO routine:
      CALL AQREAD(AQP(1,IUNI), VEC, IBLK(IUNI), NBLK, IREA(IUNI)
     ;           ,IQUEUE, ISTATU)
      IF(ISTATU.NE.0)
     ;WRITE(NOFILE,*)'AQREAD: UNIT #',UNIT(IUNI),' ISTATU=',ISTATU
      NBLRD(IUNI) = NBLRD(IUNI) + NBLK
      ELSE
      WRITE(NOFILE,*)'No such unit id.', IUNI
      END IF
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE ASWR(IUNI, NWORD, VEC, IQUEUE, NBLK)

      IMPLICIT NONE

C
C     Writes NWORD words from real vector VEC on unit IUNI,
C     starting at current record IBLK(IUNI).
C
C     IQUEUE=0: asks to process asynchronous request on IUNI
C     IQUEUE=1: queues the request without initiating output.
C
C     Subroutine output: NBLK, the number of blocks of length
C     IOBLL written.
C
C     N.B.: as blocked output is used, VEC must be large enough to
C     accomodate NBLK blocks of IOBLL words, that is the smallest number
C     of I/O blocks containing NWORD words.
C
C
c     Author: Philippe Lamalle
c     1/10/92 Version
c     Modified 24/10/2000 to allow storage of element matrices

      INCLUDE 'PARDIM.COPY'
      INCLUDE 'COMZZZ2.COPY'
C
      INTEGER IUNI, NWORD, IQUEUE, CEILQ, NBLK, ISTATU
      DOUBLE PRECISION VEC(*)

      IF(IUNI.GE.1 .AND. IUNI.LE.3)THEN
      NBLK = CEILQ(NWORD,IOBLL)
      IWRI(IUNI) = IWRI(IUNI) + 1
C     Cray AQIO routine:
      CALL AQWRITE(AQP(1,IUNI), VEC, IBLK(IUNI), NBLK, IWRI(IUNI)
     ;            ,IQUEUE, ISTATU)
      IF(ISTATU.NE.0)
     ;WRITE(NOFILE,*)'AQWRITE: UNIT #',UNIT(IUNI),' ISTATU=',ISTATU
      NBLWR(IUNI) = NBLWR(IUNI) + NBLK
      ELSE
      WRITE(NOFILE,*)'No such unit id.', IUNI
      END IF

      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE ASWA(IUNI)

      IMPLICIT NONE

C
C     Waits for completion of asynchronous I/O on unit UNIT(IUNI)
C     IUNI = 1, 2 or 3
C
c     Author: Philippe Lamalle
c     1/10/92 Version
c     Modified 24/10/2000 to allow storage of element matrices

      INCLUDE 'PARDIM.COPY'
      INCLUDE 'COMZZZ2.COPY'
C
      INTEGER IUNI, ISTATU
C
      IF(IUNI.GE.1 .AND. IUNI.LE.3)THEN
C     Cray AQIO routine:
      CALL AQWAIT(AQP(1,IUNI), ISTATU)
      IF(ISTATU.NE.0)
     ;WRITE(NOFILE,*)'AQWAIT: UNIT #',UNIT(IUNI),' ISTATU=',ISTATU
      ELSE
      WRITE(NOFILE,*)'No such unit id.', IUNI
      END IF
C
      RETURN
      END
