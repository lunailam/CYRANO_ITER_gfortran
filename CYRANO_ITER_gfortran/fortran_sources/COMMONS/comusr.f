C*******************************************************************
C     Solver internal parameters:
C     1/10/92 Version
      integer   maxlda,ldat,ldbt,ldct,ldat1,ldbt1,ldct1,ldael
     ;         ,ldbel,lcael,ldbe1,luszzz
      parameter (
     ;    maxlda=maxbll*maxnbl
     ;  , ldat=(maxbll-1)*modfac+1
     ;  , ldbt=(maxbll-1)*(1+modele*(1+modfac-2*modele))+1
     ;  , ldct= maxbll
     ;  , ldat1=iobll*(((ldat**2+maxbll*modfac)*nwn)/iobll+1)
     ;  , ldbt1=iobll*((ldbt**2*nwn)/iobll+1)
     ;  , ldct1=iobll*
     ;    ((maxbll**2*nwn)/iobll+(maxbll*maxrhs*nwn)/iobll+2)
     ;  , ldael=(2-modele)*maxbll,lcael=(modele+2)*maxbll
     ;  , ldbel=ldael
     ;  , ldbe1=((ldbel*maxrhs*nwn)/iobll+1)*iobll
     ;  , luszzz = ldat1 + ldbt1 + ldct1 + ldbe1
     ;            +nwn*(maxbll + ldael*lcael)
     ;)
c previous, giving trouble:
c     ;  , ldat1=iobll*((ldat**2+maxbll*modfac)*nwn/iobll+1)
c     ;  , ldbt1=iobll*(ldbt**2*nwn/iobll+1)
c     ;  , ldct1=iobll*(maxbll**2*nwn/iobll+(maxbll*maxrhs)*nwn/iobll +2)
c     ;  , ldael=(2-modele)*maxbll,lcael=(modele+2)*maxbll
c     ;  , ldbel=ldael
c     ;  , ldbe1=((ldbel*maxrhs)*nwn/iobll+1)*iobll
c
c*******************************************************************
C     Frontal block-tridiagonal solver common block CCOMUS
C     Accessed by user program:
C     At each step, the user builds elementary blocks of the system
C     matrix in AEL(,), and the corresponding contributions to
C     right-hand-sides in BEL(,).
C     (The block storage mode is selected with switch MODELE
C     (=0 or 1) before compiling the solver library).
C     The other arrays are solver workspace and should not be modified.
C     The whole common block memory can of course be re-employed after
C     solver exits.

c Added ZZZZZZ and several equivalences: 
      double precision 
     ;  atom1(ldat1), btom1(ldbt1), ctom1(ldct1), bel1(ldbe1)
     ;, zzzzzz(luszzz)

      complex*16
     ;  ael(ldael,lcael), bel(ldbel,maxrhs)
     ;, atom(ldat,ldat), btom(ldbt,ldbt)
     ;, ctom(ldct,maxbll+maxrhs), wkczzz(maxbll)
 
      equivalence 
     ;  (atom1,atom), (btom1,btom), (ctom1,ctom), (bel1,bel)
     ;, (ctom1,zzzzzz), (btom1,zzzzzz(ldct1+1))
     ;, (atom1,zzzzzz(ldct1+ldbt1+1))
     ;, (ael,zzzzzz(ldct1+ldbt1+ldat1+1))
     ;, (bel1,zzzzzz(ldct1+ldbt1+ldat1+nwn*ldael*lcael+1))
     ;, (wkczzz,zzzzzz(ldct1+ldbt1+ldat1+nwn*ldael*lcael+ldbe1+1))
 
c      common/ccomus/ctom1,btom1,atom1,ael,bel1,wkczzz
      common/ccomus/zzzzzz
