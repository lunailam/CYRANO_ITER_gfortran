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
     ;    ((maxbll**2*nwn)/iobll+(maxbll*maxrhs*nwn)/iobll + 2)
     ;  , ldael=(2-modele)*maxbll,lcael=(modele+2)*maxbll
     ;  , ldbel=ldael
     ;  , ldbe1=((ldbel*maxrhs*nwn)/iobll+1)*iobll
     ;  , luszzz = ldat1 + ldbt1 + ldct1 + ldbe1
     ;            +nwn*( maxbll + ldael*lcael )
     ;)
c Previous, giving trouble:
c     ;  , LDAT1=IOBLL*((LDAT**2+MAXBLL*MODFAC)*NWN/IOBLL+1)
c     ;  , LDBT1=IOBLL*(LDBT**2*NWN/IOBLL+1)
c     ;  , LDCT1=IOBLL*(MAXBLL**2*NWN/IOBLL+(MAXBLL*MAXRHS)*NWN/IOBLL +2)
c     ;  , LDAEL=(2-MODELE)*MAXBLL,LCAEL=(MODELE+2)*MAXBLL
c     ;  , LDBEL=LDAEL
c     ;  , LDBE1=((LDBEL*MAXRHS)*NWN/IOBLL+1)*IOBLL
C*******************************************************************
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
 
      integer nbmax, bllen
      parameter(nbmax=(maxcou+1)*(2*maxpom-maxcou)+2*(maxpom-1) )
c      parameter( bllen = 6 * nbmax )
c     why was this 6?
      parameter(bllen = 2 * nbmax + 64)

c      double precision
c     ;  tact(2*maxpom-1,0:maxcou,-1:1)
c     ;, treact(2*maxpom-1,0:maxcou,-1:1)

      complex*16 
     ;  vmat(ndof,ndof,2*maxpom**2), totop(2*ndof,2*ndof,bllen)
     ;, cbc(maxbll,2*maxbll), cbctc(2*maxbll,maxbll)
cPL30/5/04: sorting the case giplas=.F.
      complex*16 totoc(2*ndof,2*ndof,bllen)

c     Trying to save space:
      EQUIVALENCE (CBC,CTOM), (CBCTC,VMAT)
c     Warning: don't change length of VMAT before relocating CBCTC!
c     They presently have exactly the same length.
c     CBC covers CTOM and (part of) BTOM!
 
      equivalence 
     ;  (atom1,atom), (btom1,btom), (ctom1,ctom), (bel1,bel)
     ;, (ctom1,zzzzzz), (btom1,zzzzzz(ldct1+1))
     ;, (atom1,zzzzzz(ldct1+ldbt1+1))
     ;, (ael,zzzzzz(ldct1+ldbt1+ldat1+1))
     ;, (bel1,zzzzzz(ldct1+ldbt1+ldat1+nwn*ldael*lcael+1))
     ;, (wkczzz,zzzzzz(ldct1+ldbt1+ldat1+nwn*ldael*lcael+ldbe1+1))
 
c      common/ccomus/ctom1,btom1,atom1,ael,bel1,wkczzz
cPL30/5/04: sorting giplas=.F. added totoc array for curl
      common/ccomus/zzzzzz
     ;, vmat, totop, totoc 
c     ;, tact, treact
c     ;, cbc, cbctc
