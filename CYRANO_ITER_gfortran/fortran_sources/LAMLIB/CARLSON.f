c      In this file: routines for Carlson's elliptic integrals
c      rf, rd, rj, rc
c      and for Legendre's elliptic integrals ellf, elle, ellpi
c      in terms of the former
c
c      RF function: Carlson's elliptic integral of the first kind
c      "Numerical Recipes" - second edition, p.257
c
       function rf(x,y,z)
       implicit none
       real*8 rf,x,y,z,errtol,tiny,big,third,c1,c2,c3,c4
       integer maxit, nit

       parameter (maxit=50)
c      parameter (errtol=0.08)
       parameter (errtol=0.0025)
       parameter (tiny=1.5e-38,big=3.e37,third=1./3.,
     +            c1=1./24.,c2=0.1,c3=3./44.,c4=1./14.)

c      errtol=0.08 guarantees correct results in single precision (7 digits)
c      errtol=0.0025 guarantees correct results in double precision (16 digits)
c      tiny is at least 5 times machine underflow limit
c      big is at most 1/5 of the machine overflow limit
       real*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
       if(min(x,y,z).lt.0. .or. min(x+y,x+z,y+z).lt.tiny .or.
     +    max(x,y,z).gt.big)then
        write(6,*)'invalid arguments in function rf; x,y,z: ', x, y, z
        stop
       endif
       
       xt=x
       yt=y
       zt=z
       nit = 0
       
 1     continue
       sqrtx=dsqrt(xt)
       sqrty=dsqrt(yt)
       sqrtz=dsqrt(zt)
       alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
       xt=0.25*(xt+alamb)
       yt=0.25*(yt+alamb)
       zt=0.25*(zt+alamb)
       ave=third*(xt+yt+zt)
       delx=(ave-xt)/ave
       dely=(ave-yt)/ave
       delz=(ave-zt)/ave
       nit = nit + 1
       if(max(dabs(delx),dabs(dely),dabs(delz)).gt.errtol .and. nit.lt.maxit)
     ; goto 1
       if(nit.ge.maxit)write(6,*)'Carlson rf:', nit, ' iterations', x, y, z
       
       e2=delx*dely-delz**2
       e3=delx*dely*delz
       rf=(1.+(c1*e2-c2-c3*e3)*e2+c4*e3)/dsqrt(ave)
       
       return
       end
c
c      RD function
c
       function rd(x,y,z)
       implicit none
       real*8 rd,x,y,z,errtol,tiny,big,c1,c2,c3,c4,c5,c6
       integer maxit, nit

       parameter (maxit=50)
c      parameter (errtol=0.05)
       parameter (errtol=0.0015)
       parameter (tiny=1.e-25,big=4.5e21,
     +            c1=3./14.,c2=1./6.,c3=9./22.,c4=3./26.,
     +            c5=0.25*c3,c6=1.5*c4)

c      errtol=0.05 guarantees correct results in single precision (7 digits)
c      errtol=0.0015 guarantees correct results in double precision (16 digits)
c      tiny is at least twice the cube root of the machine underflow limit
c      big is at most 1/5 the cube root of the machine overflow limit
       real*8 alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx
     + , sqrty
     + , sqrtz, sum, xt, yt, zt
         if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.
     +      max(x,y,z).gt.BIG)then
         write(6,*)'invalid arguments in rd: ', x, y, z
         stop
         end if
       xt = x
       yt = y
       zt = z
       sum = 0.
       fac = 1.
       nit = 0
   1   continue
           sqrtx=dsqrt(xt)
           sqrty=dsqrt(yt)
           sqrtz=dsqrt(zt)
           alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
           sum = sum + fac / (sqrtz*(zt+alamb))
           fac = 0.25*fac
           xt=0.25*(xt+alamb)
           yt=0.25*(yt+alamb)
           zt=0.25*(zt+alamb)
           ave=0.2*(xt+yt+3.*zt)
           delx = (ave-xt)/ave
           dely = (ave-yt)/ave
           delz = (ave-zt)/ave
           nit = nit + 1
       if(max(dabs(delx),dabs(dely),dabs(delz)).gt.ERRTOL .and. nit.lt.maxit)
     + goto 1 
       if(nit.ge.maxit)write(6,*)'Carlson rd:', nit, ' iterations', x, y, z
       ea = delx*dely
       eb = delz*delz
       ec = ea - eb
       ed = ea - 6.*eb
       ee = ed+ec+ec
       rd = 3.*sum+fac*(1.+ed*(-c1+c5*ed-c6*delz*ee)
     +        +delz*(c2*ee+delz*(-c3*ec+delz*c4*ea))) / (ave*dsqrt(ave))  
       return
       end  
c
c      RJ function: Carlson's elliptic integral of the third kind
c      "Numerical Recipes" - second edition, p.257
c
       function rj(x,y,z,p)
       implicit none
       real*8 rj,p,x,y,z,errtol,tiny,big,c1,c2,c3,c4,c5,c6,c7,c8
       integer maxit, nit

       parameter (maxit=50)
c      parameter (errtol=0.05)
       parameter (errtol=0.0015)
       parameter (tiny=2.5e-13,big=9.e11,
     +            c1=3./14.,c2=1./3.,c3=3./22.,c4=3./26.,
     +            c5=0.75*c3,c6=1.5*c4,c7=0.5*c2,c8=c3+c3)

c      errtol=0.05 guarantees correct results in single precision (7 digits)
c      errtol=0.0015 guarantees correct results in double precision (16 digits)
c      tiny is at least twice the cube root of the machine underflow limit
c      big is at most 1/5 the cube root of the machine overflow limit

       real*8 a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,
     +      ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,
     +      xt,yt,zt,rc,rf
       if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.tiny.or.
     +    max(x,y,z,abs(p)).gt.big)then
        write(6,*)'invalid arguments in function rj: ', x, y, z, p
        stop
       endif

       sum=0.
       fac=1.
       if(p.gt.0.)then
        xt=x
        yt=y
        zt=z
        pt=p
       else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1./(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
       endif
       nit = 0

 1     continue
       sqrtx=dsqrt(xt)
       sqrty=dsqrt(yt)
       sqrtz=dsqrt(zt)
       alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
       alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
       beta=pt*(pt+alamb)**2
       sum=sum+fac*rc(alpha,beta)
       fac=0.25*fac
       xt=0.25*(xt+alamb)
       yt=0.25*(yt+alamb)
       zt=0.25*(zt+alamb)
       pt=0.25*(pt+alamb)
       ave=0.2*(xt+yt+zt+pt+pt)
       delx=(ave-xt)/ave
       dely=(ave-yt)/ave
       delz=(ave-zt)/ave
       delp=(ave-pt)/ave
       nit = nit + 1
       if(max(dabs(delx),dabs(dely),dabs(delz),dabs(delp))
     + .gt.errtol .and. nit.lt.maxit) goto 1
       if(nit.ge.maxit)write(6,*)'Carlson rj:', nit, ' iterations', x, y, z, p
       
       ea=delx*(dely+delz)+dely*delz
       eb=delx*dely*delz
       ec=delp**2
       ed=ea-3.*ec
       ee=eb+2.*delp*(ea-ec)
       rj=3.*sum+fac*(1.+ed*(-c1+c5*ed-c6*ee)
     +    +eb*(c7+delp*(-c8+delp*c4))
     +    +delp*ea*(c2-delp*c3)-c2*delp*ec)
     +    /(ave*dsqrt(ave))
       if(p.lt.0.)rj=a*(b*rj+3.*(rcx-rf(xt,yt,zt)))

       return
       end

c
c      Carlson's degenerate elliptic integral Rc
c      "Numerical Recipes" - second edition p. 259
c      Extension to 7th order: B.C.Carlson, Numerical Algorithms 10 (1995) 13-26.
c
       function rc(x,y)
       implicit none
       real*8 rc,x,y,errtol,tiny,sqrtny,big,tnbg,comp1,comp2,third,
     +      c1,c2,c3,c4,c5,c6
       integer maxit, nit

       parameter (maxit=50)
c      parameter (errtol=0.04)
       parameter (errtol=0.0012)
       parameter (tiny=1.69e-38,sqrtny=1.3e-19,big=3.e37,
     +            tnbg=tiny*big,comp1=2.236/sqrtny,comp2=tnbg*tnbg/25.,
     +            third=1./3.,c1=0.3,c2=1./7.,c3=0.375,c4=9./22.
     +           ,c5=159./208.,c6=1.125)

c      errtol=0.05 guarantees correct results in single precision (7 digits)
c      errtol=0.0012 guarantees correct results in single precision (16 digits)
c      tiny is at least 5 times the machine underflow limit
c      big is at most 1/5 the machine overflow limit

       real*8 alamb,ave,s,w,xt,yt
       if(x.lt.0..or.y.eq.0..or.(x+dabs(y)).lt.tiny
     +    .or.(x+dabs(y)).gt.big
     +    .or.(y.lt.-comp1.and.x.gt.0..and.x.lt.comp2)) then
        write(6,*)'invalid arguments in rc: ', x, y
        stop
       endif

       if(y.gt.0.)then
        xt=x
        yt=y
        w=1.
       else
        xt=x-y
        yt=-y
        w=dsqrt(x)/dsqrt(xt)
       endif
       nit = 0

 1     continue
       alamb=2.*dsqrt(xt)*dsqrt(yt)+yt
       xt=0.25*(xt+alamb)
       yt=0.25*(yt+alamb)
       ave=third*(xt+yt+yt)
       s=(yt-ave)/ave
       nit = nit + 1
       if(dabs(s).gt.errtol .and. nit.lt.maxit) goto 1
       if(nit.ge.maxit)write(6,*)'Carlson rc:', nit, ' iterations', x, y

c      Extended to order 7:
       rc=w*(1.+s*s*(c1+s*(c2+s*(c3+s*(c4+s*(c5+s*c6)))))) / dsqrt(ave)
       
       return
       end
c
c      Legendre's elliptic integral of the first kind
c      integral(1/sqrt(1-(ak*sint)^2)dt from t=0 to t=phi
c      (Numerical Recipes, 2nd Edition p.260)
c
       function ellf(phi,ak)
       implicit none
       real*8 ellf,ak,phi
c      uses rf
c      argument ranges: 0<=phi<=pi/2; 0<=k*sin(phi)<=1.
       real*8 s, rf, sak
       s = dsin(phi)
       sak = s * ak
       ellf = s*rf(dcos(phi)**2,(1.d0-sak)*(1.d0+sak),1.d0)
       return
       end
c
c      Legendre's elliptic integral of the second kind
c      integral(sqrt(1-(ak*sint)^2)dt from t=0 to t=phi
c      (Numerical Recipes, 2nd Edition p.261)
c
       function elle(phi,ak)
       implicit none
       real*8 elle,ak,phi
c      uses rd,rf
c      argument ranges: 0<=phi<=pi/2; 0<=k*sin(phi)<=1.
       real*8 cc, q, s, sak, third, rd, rf
       parameter(third=1./3.)
       s = dsin(phi)
       sak = s * ak
       cc = dcos(phi)**2
       q = (1.d0-sak)*(1.d0+sak)
       elle = s*(rf(cc,q,1.d0)-(sak**2)*rd(cc,q,1.d0)*third)
       return
       end
c
c      Legendre's elliptic integral of the third kind
c      Careful with the notations: (en opposite to Abramowitz & Stegun)
c      integral(1./((1+en*(sint)^2)*sqrt(1-(ak*sint)^2))dt from t=0 to t=phi
c      (Numerical Recipes, 2nd Edition p.261)
c
       function ellpi(phi,en,ak)
       implicit none
       real*8 ellpi,ak,en,phi
c      uses rf,rj
c      argument ranges: 0<=phi<=pi/2; 0<=k*sin(phi)<=1.
       real*8 cc, enss, q, s, sak, third, rf, rj
       parameter(third=1./3.)
       s = dsin(phi)
       sak = s * ak
       enss = en*s*s
       cc = dcos(phi)**2
       q = (1.d0-sak)*(1.d0+sak)
       ellpi = s*(rf(cc,q,1.d0)-enss*rj(cc,q,1.d0,1.d0+enss)*third)
       return
       end
c
c      Legendre's complete elliptic integral of the first kind
c      integral(1/sqrt(1-(ak*sint)^2)dt from t=0 to t=pi/2
c
       function ellfc(ak)
       implicit none
       real*8 ellfc, ak
c      uses rf
c      argument ranges: 0<=k<=1.
       real*8 rf
       ellfc = rf(0.d0,(1.d0-ak)*(1.d0+ak),1.d0)
       return
       end
c
c      Legendre's complete elliptic integral of the second kind
c      integral(sqrt(1-(ak*sint)^2)dt from t=0 to t=pi/2
c
       function ellec(ak)
       implicit none
       real*8 ellec, ak
c      uses rd,rf
c      argument ranges: 0<=k<=1.
       real*8 q, third, rd, rf
       parameter(third=1./3.)
         if(ak.ge.1.d0)then
         ellec = 0.
         else
       q = (1.d0-ak)*(1.d0+ak)
       ellec = rf(0.d0,q,1.d0) - (ak**2) * rd(0.d0,q,1.d0)*third
         end if
       return
       end
c
c      Legendre's complete elliptic integral of the first and second kinds,
c      grouped:
c      integral(1/sqrt(1-(ak*sint)^2)dt from t=0 to t=pi/2
c      integral(sqrt(1-(ak*sint)^2)dt from t=0 to t=pi/2
c
       subroutine ellfec(ak, efc, eec)
       implicit none
       real*8 ak, efc, eec
c      uses rd,rf
c      argument ranges: 0<=k<=1.
       real*8 q, third, rd, rf
       parameter(third=1./3.)
         if(ak.ge.1.d0)then
         efc = 0.
         eec = 0.
         else
         q = (1.d0-ak)*(1.d0+ak)
         efc = rf(0.d0,q,1.d0)
         eec = efc-(ak**2)*rd(0.d0,q,1.d0)*third
         end if
       return
       end
c
c      Legendre's complete elliptic integral of the third kind
c      Careful with the notations: (en opposite to Abramowitz & Stegun)
c      integral(1./((1+en*(sint)^2)*sqrt(1-(ak*sint)^2))dt from t=0 to t=pi/2
c
       function ellpic(en,ak)
       implicit none
       real*8 ellpic,ak,en
c      uses rf,rj
c      argument ranges: 0<=k<=1.
       real*8 q, third, rf, rj
       parameter(third=1./3.)
         if(ak.ge.1.d0)then
         ellpic = 0.
         else
       q = (1.d0-ak)*(1.d0+ak)
       ellpic = rf(0.d0,q,1.d0)-en*rj(0.d0,q,1.d0,1.d0+en)*third
         end if
       return
       end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           Extensions - being tested!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      Carlson's degenerate elliptic integral Rc
c      "Numerical Recipes" - second edition p. 259 &
c      B.C.Carlson, Numerical Algorithms 10 (1995) 13-26.
c      Case of complex arguments
c
       function cdrc(x,y)
       implicit none
       real*8 errtol,tiny,sqrtny,big,tnbg,comp1,comp2,third,
     +      c1,c2,c3,c4,c5,c6
       complex*16 cdrc, x, y, z
       integer maxit, nit

       parameter (maxit=50)
c      parameter (errtol=0.04)
       parameter (errtol=0.0012)
       parameter (tiny=1.69e-38,sqrtny=1.3e-19,big=3.e37,
     +            tnbg=tiny*big,comp1=2.236/sqrtny,comp2=tnbg*tnbg/25.,
     +            third=1./3.,c1=0.3,c2=1./7.,c3=0.375,c4=9./22.
     +           ,c5=159./208.,c6=1.125)

c      errtol=0.05 guarantees correct results in single precision (7 digits)
c      errtol=0.0012 guarantees correct results in single precision (16 digits)
c      tiny is at least 5 times the machine underflow limit
c      big is at most 1/5 the machine overflow limit

       real*8 cdabs1
       complex*16 xt, yt, w, alamb, ave, s

       cdabs1(z) = dabs(dreal(z)) + dabs(dimag(z))

c?:
c       if(x.lt.0..or.y.eq.0..or.(x+dabs(y)).lt.tiny
c     +    .or.(x+dabs(y)).gt.big
c     +    .or.(y.lt.-comp1.and.x.gt.0..and.x.lt.comp2)) then
c        write(6,*)'invalid arguments in rc'
c        stop
c       endif

c?:
       if(dreal(y).gt.0.)then
        xt=x
        yt=y
        w=1.
       else
        xt=x-y
        yt=-y
        w=cdsqrt(x)/cdsqrt(xt)
       endif
       nit = 0

 1     continue
       alamb=2.*cdsqrt(xt)*cdsqrt(yt)+yt
       xt=0.25*(xt+alamb)
       yt=0.25*(yt+alamb)
       ave=third*(xt+yt+yt)
       s=(yt-ave)/ave
       nit = nit + 1
       if(cdabs1(s).gt.errtol .and. nit.lt.maxit) goto 1
       if(nit.ge.maxit)write(6,*)'Carlson cdrc:', nit, ' iterations', x, y

       cdrc=w*(1.+s*s*(c1+s*(c2+s*(c3+s*(c4+s*(c5+s*c6)))))) / cdsqrt(ave)
       
       return
       end
c
c      RJ function: Carlson's elliptic integral of the third kind
c      "Numerical Recipes" - second edition, p.257 &
c      B.C.Carlson, Numerical Algorithms 10 (1995) 13-26.
c      Case of a complex parameter p - Experimental!
c
       function cdrj(x,y,z,p)
       implicit none
       real*8 x,y,z,errtol,tiny,big,c1,c2,c3,c4,c5,c6,c7,c8
       complex*16 cdrj, p
       integer maxit, nit

       parameter (maxit=50)
c      parameter (errtol=0.05)
       parameter (errtol=0.0015)
       parameter (tiny=2.5e-13,big=9.e11,
     +            c1=3./14.,c2=1./3.,c3=3./22.,c4=3./26.,
     +            c5=0.75*c3,c6=1.5*c4,c7=0.5*c2,c8=c3+c3)

c      errtol=0.05 guarantees correct results in single precision (7 digits)
c      errtol=0.0015 guarantees correct results in double precision (16 digits)
c      tiny is at least twice the cube root of the machine underflow limit
c      big is at most 1/5 the cube root of the machine overflow limit

       real*8 alamb,ap,
     +      fac,sqrtx,sqrty,sqrtz,
     +      xt,yt,zt,rf
     +, cdabs1

       complex*16 a, b, pt, tau, alpha, beta, rcx, rho, sum
     +, ave, delx, dely, delz, delp, ea,eb,ec,ed,ee
     +, cdrc
       external cdrc, rf

       cdabs1(a) = dabs(dreal(a)) + dabs(dimag(a))

       ap = cdabs1(p)
       if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,ap).lt.tiny.or.
     +    max(x,y,z,ap).gt.big)then
        write(6,*)'invalid arguments in function cdrj: ', x, y, z, p
        stop
       endif

       sum=0.d0
       fac=1.d0
c?:
       if(dreal(p).gt.0.)then
        xt=x
        yt=y
        zt=z
        pt=p
       else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1.d0/(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=cdrc(rho,tau)
       endif
       nit = 0

 1     continue
       sqrtx=dsqrt(xt)
       sqrty=dsqrt(yt)
       sqrtz=dsqrt(zt)
       alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
       alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
       beta=pt*(pt+alamb)**2
       sum=sum+fac*cdrc(alpha,beta)
       fac=0.25*fac
       xt=0.25*(xt+alamb)
       yt=0.25*(yt+alamb)
       zt=0.25*(zt+alamb)
       pt=0.25*(pt+alamb)
       ave=0.2*(xt+yt+zt+pt+pt)
       delx=(ave-xt)/ave
       dely=(ave-yt)/ave
       delz=(ave-zt)/ave
       delp=(ave-pt)/ave
       nit = nit + 1
       if(max(cdabs1(delx),cdabs1(dely),cdabs1(delz),cdabs1(delp))
     + .gt.errtol .and. nit.lt.maxit) goto 1
       if(nit.ge.maxit)write(6,*)'Carlson cdrj:', nit, ' iterations', x, y, z, p
       
       ea=delx*(dely+delz)+dely*delz
       eb=delx*dely*delz
       ec=delp**2
       ed=ea-3.*ec
       ee=eb+2.*delp*(ea-ec)
       cdrj=3.*sum+fac*(1.+ed*(-c1+c5*ed-c6*ee)
     +    +eb*(c7+delp*(-c8+delp*c4))
     +    +delp*ea*(c2-delp*c3)-c2*delp*ec)
     +    /(ave*cdsqrt(ave))

       if(dreal(p).lt.0.)cdrj=a*(b*cdrj+3.*(rcx-rf(xt,yt,zt)))

       return
       end
c
c      Legendre's complete elliptic integral of the third kind
c      case of a complex parameter.
c      Careful with the notations: (en opposite to Abramowitz & Stegun)
c      integral(1./((1+en*(sint)^2)*sqrt(1-(ak*sint)^2))dt from t=0 to t=pi/2
c      Based on Carlson's function cdrj, tentative complex translation.
c
       function cdellpic(en,ak)
       implicit none
       complex*16 cdellpic, en, enp1, z
       real*8 ak, cdabs1, tiny
c      uses cdrf,cdrj
c      argument ranges: 0<=k<=1.
       real*8 q, third, rf
       complex*16 cdrj
       parameter(third=1./3.,tiny=2.5e-13)
       cdabs1(z) = dabs(dreal(z)) + dabs(dimag(z))

       q = (1.d0-ak)*(1.d0+ak)
       enp1 = en + 1.d0
         if(cdabs1(enp1).lt.tiny)then
         write(6,*)'cdellpic: invalid first argument=',en
         return
         end if
       cdellpic = rf(0.d0,q,1.d0)-en*cdrj(0.d0,q,1.d0,1.d0+en)*third
       return
       end
c
c     Legendre's complete elliptic integral of the third kind,
c     case of a complex parameter.
c     Careful with the notations: (en opposite to Abramowitz & Stegun)
c     integral(1./((1+en*(sint)^2)*sqrt(1-(ak*sint)^2))dt from t=0 to t=pi/2
c
c     Based on Byrd & Friedman 416.00, 417.00 p.231 + 2 errata:
c     Lang & Stevens, Math.Comp.14 (1960) 195, and
c     Sutherland, Math.Comp.19 (1965) 132.
c     (Thanks to Bill Carlson!)
c
      function cdelpc(en,ak)
      implicit none
      real*8 ak
      complex*16 cdelpc, en
c     uses ellpic
c     argument ranges: 0<=k<=1.
      real*8 ak2, r2, al12, al22, b1, b2
     ;, ellfc, ellpic
     ;, g1, g2, r2i, a, bp, c, del, m1, m2, t, n1, n2, h2i
     ;, s1, s2, t1, t2, pi, tau2
     
      parameter(pi=3.14159265359d0)
      
      g1 = dreal(en)
      g2 = dimag(en)
      r2 = g1 * g1 + g2 * g2

        if(r2 .eq. 0.d0)then
        cdelpc = dcmplx(ellfc(ak), 0.d0)
        return
        else if(g2 .eq. 0.d0)then
        cdelpc = dcmplx(ellpic(g1,ak), 0.d0)
        return
        end if

      ak2 = ak * ak
      r2i = 1.d0 / r2

      a = r2 + ak2 * (2.*g1 + 1.)
      bp = (1. - ak2) * r2
      c = - r2 * (r2 + 2.*g1 + ak2)
      del = dsqrt(bp * bp - a * c)
      m1 = (- bp + del) / a
      m2 = (- bp - del) / a
        if(dabs(m1) .gt. dabs(m2))then
        t = m1
        m1 = m2
        m2 = t
        end if       
c        if(ak .eq. 0.d0)then
c        end if       
      b1 = m1 * m1 * r2i
      b2 = m2 * m2 * r2i
      
      al12 = ak2 * b1
      al22 = ak2 * b2
c     Lang & Stevens's notations:
      n1 = (m1*(al12*(al12-2.-m1)+(1.+2.*m1)*ak2) - r2) 
     ;     / (r2+(2.*g1+al12)*al12)
      n2 = (m2*(al22*(al22-2.-m2)+(1.+2.*m2)*ak2) - r2) 
     ;     / (r2+(2.*g1+al22)*al22)
c Alternately,
c      n1 = (m1*(b1*(al12-2.-m1)+1.+2.*m1)*ak2 - r2) 
c     ;     / (r2+(2.*g1+al12)*al12)
c      n2 = (m2*(b2*(al22-2.-m2)+1.+2.*m2)*ak2 - r2) 
c     ;     / (r2+(2.*g1+al22)*al22)

      s1 = m1 - n1 - 1.
      s2 = m2 - n2 - 1.
       
      t1 = ((m1 + g1 + 2. - al12) * m1 + n1 * (g1 + al12) + g1) / g2
      t2 = ((m2 + g1 + 2. - al22) * m2 + n2 * (g1 + al22) + g1) / g2

c      q = ak2 * (1. - ak2)
c      o = (ak2*(ak2+2.*g1) + r2) / q
c      h1 = (b1 - 1.) * ak2 * o

c     Suitable for k=0:
      h2i = (1. - ak2) / ((b2 - 1.) * (ak2*(ak2+2.*g1) + r2))

c     Sutherland's correction:
        if(m2.lt.-1.)then
        tau2 = pi * m2 * dsqrt(h2i)
        else if(m2 .eq. -1.)then
        tau2 = 0.5 * pi * m2 * dsqrt(h2i)
        else
        tau2 = 0.
        end if  

      cdelpc = (
     ; dcmplx(t1-t2, s2-s1) * ellfc(ak) + 
     ; dcmplx(t1, -s1) * (n2 * ellpic(-al22,ak) - tau2) + 
     ; n1 * dcmplx(-t2,s2) * ellpic(-al12,ak)
     ; ) / (s1 * t2 - s2 * t1)

      return
      end
