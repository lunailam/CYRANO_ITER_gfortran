      complex*16 function zif4(zarg, zkll, onlyab)

      implicit none
      logical onlyab
      real*8 zarg, zkll
c
c version of the zif1 function optimized for cray
c plasma dispersion function for real values of k// and omega.
c k//=0 limit to be treated elsewhere.
c zif2 is i=sqrt(-1) times the f0 function of stix (theory of plasma
c waves 1962 , eq. (32), sect.8.8). in the case k// is positive
c zif2 is the fried&conte dispersion function.
c input arguments:
c     raymond's convention on vll:
c                 zarg     argument (omega-n*omc)/sqrt(2)*kll*vll
c     (in cyrano, vth = sqrt(2)*vll )
c                 zkll     k//
c                 onlyab   if .true., only compute the imaginary part
c                          
      integer ini, i, irange, n
     ;, isrchfgt 
      real*8 x, x2, brkpts(6), sq2, fn, term, dd, rx2, xabs
      complex*16
     ;  facps(63), facas(28), powers(63), zdotu
     ;, irtpi, irtpi2, sum

      external isrchfgt, zdotu
c
      data irtpi/ (0.d0,1.7724 53851d0)/,sq2/1.4142 13562 37309 5d0/
     ;,    irtpi2/ (0.d0,3.5449 07702d0)/
      data ini/0/
      data brkpts/4.d0,5.3,10.,50.,1.d3,1.d6/
      save ini

        if(ini.eq.0) then
        ini=1
        fn=0.
        term=1.
          do 180 i=1,63
          fn=fn+1.d0
          term=term/fn
180       facps(i)= term/(fn+fn+1.d0)
        fn=0.
        term=1.
          do 190 i=1,28
          fn=fn+1.d0
          term=term*(fn+fn-1.d0)
190       facas(i)= term
        end if
c should write another function for easy limit k//=0. 
        if( zkll.eq.0.d0 ) then
        zif4=0.d0
 
        else
        x = zarg
        x2 = x * x
        xabs = dabs(x)
        irange = isrchfgt(6, brkpts, 1, xabs)
 
          if(onlyab)then
c         --------------
c         Only compute imaginary part
            if( irange .eq. 1 )then
            dd = dexp(- x2)
            zif4 = dsign(dd, zkll) * irtpi
c              if( zkll .lt. 0.d0 )then
c               zif4=-dd*irtpi
c               else
c               zif4=dd*irtpi
c               end if
            return
            end if
            if(x2 .gt. 170.d0) then
            dd = 0.d0
            else
            dd = dexp(- x2)
            end if
          zif4 = dsign(dd, zkll) * irtpi
          else
c         ----
c         Normal case: compute dispersion function
            if( irange .eq. 1 )then
            dd = dexp(- x2)

            n=6.2+xabs*(5.75+2.125*xabs)
            powers(1)=x2
              do 220 i=2,n
220           powers(i)=x2*powers(i-1)
            sum=(1.d0,0.d0)+zdotu(n,facps,1,powers,1)
 
c             choose the determination of zif4
              if( zkll .lt. 0.d0 )then
              zif4=-dd*(irtpi+2.d0*x*sum)
              else
              zif4=dd*(irtpi-2.d0*x*sum)
              end if
c            write(6,255)irange,n,sum,     zarg,zif4
c255         format(' irange=',i1,i3  ,1p,8d11.3)
            return
c
            else if( irange .eq. 2 )then
            n=xabs*xabs
            else if( irange .eq. 3 )then
            n=54.5/(xabs-3.243)
            else if( irange .eq. 4 )then
            n=8
            else if( irange .eq. 5 )then
            n=4
            else if( irange .eq. 6 )then
            n=2
            else if( irange .eq. 7 )then
            sum=1.
            n=0
            end if

            if(x2 .gt. 170.d0) then
            dd = 0.d0
            else
            dd = dexp(- x2)
            end if
 
            if(n.ne.0) then
            x2=0.5d0/x2
            powers(1)=x2
              do 390 i=2,n
390           powers(i)=x2*powers(i-1)
            sum = (1.,0.d0) + zdotu(n, facas, 1, powers, 1)
            end if
          zif4=-sum/x
          zif4=zif4+dsign(1.d0,zkll)*irtpi*dd
          end if
        end if
c       ------

      return
 
901      format(' ',10('='),'    default zarg**2=',1p,e11.3,
     1          '  has been taken: potentially meaningless result',
     2          /' ',100('=')/)
      end
