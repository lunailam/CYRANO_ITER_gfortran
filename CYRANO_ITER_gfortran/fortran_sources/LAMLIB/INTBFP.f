      subroutine intbfp

      implicit none
c
c     Computes analytical integrals of products of 2 basis functions by a power
c     of reduced variable ksi over (0,1). For power -1, the singular entries are
c     set to 1.d20.
c     Results stored in PRODIN as follows:
c     prodin(...,-1): integrals of f1 * f2 / ksi
c     prodin(..., 0): integrals of f1 * f2 * (1 - ksi)
c     prodin(..., 1): integrals of f1 * f2 * ksi
c
c     Also computes analytical integrals of a basis function f times a power of
c     ksi over (0,1): 
c     Results stored in BAFINT as follows:
c     bafint(..., 0): integrals of f * (1 - ksi)
c     bafint(..., 1): integrals of f * ksi
c
      include 'pardim.copy'

      include 'comfin.copy'
      include 'comin2.copy'
c      
c     Basis functions library: NTYP=2 available types are
c     (1) Hermite cubic, (2) Lagrange quadratic
c      
      integer it1, j1, ip1, it2, j2, ip2, i, j, l, id1, id2, ipo

      double precision
     ;  coe(maxdeg+1,maxbaf,ntyp), t1, u1, t2, u2
      
      save coe

c     Duplicates information in BAFADN:      
c     Coefficients of basis functions, for each type:
c     maximum of MAXBAF basis functions for each basis function set type, with maximum
c     degree MAXDEG
c     NTYP types of basis functions, NBAF basis functions for each type;
c     NELTYP types of finite element
      data
     ;  coe/ 1.d0, 0.d0,-3.d0, 2.d0
     ;     , 0.d0, 1.d0,-2.d0, 1.d0
     ;     , 0.d0, 0.d0, 3.d0,-2.d0
     ;     , 0.d0, 0.d0,-1.d0, 1.d0
     ;
     ;     , 1.d0,-3.d0, 2.d0, 0.d0
     ;     , 0.d0, 4.d0,-4.d0, 0.d0
     ;     , 0.d0,-1.d0, 2.d0, 0.d0
     ;     , 4*0.d0
     ;  /

c      call dset((maxbaf*2*ntyp)**2 * 3, 0.d0, prodin(1,0,1,1,0,1,-1), 1)
c      call dset((maxbaf*2*ntyp) * 2, 0.d0, bafint(1,0,1,0), 1)
      call dset((maxbaf*2*ntyp)**2 * 3, 0.d0, prodin, 1)
      call dset((maxbaf*2*ntyp) * 2, 0.d0, bafint, 1)


      do it1 = 1, ntyp                  ! For all basis function types
        do j1 = 1, nbaf(it1)            ! For all basis function belonging to this type
          do ip1 = 0, deg(it1)          ! Loop over terms of polynomial, by ascending degree
          t1 = coe(ip1+1,j1,it1)
          u1 = dfloat(ip1) * coe(ip1+1,j1,it1)

            do it2 = 1, ntyp            ! For all basis function types
              do j2 = 1, nbaf(it2)      ! For all basis function belonging to this type
                do ip2 = 0, deg(it2)    ! Loop over terms of polynomial, by ascending degree
                t2 = coe(ip2+1,j2,it2)
                u2 = dfloat(ip2) * coe(ip2+1,j2,it2)
                  do ipo = -1, 1        ! Loop over powers of reduced variable ksi present in the integrand
                  j = ip1 + ip2 + ipo   ! Degree of power of ksi in current monomial (j>=-1)
                  i = j + 1
                  l = j - 1
                    if(i.gt.0)then  ! add contribution to integral of product of two b.f. times ksi^ipo
                    prodin(j1,0,it1,j2,0,it2,ipo) = prodin(j1,0,it1,j2,0,it2,ipo) + t1 * t2 / dfloat(i)
                    else            ! the case j=-1, i=0: the integral does not exist
                    prodin(j1,0,it1,j2,0,it2,ipo) = 1.d20
                    end if
                    if(j.gt.0)then  ! add contribution to integral of product of a b.f. by derivative of another times ksi^ipo
                    prodin(j1,1,it1,j2,0,it2,ipo) = prodin(j1,1,it1,j2,0,it2,ipo) + u1 * t2 / dfloat(j)
                    prodin(j1,0,it1,j2,1,it2,ipo) = prodin(j1,0,it1,j2,1,it2,ipo) + t1 * u2 / dfloat(j)
                    else            ! the cases j=-1 and j=0
                    if(ip1.ne.0)prodin(j1,1,it1,j2,0,it2,ipo) = 1.d20
                    if(ip2.ne.0)prodin(j1,0,it1,j2,1,it2,ipo) = 1.d20
                    end if
                    if(l.gt.0)then  ! add contribution to integral of product of two b.f. derivatives times ksi^ipo
                    prodin(j1,1,it1,j2,1,it2,ipo) = prodin(j1,1,it1,j2,1,it2,ipo) + u1 * u2 / dfloat(l)
c                    prodin(j1,1,it1,j2,1,it2,ipo) = prodin(j1,1,it1,j2,1,it2,ipo) + u1 * u2 / dfloat(max0(1,l))
                    else if(ip1*ip2 .ne.0)then ! the case l=0
                    prodin(j1,1,it1,j2,1,it2,ipo) = 1.d20
                    end if
                  end do      
                end do      
              end do      
            end do      
          end do      
        end do      
      end do      

c     Now modify prodin(...,0) to store integral of products of b.f. times (1-ksi) there:
      do it1 = 1, ntyp
        do j1 = 1, nbaf(it1)
        do id1 = 0, 1
          do it2 = 1, ntyp
            do j2 = 1, nbaf(it2)
            do id2 = 0, 1
            prodin(j1,id1,it1,j2,id2,it2,0) = prodin(j1,id1,it1,j2,id2,it2,0) - prodin(j1,id1,it1,j2,id2,it2,1)
            end do      
            end do      
          end do      
        end do      
        end do      
      end do      

      do it1 = 1, ntyp
        do j1 = 1, nbaf(it1)
          do ip1 = 0, deg(it1)
          t1 = coe(ip1+1,j1,it1)
          u1 = dfloat(ip1) * coe(ip1+1,j1,it1)
            do ipo = 0, 1
            j = ip1 + ipo
            i = j + 1
            bafint(j1,0,it1,ipo) = bafint(j1,0,it1,ipo) + t1 / dfloat(i)
            bafint(j1,1,it1,ipo) = bafint(j1,1,it1,ipo) + u1 / dfloat(max(j,1))
            end do      
          end do      
        end do      
      end do      

      do it1 = 1, ntyp
        do j1 = 1, nbaf(it1)
          do id1 = 0, 1
          bafint(j1,id1,it1,0) = bafint(j1,id1,it1,0) - bafint(j1,id1,it1,1)
          end do      
        end do      
      end do      
      
      return
      end
