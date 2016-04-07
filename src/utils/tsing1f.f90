    !
    ! *************************************************************
    ! *                                                           *
    ! *  The source point is in the mesh                          *
    ! *                                                           *
    ! *************************************************************
    !
    subroutine sgbd0_1(is,ielem,nodn,p0,valg,valdg) 
        use kinds
        use mesh,only:xyz,dxyz,ncn,ncon,ncond
        use shape_funcs
        use proj_cnst,only:ex,ey,xiqet,xiqsi
        use linalg,only:cross_product
        use green_funcs,only:gcombo1_1

        implicit none

        integer,intent(in)      ::  is,ielem,nodn   
        real(rk),intent(in)     ::  p0(3) 
        real(rk),intent(out)    ::  valg(8),valdg(8)

        integer loop1,loop2,nsamb
        integer i,j

        real(rk)  sf0(8),dsf0(2,8),ddsf0(3,8) 
        real(rk)  sf(8),dsf(2,8)
        real(rk)  xj(2,3),xjp(2,3),xxj(3,3)

        real(rk) :: det,ans
        real(rk)  xxx(3,8),xxd(3,8)
        real(rk) ::xi0(2),p(3),np(3),xi(2)
        real(rk)  tot(8),totf(8)
        real(rk)  f,f1,f2,vs
        real(rk)  plo,theta(2)
        real(rk)  gxf0(4),dgn0,n01(3),jk01(3,3),param(2)

        !< --- parameter
        real(rk)  xitsi(6),xitet(6)
        real(rk)  xiq(8),wiq(8),xit(7),wit(7)

        data xiq/ 0.960289856497536d+00, 0.796666477413626d+00, &
            0.525532409916329d+00, 0.183434642495650d+00, &
            -0.183434642495650d+00,-0.525532409916329d+00, &
            -0.796666477413626d+00,-0.960289856497536d+00/

        data wiq/ 0.101228536290376d+00, 0.222381034453374d+00, &
            0.313706645877887d+00, 0.362683783378362d+00, &
            0.362683783378362d+00, 0.313706645877887d+00, &
            0.222381034453374d+00, 0.101228536290376d+00/     

        !<---------------------------------------------------

        valg=0.0d0 
        valdg=0.0d0
        tot=0.0d0
        totf=0.0d0

        if(ncn(ielem).eq.8)  then 

            xi0=[xiqsi(nodn),xiqet(nodn)]
            call spfunc8_1(xi0(1),xi0(2),sf0,dsf0,ddsf0) 

            do i =1,8
                xxx(1:3,i) = xyz(1:3,ncon(ielem,i))
                xxd(1:3,i) = dxyz(1:3,ncond(ielem,i))
                print '(3f10.5)',xxd(1:3,i)
            end do

        endif

        !< shape func 1st,2nd order derivatives at xi0
        xj(1:2,1:3) = matmul(dsf0(1:2,1:8),transpose(xxx))
        xxj(1:3,1:3) = matmul(ddsf0(1:3,1:8),transpose(xxx))
        !< jk0,jk01c,jk01c
        jk01(1,:) = cross_product(xj(1,:),xj(2,:))
        jk01(2,:) = cross_product(xxj(1,:),xj(2,:))+cross_product(xj(1,:),xxj(3,:))
        jk01(3,:) = cross_product(xxj(3,:),xj(2,:))+cross_product(xj(1,:),xxj(2,:))

        !< line integral part
        do  j=1, ncn(ielem)
            n01=[sf0(j),dsf0(1,j),dsf0(2,j)]
            call cirbod_2(nodn,ncn(ielem),xi0,jk01,n01,xj,xxj,ans)
            valdg(j) = ans
        enddo
        

        !< area integral part
        nsamb=0
        if(ncn(ielem).eq.8)  then 
            do  loop1=1, 8 ;do  loop2=1, 8      
                nsamb=nsamb+1  

                xi=[xiq(loop1),xiq(loop2)]
                call spfunc8(xi(1),xi(2),sf,dsf)
                plo=norm2(xi-xi0)
                theta=(xi-xi0)/plo

                ! < gp global pos, gp corresponding normal
                p=matmul(sf,transpose(xxx))
                np=matmul(sf,transpose(xxd))

                call gcombo1_1(-1.0d0,p,p0,gxf0)
                dgn0=dot_product(gxf0(2:4),np)

                xjp(1:2,1:3)=matmul(dsf,transpose(xxx))
                det =norm2(cross_product(xjp(1,:),xjp(2,:)))

                do  j=1, ncn(ielem)
                    n01=[sf0(j),dsf0(1,j),dsf0(2,j)]
                    call comp_coef(theta,xj,xxj,jk01,n01,f1,f2,param)
                    print '(2f10.5)',f1,f2

                    f=f2/plo/plo+f1/plo                       ! (c19)

                    tot(j)=tot(j)+f/plo*wiq(loop1)*wiq(loop2)
                    !todo added minus
                    totf(j)=totf(j)+dgn0*det*wiq(loop1)*wiq(loop2)*sf(j)


                enddo
            

            enddo; enddo
            ! fixme add vs to compensate sign flip
            vs=sign(1.0d0,totf(nodn)*tot(nodn))
            write (*,'(a10,8f10.5)') "line part",valdg(:)
            write (*,'(a10,8f10.5)') "tot part",tot(:)
            write (*,'(a10,8f10.5)') "totf part",totf(:)

            valdg=vs*valdg+totf-vs*tot

       endif

   end subroutine

   !
   !C
   !C **************************************************************
   !C *                                                            *
   !C *  Calculating the one-dimension integral part of the        *
   !C *  singular integral                                         * 
   !C *  The singularity is on the body surface                    *
   !C *     by M. Guiggiani and A. Gigante's Method (Dec. 1990)    *
   !C *                                                            *
   !C **************************************************************
   !C
   subroutine cirbod_2(nodn,ncne,xi0,jk01,n01,xj,xxj,ans)
       use kinds
       implicit none 

       integer newnum,nodn,in,i,ncne
       integer new6(6),new8(8)

       real(rk) xiqsi(32),xiqet(32),wiq1(17),wiq2(25)
       real(rk) xitsi(24),xitet(24),wit1(9),wit2(17)
       real(rk) plo,ans
       real(rk) f1,f2       
       real(rk) xi(2),xi0(2),theta(2),xj(2,3),xxj(3,3)
       real(rk) jk01(3,3),n01(3),param(2)
       !    
       ! ** sampling points and weighting factors for quadrilateral element
       !
       !   8 nodes at each side
       !       
       data new8/1, 5, 9, 13, 17, 21, 25, 29/
       data xiqsi/-1.d0,-.75d0,-.5d0,-.25d0,0.d0,.25d0,.5d0,.75d0,1.d0,   &
           1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,    &
           .75d0,.5d0,.25d0,0.d0,-.25d0,-.5d0,-.75d0,-1.d0,   &
           -1.d0, -1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0/
       data xiqet/-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,  &
           -.75d0,-.5d0,-.25d0,0.d0,.25d0,.5d0,.75d0,1.d0,    &
           1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,   &
           .75d0,.5d0,.25d0,0.d0,-.25d0,-.5d0,-.75d0/ 

       data wiq1/0.062177497273381,0.122489331563432,0.117207837861905,  &
           0.109334472936971,0.099914322536495,0.089926749896239,  &
           0.080115342139031,0.070948527302082,0.066568163775824,  &
           0.070948527302082,0.080115342139031,0.089926749896239,  &
           0.099914322536495,0.109334472936971,0.117207837861905, &
           0.122489331563432,0.062177497273381/


       data wiq2/0.122489331563432,0.231823804500403,0.199261222833210,  &
           0.160875277198321,0.126277137889030,0.098697779924940,  &
           0.077797413988515,0.062177497273381,0.080187721987975,  &
           0.109334472936971,0.117207837861905,0.122489331563432,  &
           0.124354994546761,0.122489331563432,0.117207837861905,  &
           0.109334472936971,0.080187721987975,0.062177497273381,  &
           0.077797413988515,0.098697779924940,0.126277137889030,  &
           0.160875277198321,0.199261222833210,0.231823804500403,  &
           0.122489331563432/
       !    
       ! ** sampling points and weighting factors for triangular element
       !   8 nodes at each side
       !
       data new6/1, 9, 17, 5, 13, 21/     
       data xitsi/0.0d0, 0.125d0, 0.25d0, 0.375d0, 0.5d0, 0.615d0, 0.75d0, 0.875d0,   & 
           1.0d0, 0.875d0, 0.75d0, 0.615d0, 0.5d0, 0.375d0, 0.25d0, 0.125d0,   &
           0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0  /
       data xitet/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,       &
           0.125d0, 0.25d0, 0.375d0, 0.5d0, 0.615d0, 0.75d0, 0.875d0, 1.0d0,   &
           0.875d0, 0.75d0, 0.615d0, 0.5d0, 0.375d0, 0.25d0, 0.125d0/


       ans=0.0d0


       !  =====================================================================         
       ! **** FOR QUADRILATERIAL ELEMENTS  **************************
       !       
       if (ncne .eq. 8) then 
           newnum=new8(nodn) 
           select case(nodn)
           case(1,3,5,7)
               do   i=1, 17!two edges

                   in=mod(newnum+i+6, 32) + 1 
                   xi=[xiqsi(in),xiqet(in)]
                   plo=norm2(xi-xi0)
                   theta = (xi-xi0)/plo
                   call comp_coef(theta,xj,xxj,jk01,n01,f1,f2,param)

                   ans=ans-f2*(param(1)/param(2)**2+1.0/plo)*wiq1(i)
                   ans=ans+f1*dlog(plo/param(2))*wiq1(i) 

               end do

           case(2,4,6,8)
               do  i=1, 25!3edges

                   in=mod(newnum+i+2, 32) + 1 
                   xi=[xiqsi(in),xiqet(in)]
                   plo=norm2(xi-xi0)
                   theta = (xi-xi0)/plo
                   call comp_coef(theta,xj,xxj,jk01,n01,f1,f2,param)

                   ans=ans-f2*(param(1)/param(2)**2+1.0/plo)*wiq2(i)
                   ans=ans+f1*dlog(plo/param(2))*wiq2(i)
               enddo
           end select

       end if             


   end subroutine



   subroutine comp_coef(theta,xj,xxj,jk0,n0,f1,f2,param) 

       use proj_cnst,only :pi4

       implicit none

       real(8),intent(in) :: theta(2),xj(2,3),xxj(3,3),jk0(3,3),n0(3)
       real(8),intent(out) :: f1,f2,param(2)

       real(8) :: jk1(3),n1  
       real(8) :: a(3),b(3),ast,bst,as_3,as_2,g31,b30,a30,b31,a31

       !theta=[csst,snst]
       print '(2f10.5)',theta
       print '(3f10.5)',xj(1,:)
       print '(3f10.5)',xj(2,:)
       print '(3f10.5)',xxj(1,:)
       print '(3f10.5)',xxj(2,:)
       print '(3f10.5)',xxj(3,:)
       print '(3f10.5)',n0(:)
       print '(3f10.5)',jk0(1,:)
       print '(3f10.5)',jk0(2,:)
       print '(3f10.5)',jk0(3,:)

       !jk1(1:3) = jk0(2,1:3)*csst+jk0(3,1:3)*snst
       jk1(1:3) = matmul(theta,jk0(2:3,:))

       !n1=n0(2)*csst+n0(3)*snst
       n1=dot_product(n0(2:3),theta)

       !a(1:3) = xj(1,:)*theta(1)+xj(2,:)*theta(2)
       !b(1:3) = xxj(1,:)*theta(1)**2*0.5+xxj(2,:)*theta(2)**2*0.5+xxj(3,:)*theta(1)*theta(2)
       a(1:3) = matmul(theta,xj(1:2,:))
       b(1:3) = matmul([theta(1)**2*0.5,theta(2)**2*0.5,theta(1)*theta(2)],xxj(1:3,:))

       ast=norm2(a)
       bst=norm2(b)

       print '(3f15.5)',a
       print '(3f15.5)',b
       print '(2f15.5)',bst,ast

       as_3=1.0d0/ast**3
       as_2=-3.0d0*(dot_product(a,b))/ast**5
       g31=dot_product(b,jk0(1,:))+dot_product(a,jk1)
       g31=(a(3)/ast**2)*g31

       b30=-jk0(1,3)
       a30=b30*n0(1)

       b31=3*g31-jk1(3)
       a31=b31*n0(1)+b30*n1

       !!todo neg removed
       !f1=-(as_2*a30+as_3*a31)/pi4
       !f2=-(as_3*a30)/pi4
       f1=-(as_2*a30+as_3*a31)/pi4
       f2=-(as_3*a30)/pi4

       param(1) = -dot_product(a,b)/ast**4
       param(2) = 1.0d0/ast
       stop

   end subroutine
