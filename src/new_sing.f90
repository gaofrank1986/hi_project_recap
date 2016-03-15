    subroutine guig_sing(is,ielem,nodj,xp,yp,zp,valdg)
        use mvar_mod
        use proj_cnst,only:ex,ey,xiqet,xiqsi,cross_product
        use green_funcs,only:Dy3DGFunc
        use mfunc_mod,only:spfunc8_1

        implicit none

        integer,intent(in) :: is,ielem,nodj
        real(8),intent(in) :: xp,yp,zp
        real(8),intent(out) :: valdg(8)

        real(8) :: nodes(8,3),nrmls(8,3),eta0,si0
        real(8) :: sf0(8),dsf0(2,8),ddsf0(3,8)
        real(8) :: sf(8),dsf(2,8),ddsf(3,8)
        real(8) :: jk0(3,3),n0(3),xj(2,3),xxj(3,3)
        real(8) :: v1,xjp(2,3),p(3),np(3)
        real(8) :: total(8),fitg,si,eta
        real(8) :: Xiq(8),Wiq(8),rho1,csst,snst,p0(3)
        real(8) :: det,f1,f2,param(2),f
        
        integer :: nsamb,loop1,loop2,j,i

        DATA XIQ/ 0.960289856497536D+00, 0.796666477413626D+00, &
            0.525532409916329D+00, 0.183434642495650D+00, &
            -0.183434642495650D+00,-0.525532409916329D+00, &
            -0.796666477413626D+00,-0.960289856497536D+00/

        DATA WIQ/ 0.101228536290376D+00, 0.222381034453374D+00, &
            0.313706645877887D+00, 0.362683783378362D+00, &
            0.362683783378362D+00, 0.313706645877887D+00, &
            0.222381034453374D+00, 0.101228536290376D+00/     
        !nrml_id = ncond(ielem,nodj)

        do i =1,ncn(ielem)
            nodes(i,:)=xyz(:,ncon(ielem,i))
            nrmls(i,:)=dxyz(:,ncond(ielem,i))
        end do
        
   

        si0=xiqsi(nodj)
        eta0=xiqet(nodj)
        p0=(/xp,yp,zp/)
        call spfunc8_1(si0,eta0,sf0,dsf0,ddsf0)
        xj = matmul(dsf0,nodes)
        xxj = matmul(ddsf0,nodes)

        jk0(1,:)=cross_product(xj(1,:),xj(2,:))
        jk0(2,:) = cross_product(xxj(1,:),xj(2,:))+cross_product(xj(1,:),xxj(3,:))
        jk0(3,:) = cross_product(xxj(3,:),xj(2,:))+cross_product(xj(1,:),xxj(2,:))
  !      do i=1,3
        !write(*,'(3f14.8)') jk0(i,:) 
        !enddo
        !pause


        do j=1,ncn(ielem)
            n0=(/sf0(j),dsf0(1,j),dsf0(2,j)/)
            call line_part(nodj,ncn(ielem),si0,eta0,xj,xxj,jk0,n0,fitg)!jk0,jk1c,jk1s,n0,n1c,n1s,xj,xxj)
            valdg(j)=fitg
        end do
        write (*,'(8f14.8)') valdg
        
        pause
        nsamb=0

        do loop1=1,8
            do loop2=1,8
                nsamb=nsamb+1
                si = xiq(loop1)
                eta = xiq(loop2)

                rho1 = norm2((/si0,eta0/)-(/si,eta/))

                csst = (si-si0)/rho1
                snst = (eta-eta0)/rho1

                call spfunc8_1(si,eta,sf,dsf,ddsf)

                p=matmul(sf,nodes)
                np=matmul(sf,nrmls)
                v1=dot_product(np,Dy3DGFunc(p,p0))

                xjp = matmul(dsf,nodes)
                det = norm2(cross_product(xjp(1,:),xjp(2,:)))

                !jkt1:3) = jk0(1:3)+rho*(jk1c(1:3)*csst+jk1s(1:3)*snst)
                do j =1,ncn(ielem) ! 8
                    !n0=sf0(j)
                    !n1c=dsf0(1,j)
                    !n1s=dsf0(2,j)
                    n0=(/sf0(j),dsf0(1,j),dsf0(2,j)/)

                    call comp_coef(csst,snst,xj,xxj,jk0,n0,f1,f2,param)


                    f=f2/rho1**2+f1/rho1
                    !fixme
                    !valdg(j)=(v1*det*sf(j))*wiq(loop1)*wiq(loop2) !/rho to change back to xi coordinate
                    !valdg(j)=valdg(j)+(v1*det*sf(j)-f/rho1)*wiq(loop1)*wiq(loop2) !/rho to change back to xi coordinate
                    valdg(j)=(v1*det*sf(j))*wiq(loop1)*wiq(loop2) !/rho to change back to xi coordinate
                end do
            end do
        end do

        write (*,'(8f14.8)') valdg
        pause
        !valdg=total

        end subroutine

        subroutine line_part(nodj,ncne,si0,eta0,xj,xxj,jk0,n0,fitg) 

            use proj_cnst,only: new8
            implicit none

            integer, intent(in) :: nodj,ncne
            real(8),intent(in) :: si0,eta0,jk0(3,3),n0(3),xj(2,3),xxj(3,3)
            real(8),intent(out) :: fitg 

            real(8) ::  xiqsi(32),xiqet(32),wiq1(17),wiq2(25),f1,f2,param(2)

            real(8) :: si,eta,gammas,betas,fitg1,fitg2,csst,snst,rho1
            integer :: newnum,id,ind,i

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
                0.122489331563432,0.062177497273381/! to be d\theta


            data wiq2/0.122489331563432,0.231823804500403,0.199261222833210,  &
                0.160875277198321,0.126277137889030,0.098697779924940,  &
                0.077797413988515,0.062177497273381,0.080187721987975,  &
                0.109334472936971,0.117207837861905,0.122489331563432,  &
                0.124354994546761,0.122489331563432,0.117207837861905,  &
                0.109334472936971,0.080187721987975,0.062177497273381,  &
                0.077797413988515,0.098697779924940,0.126277137889030,  &
                0.160875277198321,0.199261222833210,0.231823804500403,  &
                0.122489331563432/

            newnum = new8(nodj)! give a new number

            select case(nodj)

            case(1,3,5,7)!90deg
                do i=1,17
                    ind=mod(newnum+i+6,32)+1
                    si = xiqsi(ind)
                    eta = xiqet(ind)
                    

                    rho1=norm2((/si,eta/)-(/si0,eta0/))!si,eta is on the edge
                    csst=(si-si0)/rho1
                    snst=(eta-eta0)/rho1
                    

                 !   print *,si,eta,si0,eta0
                    !print *,rho1
                    !print *,csst,snst

                    call comp_coef(csst,snst,xj,xxj,jk0,n0,f1,f2,param)
                   ! print *,"checking"
                    !print *,f1,f2
                    !print *, param
                    !pause
                    gammas=param(1)
                    betas=param(2)

                    fitg2=fitg2-f2*(gammas/betas**2+1.0/rho1)*wiq1(i)
                    fitg1=fitg1+f1*dlog(rho1/betas)*wiq1(i)
                    !sum1=sum1+ f1*wiq1(i)
                    !sum2=sum2+f2/betas*wiq1(i)
                end do
            case (2,4,6,8)!180 deg
                do i =1,25
                    ind=mod(newnum+i+2,32)+1
                    si=xiqsi(ind)
                    eta=xiqet(ind)

                    rho1=norm2((/si,eta/)-(/si0,eta0/))!si,eta is on the edge
                    csst=(si-si0)/rho1
                    snst=(eta-eta0)/rho1
                 
                    call comp_coef(csst,snst,xj,xxj,jk0,n0,f1,f2,param)
                    gammas=param(1)
                    betas=param(2)

                    fitg2=fitg2-f2*(gammas/betas**2+1.0/rho1)*wiq2(i)
                    fitg1=fitg1+f1*dlog(rho1/betas)*wiq2(i)

                end do
            end select
            fitg =fitg1+fitg2

            !!!no negtive sign
        end subroutine

        !   @param csst
        !   @param snst
        !   @param xj
        !   @param xxj
        !   @param jk0
        !   @param n0
        !   @return f1,f2
        !   @return param
        subroutine comp_coef(csst,snst,xj,xxj,jk0,n0,f1,f2,param) 
            
            use proj_cnst,only :pi4

            implicit none
            
            real(8),intent(in) :: csst,snst,xj(2,3),xxj(3,3),jk0(3,3),n0(3)
            real(8),intent(out) :: f1,f2,param(2)

            real(8) :: jk1(3),n1
            real(8) :: a(3),b(3),ast,bst,as_3,as_2,g31,b30,a30,b31,a31

            jk1(1:3) = jk0(2,1:3)*csst+jk0(3,1:3)*snst
            n1=n0(2)*csst+n0(3)*snst

            a(1:3) = xj(1,:)*csst+xj(2,:)*snst
            b(1:3) = xxj(1,:)*csst**2*0.5+xxj(2,:)*snst**2*0.5+xxj(3,:)*csst*snst

            ast=norm2(a)
            bst=norm2(b)

            as_3=1.0d0/ast**3
            as_2=-3.0d0*(dot_product(a,b))/ast**5
            g31=dot_product(b,jk0(1,:))+dot_product(a,jk1)
            g31=(a(3)/ast**2)*g31

            b30=-jk0(1,3)
            a30=b30*n0(1)

            b31=3*g31-jk1(3)
            a31=b31*n0(1)+b30*n1

            f1=-(as_2*a30+as_3*a31)/pi4
            f2=-(as_3*a30)/pi4

            !gammas = -dot_product(a,b)/ast**4
            !betas=1.0d0/ast

            param(1) = -dot_product(a,b)/ast**4
            param(2) = 1.0d0/ast

        end subroutine
