    !< @file 
    !
    !
    !!
    subroutine norm_elem_wrapper(inode,ielem,xp,yp,zp,aval,bval,hi) 
        use mvar_mod
        !use sieppem,only:sieppema
        use hi_intg,only:swap_result
        implicit none
        integer :: is,ielem,hi,inode,i,sel
        real(8) :: p0(3),aval(4,8),bval(4,8),xp,yp,zp,tmp(8),tmp1(8)
        real(8) :: cd(3,8),xis(2,1)

        !xp=xyz(1,117)
        !yp=xyz(2,117)
        !zp=xyz(3,117)
        !ielem=401
        !!xis(:,1)=(/-1.0d0,-1.0d0/)
        !!do i=1,8
            !!cd(:,i) = xyz(1:3,ncon(ielem,i))
        !!end do
            !cd(:,1) = xyz(:,ncon(ielem,1))
            !cd(:,2) = xyz(:,ncon(ielem,3))
            !cd(:,3) = xyz(:,ncon(ielem,5))
            !cd(:,4) = xyz(:,ncon(ielem,7))
            !cd(:,5) = xyz(:,ncon(ielem,2))
            !cd(:,6) = xyz(:,ncon(ielem,4))
            !cd(:,7) = xyz(:,ncon(ielem,6))
            !cd(:,8) = xyz(:,ncon(ielem,8))
        !p0=(/xp,yp,zp/)
        !do i=1,8
            !write(*,'(3f14.8)') cd(:,i)
        !end do
        !print *,p0
        !sel=0

       
        
        do is=1,nsys
           call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval(is,:),bval(is,:),hi)
           !call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval(is,:),bval(is,:),1)
           !call sieppema(2.0d0,cd,p0,tmp,xis,sel)
           !call swap_result(tmp)
           
           !call direct_eval(is,ielem,ncn(ielem),xp,yp,zp,tmp,tmp1,hi)
           !write (*,'(3i5,8f16.8)') inode,ielem,1,aval(is,:)
           !write (*,'(3i5,8f16.8)') inode,ielem,2,bval(is,:)
           !write (*,'(3i5,8f16.8)') inode,ielem,1,tmp
           !write (*,'(3i5,8f16.8)') inode,ielem,2,tmp1
       end do 
    end subroutine


    subroutine sing_elem_wrapper(inode,ielem,numqua,xp,yp,zp,aval,bval,hi)
        use mvar_mod
        use mfunc_mod,only:tripol
        !use tripole_mod
        implicit   none  
        integer i,j,is,ielem,inode,nodnum,nd,np,numqua,hi
        real*8  xp,yp,zp,xyzt(3,8),dxyzt(3,8)
        real*8 aval(4,8),bval(4,8),tmp(8),tmp1(8)
        !print *,'sing',ielem

        if (hi.eq.1) then

            do  i=1,  ncn(ielem)
                xyzt(1:3, i)  =  xyze(1:3, i, ielem)  
                dxyzt(1:3, i) = dxyze(1:3, i, ielem)  
                if(inode.eq.ncon(ielem,i)) nodnum=i
            end do
            call tripol(nodnum,ncn(ielem),xyzt,dxyzt)
        elseif (hi.eq.2) then
            do  i=1,  ncn(ielem)
                if(inode.eq.ncon(ielem,i)) nodnum=i
            end do
        end if

        if(numqua.eq.0)       then
            do 100 is=1,  nsys
                if(is.eq.1) then 
                    call sing_int(is,ielem,nodnum,xp,yp,zp,aval(is,:),bval(is,:),hi)
                    !call direct_eval(is,ielem,ncn(ielem),xp,yp,zp,tmp,tmp1,hi)

                    !write (8000,'(3i5,8f16.8)') inode,ielem,1,aval(is,:)
                    !write (8000,'(3i5,8f16.8)') inode,ielem,2,bval(is,:)
                    !write (8000,'(3i5,8f16.8)') inode,ielem,1,tmp
                    !write (8000,'(3i5,8f16.8)') inode,ielem,2,tmp1
                    !stop
                else if(is.ne.1 ) then   
                    call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval(is,:),bval(is,:),hi)
                end if
            100      continue

        else if(numqua.eq.2) then
            do 200 is=1,nsys     
                if(is.eq.1.or.is.eq.2) then
                    call sing_int(is,ielem,nodnum,xp,yp,zp,aval(is,:),bval(is,:),hi)
                else if(is.eq.3.or.is.eq.4) then  
                    call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval(is,:),bval(is,:),hi)
                end if
                200      continue  

        else if(numqua.eq.4) then
            do 300  is=1,  nsys
                if(is.eq.1.or.is.eq.4) then   
                    call sing_int(is,ielem,nodnum,xp,yp,zp,aval(is,:),bval(is,:),hi)
                else if(is.eq.2.or.is.eq.3) then  
                    call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval(is,:),bval(is,:),hi)
                end if
            300      continue

        else if(numqua.eq.5) then
            do 400 is=1, nsys  
                call sing_int(is,ielem,nodnum,xp,yp,zp,aval(is,:),bval(is,:),hi)
                400      continue
        endif
    end subroutine
   
    subroutine gauss_int(is,ielem,ncne,xp,yp,zp,aval,bval,hi)


        use mvar_mod
        use green_funcs
        use proj_cnst,only:ex,ey
        implicit   none  

        integer is,ielem,n,nsamb,ncne,j,hi

        real*8  xp,yp,zp
        real(8) :: p(3),p0(3),np(3)

        real*8  v1,v2
        real*8  aval(8),bval(8)

        !print *,"enter gauss_int"
        aval=0.0d0
        bval=0.0d0

        p0 = (/ex(is)*xp,ey(is)*yp,zp/)
        
        nsamb=16
        if(ncne.eq.6)   nsamb=4

        do n=1,   nsamb     
            !todo consider h<0,namely,infinite depth
            p =sambxy(ielem,n,1:3)
            np = dsamb(ielem,n,1:3)
            if (hi.eq.1) then
                v1 = GFunc(p,p0)!+GFunc(p,mirror(h,p0))
                v2 = dot_product(np,DGFunc(p,p0))!+DGFunc(p,mirror(h,p0)))
            elseif (hi.eq.2) then
                ! add src point at mirror position about bottom plane
                !v1 = Dy3GFunc(p,p0)+Dy3GFunc(p,mirror(h,p0))
                !v2 = dot_product(np,Dy3DGFunc(p,p0)+Dy3DGFunc(p,mirror(h,p0)))
                
                ! add sink point at mirror positon about bottom plane
                v1 = Dy3GFunc(p,p0)!-Dy3GFunc(p,mirror(h,p0))
                v2 = dot_product(np,Dy3DGFunc(p,p0))!-Dy3DGFunc(p,mirror(h,p0)))

            elseif (hi.eq.7) then

                ! test sing1
                v2 = dot_product(np,Dy3DGFunc(p,p0))
                v1 = dot_product(np,Dy3DGFunc(p,p0))
            elseif (hi.eq.8) then
                v1 = Dy3GFunc(p,p0)
                v2 = 0.
            elseif (hi.eq.9) then
                !second part of sing_int1
                v2 = -Dy3GFunc(p,mirror(h,p0))
                v1 = dot_product(np,-Dy3DGFunc(p,mirror(h,p0)))
            end if

            do  j=1,   ncne
                bval(j)=bval(j)+v1*samb(ielem,n,j)
                aval(j)=aval(j)+v2*samb(ielem,n,j)
            enddo


        end do
    end subroutine





    subroutine sing_int(is,ielem,nodej,xp,yp,zp,aval,bval,hi)
        implicit none
        integer,intent(in) :: is,ielem,hi,nodej
        real(8),intent(in) :: xp,yp,zp
        real(8),intent(out) :: aval(8),bval(8)

        if (hi.eq.1) then
            call sing_int0(is,ielem,xp,yp,zp,aval,bval)
        elseif (hi.eq.2) then
            call sing_int1(is,ielem,nodej,xp,yp,zp,aval,bval)
        endif
    end subroutine

!
    subroutine sing_int0(is,ielem,xp,yp,zp,aval,bval)


        use trvar_mod
        use mvar_mod,only:ncn,h
        use green_funcs,only : GFunc,DGFunc,mirror
        use proj_cnst,only:ex,ey
        implicit   none  

        real(8) :: p(3),p0(3),np(3)
        integer is,ielem,n,j,ip       
        real*8  xp,yp,zp
        real*8  v1,v2
        real*8  aval(8),bval(8)

        p0 = (/ex(is)*xp,ey(is)*yp,zp/)
        aval(:)=0.0d0
        bval(:)=0.0d0

        do  n=1, nosamp

            p =xynod(1:3,n)
            np = dxynod(1:3,n)

            v1 = gfunc(p,p0)!+gfunc(p,mirror(h,p0))
            v2 = dot_product(np,dgfunc(p,p0))!+dgfunc(p,mirror(h,p0)))

            do j=1, ncn(ielem)
                aval(j)=aval(j)+v2*samnod(n,j)
                bval(j)=bval(j)+v1*samnod(n,j)
            enddo
        end do

        end



        subroutine sing_int1(is,ielem,nodj,xp,yp,zp,amatrix,bmatrix) 

            use mvar_mod
            use mfunc_mod
            use hi_intg
            use proj_cnst,only :ex,ey,xiqsi,xiqet

            implicit none

            integer,intent(in):: is,ielem,nodj
            real(8),intent(in)::  xp,yp,zp
            real(8),intent(out):: bmatrix(8),amatrix(8)

            real(8) :: src_lcl(2),src_glb(3),origin_offset(3)
            real(8) :: cnr_glb_mtx(3,8),cnr_glb_nrml(3,8)
            real(8) ::tmp(8),tmp1(8),tmp2(8)

            real(8) :: result0(8),result1(8)
            real(8) ::  si,eta,p0(3)

            amatrix=0
            bmatrix=0

            if(ncn(ielem).eq.8)  then 
                si =xiqsi(nodj) !get local coordinate for the src
                eta=xiqet(nodj)
            endif


            p0 = (/ex(is)*xp,ey(is)*yp,zp/)

            ! if mesh is created by using symmetrical information
            ! basically, we get the same local layout
            ! however, the glb position of the src point is different
            ! this will lead to ...
            ! when the src information is passed to the integration 
            ! both src local pos and global pos is set
            ! they has to match each other, the relationship cannot be reflected here
            ! since both NODJ,and XP-ZP are input information

            src_glb(:) = (/xp,yp,zp/)
            src_lcl(:) = (/si,eta/)
            origin_offset(:) = (/0,0,0/)

            !========switch cnr_glb_mtx order for hi module
            cnr_glb_mtx(:,1) = xyz(:,ncon(ielem,1))
            cnr_glb_mtx(:,2) = xyz(:,ncon(ielem,3))
            cnr_glb_mtx(:,3) = xyz(:,ncon(ielem,5))
            cnr_glb_mtx(:,4) = xyz(:,ncon(ielem,7))
            cnr_glb_mtx(:,5) = xyz(:,ncon(ielem,2))
            cnr_glb_mtx(:,6) = xyz(:,ncon(ielem,4))
            cnr_glb_mtx(:,7) = xyz(:,ncon(ielem,6))
            cnr_glb_mtx(:,8) = xyz(:,ncon(ielem,8))


            cnr_glb_nrml(:,1) = dxyz(:,ncond(ielem,1))
            cnr_glb_nrml(:,2) = dxyz(:,ncond(ielem,3))
            cnr_glb_nrml(:,3) = dxyz(:,ncond(ielem,5))
            cnr_glb_nrml(:,4) = dxyz(:,ncond(ielem,7))
            cnr_glb_nrml(:,5) = dxyz(:,ncond(ielem,2))
            cnr_glb_nrml(:,6) = dxyz(:,ncond(ielem,4))
            cnr_glb_nrml(:,7) = dxyz(:,ncond(ielem,6))
            cnr_glb_nrml(:,8) = dxyz(:,ncond(ielem,8))

            call preset_src(si,eta,xyz(1:3,ncon(ielem,nodj)),origin_offset)
            !call eval_singular_elem(cnr_glb_mtx,result0)
            !amatrix(:)=amatrix(:)+result0
            !write (*,'(a,8f14.8)') "evalu",result0

            !call gauss_int(is,ielem,8,xp,yp,zp,tmp,tmp1,7) 
            !write (*,'(a,8f14.8)') "direct",tmp(:)
            !call direct_eval(is,ielem,8,xp,yp,zp,tmp,tmp1,7) 
            !write (*,'(a,8f14.8)') "direc2",tmp(:)
            !write (*,'(8f14.8)') tmp1(1,:)
            !print *,"test done"
            !pause

            call SGBD0_1(IS,IELEM,NODj,XP,YP,ZP,result0,result1) 
            !write (*,'(a,8f14.8)') "sgbd0",result1
            amatrix=amatrix+result1
            
            !call gauss_int(is,ielem,8,xp,yp,zp,tmp,tmp1,7) 
            ! the mirrored source is evaluated use gauss intg
            !call guig_sing(is,ielem,nodj,xp,yp,zp,tmp2) 
            !amatrix(is,:)=amatrix(is,:)+result0
            
            !todo change variable name
            !todo bmatrix is meaningful
            !tmp1=0.
            !call gauss_int(is,ielem,8,xp,yp,zp,tmp,tmp1,8) 
            !write (*,'(a,8f14.8)') "lower",tmp1(:)
            !bmatrix(:)=bmatrix(:)+tmp1(:)
            

            !call gauss_int(is,ielem,8,xp,yp,zp,tmp,tmp1,9) 
!            write (*,'(a,8f14.8)') "mlowe",tmp1(:)
            !write (*,'(a,8f14.8)') "mhigh",tmp(:)
          !  call direct_eval(is,ielem,8,xp,yp,zp,tmp,tmp1,9) 
            !write (*,'(a,8f14.8)') "mlowe",tmp1(:)
            !write (*,'(a,8f14.8)') "mhigh",tmp(:)
            !stop
            write(9000,*) "ielem=", ielem
            write(9000,*) "amatrix="
            write(9000,'(f14.8)') amatrix(:)
            write(9000,*) "bmatrix="
            write(9000,'(f14.8)') bmatrix(:)

            !amatrix(:) = amatrix(:)+tmp(:)
            !bmatrix(:) = bmatrix(:)+tmp1(:)
            !print *,"Exiting"
            !print *,"H=",h
            !stop


        end subroutine 






    subroutine direct_eval(is,ielem,ncne,xp,yp,zp,aval,bval,hi)
        use mvar_mod
        use green_funcs
        use proj_cnst,only:ex,ey
        implicit   none  

        integer is,ielem,n,nsamb,ncne,j,hi,i

        real*8  xp,yp,zp
        real(8) :: p(3),p0(3),np(3)

        real*8  v1,v2
        real*8  aval(8),bval(8)

        aval=0.0d0
        bval=0.0d0

        p0 = (/ex(is)*xp,ey(is)*yp,zp/)
        

        do  i=1,ncn(ielem)
            p =xyz(1:3,ncon(ielem,i))
            np = dxyz(1:3,ncond(ielem,i))
            
            if (hi.eq.1) then
                v1 = GFunc(p,p0)+GFunc(p,mirror(h,p0))
                v2 = dot_product(np,DGFunc(p,p0)+DGFunc(p,mirror(h,p0)))
            elseif (hi.eq.2) then
                v1 = Dy3GFunc(p,p0)+Dy3GFunc(p,mirror(h,p0))
                v2 = dot_product(np,Dy3DGFunc(p,p0)+Dy3DGFunc(p,mirror(h,p0)))

            elseif (hi.eq.7) then

                ! test sing1
                v2 = dot_product(np,Dy3DGFunc(p,p0))
                v1 = dot_product(np,Dy3DGFunc(p,p0))
            elseif (hi.eq.8) then
                v1 = Dy3GFunc(p,p0)+Dy3GFunc(p,mirror(h,p0))
                v2 = 0.
            elseif (hi.eq.9) then
                !second part of sing_int1
                v2 = Dy3GFunc(p,mirror(h,p0))
                v1 = dot_product(np,Dy3DGFunc(p,mirror(h,p0)))
            end if

                bval(i)=v1
                aval(i)=v2


        end do
    end subroutine
