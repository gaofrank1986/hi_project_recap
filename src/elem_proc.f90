    !< @file 
    !
    !
    !!
    subroutine norm_elem_wrapper(ielem,xp,yp,zp,aval,bval,hi) 
        use mvar_mod
        implicit none
        integer :: is,ielem,hi
        real(8) :: p0(3),aval(4,8),bval(4,8),xp,yp,zp
        
        do is=1,nsys
           call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval,bval,hi)
           !print *,ielem,'norm'
           !write (*,'(f14.8)') aval(1,8)
           !write (*,'(f14.8)') bval(1,8)
       end do 
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
        real*8  aval(4,8),bval(4,8)

        !print *,"enter gauss_int"
        aval=0.0d0
        bval=0.0d0

        p0 = (/ex(is)*xp,ey(is)*yp,zp/)
        
        nsamb=16
        if(ncne.eq.6)   nsamb=4

        do n=1,   nsamb     

            p =sambxy(ielem,n,1:3)
            np = dsamb(ielem,n,1:3)
            if (hi.eq.1) then
                v1 = GFunc(p,p0)+GFunc(p,mirror(h,p0))
                v2 = dot_product(np,DGFunc(p,p0)+DGFunc(p,mirror(h,p0)))
            elseif (hi.eq.2) then
                v1 = Dy3GFunc(p,p0)+Dy3GFunc(p,mirror(h,p0))
                v2 = dot_product(np,Dy3DGFunc(p,p0)+Dy3DGFunc(p,mirror(h,p0)))
            elseif (hi.eq.9) then
                !second part of sing_int1
                v2 = Dy3GFunc(p,mirror(h,p0))
                v1 = dot_product(np,Dy3DGFunc(p,mirror(h,p0)))
            end if

            do  j=1,   ncne
                bval(is,j)=bval(is,j)+v1*samb(ielem,n,j)
                aval(is,j)=aval(is,j)+v2*samb(ielem,n,j)
            enddo


        end do
    end subroutine




        subroutine sing_elem_wrapper(inode,ielem,numqua,xp,yp,zp,aval,bval,hi)
            use mvar_mod
            use mfunc_mod,only:tripol
            !use tripole_mod
            implicit   none  
            integer i,j,is,ielem,inode,nodnum,nd,np,numqua,hi
            real*8  xp,yp,zp,xyzt(3,8),dxyzt(3,8)
            real*8 aval(4,8),bval(4,8)
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
                        call sing_int(is,ielem,nodnum,xp,yp,zp,aval,bval,hi)
                    else if(is.ne.1 ) then   
                        call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval,bval,hi)
                    end if
                100      continue

            else if(numqua.eq.2) then
                do 200 is=1,nsys     
                    if(is.eq.1.or.is.eq.2) then
                        call sing_int(is,ielem,nodnum,xp,yp,zp,aval,bval,hi)
                    else if(is.eq.3.or.is.eq.4) then  
                        call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval,bval,hi)
                    end if
                200      continue  

            else if(numqua.eq.4) then
                do 300  is=1,  nsys
                    if(is.eq.1.or.is.eq.4) then   
                        call sing_int(is,ielem,nodnum,xp,yp,zp,aval,bval,hi)
                    else if(is.eq.2.or.is.eq.3) then  
                        call gauss_int(is,ielem,ncn(ielem),xp,yp,zp,aval,bval,hi)
                    end if
                300      continue

            else if(numqua.eq.5) then
                do 400 is=1, nsys  
                    call sing_int(is,ielem,nodnum,xp,yp,zp,aval,bval,hi)
                400      continue
            endif
                !
            return
        end

    subroutine sing_int(is,ielem,nodej,xp,yp,zp,aval,bval,hi)
        implicit none
        integer,intent(in) :: is,ielem,hi,nodej
        real(8),intent(in) :: xp,yp,zp
        real(8),intent(out) :: aval(4,8),bval(4,8)

        if (hi.eq.1) then
            call sing_int0(is,ielem,xp,yp,zp,aval,bval)
        elseif (hi.eq.2) then
            call sing_int1(is,ielem,nodej,xp,yp,zp,aval,bval)
        endif
    end subroutine

!
    subroutine sing_int0(is,ielem,xp,yp,zp,aval,bval)


        use trvar_mod
        use mvar_mod
        use green_funcs,only : GFunc,DGFunc,mirror
        use proj_cnst,only:ex,ey
        implicit   none  

        real(8) :: p(3),p0(3),np(3)
        integer is,ielem,n,j,ip       
        real*8  xp,yp,zp
        real*8  v1,v2
        real*8  aval(4,8),bval(4,8)

        p0 = (/ex(is)*xp,ey(is)*yp,zp/)
        aval=0.0d0
        bval=0.0d0

        do  n=1, nosamp

            p =xynod(1:3,n)
            np = dxynod(1:3,n)

            v1 = gfunc(p,p0)+gfunc(p,mirror(h,p0))
            v2 = dot_product(np,dgfunc(p,p0)+dgfunc(p,mirror(h,p0)))

            do j=1, ncn(ielem)
                aval(is,j)=aval(is,j)+v2*samnod(n,j)
                bval(is,j)=bval(is,j)+v1*samnod(n,j)
            enddo
        end do

        end



        subroutine sing_int1(is,ielem,nodj,xp,yp,zp,amatrix,bmatrix) 

            use mvar_mod
            use trvar_mod    
            use mfunc_mod
            use hi_intg
            use proj_cnst,only :ex,ey,xiqsi,xiqet

            implicit none

            integer,intent(in):: is,ielem,nodj
            real(8),intent(in)::  xp,yp,zp
            real(8),intent(out):: bmatrix(4,8),amatrix(4,8)

            real(8) :: src_lcl(2),src_glb(3),origin_offset(3)
            real(8) :: cnr_glb_mtx(3,8),cnr_glb_nrml(3,8)

            real(8) :: result0(8),result1(8)
            real(8) ::  si,eta,p0(3)



            if(ncn(ielem).eq.8)  then 

                si =xiqsi(nodj) !get local coordinate for the src
                eta=xiqet(nodj)

            else if(ncn(ielem).eq.6)  then
                !si =xitsi(nodj)
                !eta=xitet(nodj)
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
            call eval_singular_elem(cnr_glb_mtx,result0,result1)
            ! the mirrored source is evaluated use gauss intg
            call gauss_int(is,ielem,8,xp,yp,zp,amatrix,bmatrix,9) 

            amatrix(is,:) = amatrix(is,:)+result0(:)
            bmatrix(is,:) = bmatrix(is,:)+result1(:)


        end subroutine 






