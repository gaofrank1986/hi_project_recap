!<  @file Deal with integration on element,regular or singular
    !< --------------------------------------------------
    !< norm ele1 is wrapper function for norm integration
    !!
    subroutine norm_elem_wrapper(ielem,xp,yp,zp,amatrix,bmatrix,hi)

        use mesh,only:ncn,nsys
        use green_funcs,only:gcombo1,gcombo0

        implicit   none 
        integer,intent(in) :: ielem,hi
        real(8),intent(in) ::  xp,yp,zp 
        real(8),intent(out) ::  bmatrix(4,8),amatrix(4,8)

        integer :: is

        interface
            subroutine gcombo(h,p0,p,gxf)
                real(8),intent(in) :: h,p0(3),p(3)
                real(8),intent(out) :: gxf(4)
            end subroutine
        end interface

        procedure(gcombo),pointer :: gpointer => NULL()

        !< --------------execute part----------------->
        if (hi==1) then
            gpointer => gcombo0
        elseif (hi==2) then
            gpointer => gcombo1
        end if

        do is=1,   nsys  
            call norm_int(is,ielem,ncn(ielem),(/xp,yp,zp/),&
                amatrix(is,:),bmatrix(is,:),gpointer) 
        end do
    end subroutine



    !< 
    !!  @params a_res,b_res return high,low order parts of integral on given element
    !!  @params bie_called select which gfunction combo to be used
    subroutine norm_int(is,ielem,ncne,p0,aval,bval,bie_called)
        use gaussian_info,only:sambxy,dsamb,samb
        use wave,only:h
        use proj_cnst,only : ex,ey
        implicit   none  

        integer,intent(in) :: is,ielem,ncne
        real(8),intent(in) ::  p0(3)!xp,yp,zp
        real(8),intent(out) :: aval(8),bval(8)

        interface
            subroutine bie_called(h,p0,p,gxf)
                real(8),intent(in) :: h,p0(3),p(3)
                real(8),intent(out) :: gxf(4)
            end subroutine
        end interface

        real(8) :: dgn,gxf(4),p0m(3),p(3),prefix(3),np(3)
        integer :: n,nsamb,j

        aval=0.0d0
        bval=0.0d0

        nsamb=16
        !if(ncne.eq.6)   nsamb=4
        prefix=(/ex(is),ey(is),1.0d0/)
        p0m = prefix*p0!(/xp,yp,zp/)

        do  n=1,   nsamb     

            !< get gaussian point
            p=sambxy(ielem,n,1:3)

            call bie_called(h,p,p0,gxf)

            !< get nrml at gaussian point
            np = prefix*dsamb(ielem,n,1:3) 
            dgn = dot_product(gxf(2:4),np)

            bval(:)=bval(:)+gxf(1)*samb(ielem,n,1:8)
            aval(:)=aval(:)+dgn*samb(ielem,n,1:8)

        enddo
    end subroutine


    subroutine sing_elem_wrapper(inode,ielem,numqua,xp,yp,zp,amatrix,bmatrix,hi)
        use mvar_mod
        use green_funcs,only:gcombo1,gcombo0
        use mfunc_mod,only:tripol

        implicit   none  

        integer,intent(in) :: inode,ielem,numqua,hi

        integer i,j,is,nodnum,nd,np,nsamb
        real*8  xp,yp,zp,xyzt(3,8),dxyzt(3,8)
        real*8 bmatrix(4,8),amatrix(4,8)

        interface
            subroutine gcombo(h,p0,p,gxf)
                real(8),intent(in) :: h,p0(3),p(3)
                real(8),intent(out) :: gxf(4)
            end subroutine
        end interface

        procedure(gcombo),pointer :: gpointer => NULL()

        !< --------------execute part----------------->
        if (hi==2) then
            gpointer => gcombo1
        elseif (hi==1) then
            gpointer => gcombo0
        end if

        !< compute xynod,dxynod for weak singular elem
        !
        if (hi==1) then
            do  i=1,  ncn(ielem)
                xyzt(1:3, i)  =  xyze(1:3, i, ielem)  
                dxyzt(1:3, i) = dxyze(1:3, i, ielem)  

                if(inode.eq.ncon(ielem,i)) nodnum=i
            enddo
            call tripol(nodnum,ncn(ielem),xyzt,dxyzt)
        end if
            
        do   i=1,  ncn(ielem)
            if(inode.eq.ncon(ielem,i)) nodnum=i ! get node num of inode
        enddo



        !!!!attention !!!!!! sing_int1 didn't affect matrix A
        ! it is ok, bmatrix is almost 0 fort surface boudary

        if(numqua.eq.0)       then
            do is=1,  nsys
                if(is.eq.1) then 
                    call sing_int(is,ielem,nodnum,xp,yp,zp,amatrix,bmatrix,hi)
                else if( is.ne.1 ) then   
                    !call norm_int(is,ielem,ncn(ielem),xp,yp,zp,amatrix,bmatrix,gpointer)
                    ! write(*,*) 'after subroutine sgwp0_1'
                end if
            end do

        else if(numqua.eq.2) then
            do is=1,nsys     
                if(is.eq.1.or.is.eq.2) then  
                    call sing_int(is,ielem,nodnum,xp,yp,zp,amatrix,bmatrix,hi)
                else if(is.eq.3.or.is.eq.4) then  
                    !call norm_int(is,ielem,ncn(ielem),xp,yp,zp,amatrix,bmatrix,gpointer) 
                end if
            end do

        else if(numqua.eq.4) then
            do   is=1,  nsys
                if(is.eq.1.or.is.eq.4) then   
                    call sing_int(is,ielem,nodnum,xp,yp,zp,amatrix,bmatrix,hi)
                else  if(is.eq.2.or.is.eq.3) then  
                    !call norm_int(is,ielem,ncn(ielem),xp,yp,zp,amatrix,bmatrix,gpointer) 
                end if
            end do
        else if(numqua.eq.5) then
            do  is=1, nsys  
                call sing_int(is,ielem,nodnum,xp,yp,zp,amatrix,bmatrix,hi)
            end do
        endif

        end subroutine

    subroutine sing_int(is,ielem,nodnum,xp,yp,zp,amatrix,bmatrix,hi)
        use mvar_mod,only:ncn,xyze,dxyze
        implicit none
        integer,intent(in) :: is,ielem,nodnum,hi
        real(8),intent(in) :: xp,yp,zp
        real(8) :: amatrix(4,8),bmatrix(4,8)
        
        real(8) :: p0(3)
        
        p0=(/xp,yp,zp/)
        if (hi.eq.2) then
            call sing_int1(is,ielem,nodnum,p0,amatrix(is,:),bmatrix(is,:)) 
        elseif (hi.eq.1) then
            call sing_int0(is,ielem,nodnum,p0,amatrix(is,:),bmatrix(is,:)) 
        end if

    end subroutine


    !<  -----------------------------------------------
    !<  element integration with src point on one node
    !!  method only capable of Gfunc desingularizaito


    subroutine sing_int0(is,ielem,nodnum,p0,aval,bval)
        use mvar_mod
        use green_funcs,only:gcombo0
        use proj_cnst,only:ex,ey
        use trvar_mod,only:xynod,dxynod,nosamp,samnod
        implicit   none  

        integer is,ielem,n,j,ip,nodnum
        real*8  xp,yp,zp!,ex(4,4),ey(4,4)
        real*8  x,y,z,x0,y0,z0      
        real*8  nx,ny,nz,dgn
        real*8  aval(8),bval(8),gxf(4),p(3),p0(3),np(3),p0m(3),prefix(3)

        prefix=[ex(is),ey(is),1.0d0]
        !< mapped p0
        p0m = prefix*p0

        aval=0.0d0
        bval=0.0d0

        do n=1, nosamp

            p= xynod(1:3,n)
            np=dxynod(1:3,n)

            call gcombo0 (h,p,p0m,gxf) 
            dgn = dot_product(gxf(2:4),np)

            aval(:)=aval(:)+dgn*samnod(n,1:8)
            bval(:)=bval(:)+gxf(1)*samnod(n,1:8)
        enddo
            !if(ielem>nelemf) then
            !write(9011,'(2i6,8f12.6)') ielem,nodnum,aval(is,:)
            !call sing_int1(is,ielem,nodnum,xp,yp,zp,aval,bval) 
            !end if
            !
    end subroutine


    !   @params p0 input src point
    subroutine sing_int1(is,ielem,nodj,p0,aval,bval) 

        use mvar_mod
        use mfunc_mod
        use hi_intg
        use green_funcs,only:Gcombo1_2
        use proj_cnst,only: ex,ey,xiqsi,xiqet

        implicit none

        integer,intent(in):: is,ielem,nodj
        real(8),intent(in):: p0(3)
        real(8),intent(out):: bval(8),aval(8)

        integer ::inode,inodd,j,pwr_g

        real(8) :: src_lcl(2),src_glb(3),origin_offset(3),xi0(2),p0m(3)
        real(8) :: cnr_glb_mtx(3,8)
        real(8) :: passed_nrml(3,8)

        real(8) :: result0(8),result1(8),result2(8),prefix(3)
        real(8) ::  x0,y0,z0,si,eta,xp,yp,zp



        inode=ncon(ielem,nodj) ! corresponding node id
        inodd=ncond(ielem,nodj)!                 normal id
        aval=0.0d0
        bval=0.0d0

        if(ncn(ielem).eq.8)  then 

            si =xiqsi(nodj) !get local coordinate for the src
            eta=xiqet(nodj)

        else if(ncn(ielem).eq.6)  then
            !si =xitsi(nodj)
            !eta=xitet(nodj)
        endif

        ! get src point in sysmetric mesh
        prefix=(/ex(is),ey(is),1.0d0/)
        p0m=prefix*p0

        
        ! if mesh is created by using symmetrical information
        ! basically, we get the same local layout
        ! however, the glb position of the src point is different
        ! this will lead to ...
        ! when the src information is passed to the integration 
        ! both src local pos and global pos is set
        ! they has to match each other, the relationship cannot be reflected here
        ! since both NODJ,and XP-ZP are input information

        xi0 = (/si,eta/)
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

        passed_nrml(:,1) = dxyz(:,ncond(ielem,1))
        passed_nrml(:,2) = dxyz(:,ncond(ielem,3))
        passed_nrml(:,3) = dxyz(:,ncond(ielem,5))
        passed_nrml(:,4) = dxyz(:,ncond(ielem,7))
        passed_nrml(:,5) = dxyz(:,ncond(ielem,2))
        passed_nrml(:,6) = dxyz(:,ncond(ielem,4))
        passed_nrml(:,7) = dxyz(:,ncond(ielem,6))
        passed_nrml(:,8) = dxyz(:,ncond(ielem,8))


        ! add mirrored sink
        ! Gcombo1_2 is the mirrored sink part only
        ! attention p0 is used here,not p0m
        call norm_int(is,ielem,8,p0,aval,bval,Gcombo1_2)
        !write(9013,'(2i6,8f12.6)') ielem,nodj,amatrix(is,:)

        !< note p03 is used here
        call preset_src(si,eta,p0m,origin_offset)
        call eval_singular_elem(cnr_glb_mtx,passed_nrml,result0,result1,result2)

        !aval(:) = amatrix+result0 !< elemental plus
        !bmatrix(:) = bmatrix(:)+result1!elemental plus
        
        aval(:) = aval+result0 !< elemental plus
        bval(:) = bval+result1 !< elemental plus
            
        write(9010,'(2i6,8f12.6)') ielem,nodj,result2


    end subroutine 
         

