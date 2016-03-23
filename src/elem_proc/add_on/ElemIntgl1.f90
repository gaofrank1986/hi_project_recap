!
!  NORM_ELE1+SING_ELE1
!  NORM_INT1+NORM_INT1
!
! ======================================================
!
!   Integration on an element without source in itself
!   nor its symmetrical ones
!
! ======================================================
!
    !< norm ele1 is wrapper function for norm integration
    !!
    subroutine norm_ele1(ielem,xp,yp,zp,amatrix,bmatrix)

        use mvar_mod,only:ncn,nsys
        use green_funcs,only:gcombo1

        implicit   none 
        integer,intent(in) :: ielem
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
        gpointer => gcombo1

        bmatrix=0.0d0
        amatrix=0.0d0

        do is=1,   nsys  
            call norm_int(is,ielem,ncn(ielem),xp,yp,zp,&
                amatrix(is,:),bmatrix(is,:),gpointer) 
        end do
        end subroutine

!
! ======================================================
!   Integration on an element with source in itself
!   or its mirror ones about any symmetrical axis
!
! ======================================================
!
    subroutine sing_ele1(inode,ielem,numqua,xp,yp,zp,amatrix,bmatrix)
        use mvar_mod
        use green_funcs,only:gcombo1,gcombo0

        implicit   none  

        integer i,j,is,ielem,inode,nodnum,nd,np,nsamb,numqua,hi
        real*8  xp,yp,zp,xyzt(3,8),dxyzt(3,8)
        real*8 bmatrix(4,8),amatrix(4,8)

        interface
            subroutine gcombo(h,p0,p,gxf)
                real(8),intent(in) :: h,p0(3),p(3)
                real(8),intent(out) :: gxf(4)
            end subroutine
        end interface

        procedure(gcombo),pointer :: gpointer => NULL()

        hi=2
        if (hi==2) then
            gpointer => gcombo1
        elseif (hi==1) then
            gpointer => gcombo0
        end if

        do   i=1,  ncn(ielem)
            if(inode.eq.ncon(ielem,i)) nodnum=i ! get node num of inode
        enddo


        bmatrix= 0.0d0     
        amatrix= 0.0d0  

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
                !call sing_int(is,ielem,nodnum,xp,yp,zp,amatrix,bmatrix,hi)
            end do
        endif

        end subroutine

    subroutine sing_int(is,ielem,nodnum,xp,yp,zp,amatrix,bmatrix,hi)
        implicit none
        integer,intent(in) :: is,ielem,nodnum,hi
        real(8),intent(in) :: xp,yp,zp
        real(8) :: amatrix(4,8),bmatrix(4,8)
        
        real(8) :: p0(3)
        
        p0=(/xp,yp,zp/)
        if (hi.eq.2) then
            call sing_int1(is,ielem,nodnum,p0,amatrix,bmatrix) 
        elseif (hi.eq.1) then
            !call sing_int0(is,ielem,nodnum,p0,amatrix,bmatrix) 
        end if

    end subroutine

        
    !   < 
    !!  @params bie_called select which gfunction combo to be used
    subroutine norm_int(is,ielem,ncne,xp,yp,zp,a_res,b_res,bie_called)
        use mvar_mod
        use proj_cnst,only : ex,ey
        implicit   none  

        integer,intent(in) :: is,ielem,ncne
        real(8),intent(in) ::  xp,yp,zp
        real(8),intent(out) :: a_res(8),b_res(8)

        interface
            subroutine bie_called(h,p0,p,gxf)
                real(8),intent(in) :: h,p0(3),p(3)
                real(8),intent(out) :: gxf(4)
            end subroutine
        end interface

        real(8) :: dgn,gxf(4),p0(3),p(3),prefix(3),np(3)
        integer :: n,nsamb,j

        a_res=0.0d0
        b_res=0.0d0

        nsamb=16
        !if(ncne.eq.6)   nsamb=4
        prefix=(/ex(is),ey(is),1.0d0/)
        p0 = prefix*(/xp,yp,zp/)

        do  n=1,   nsamb     

            p=sambxy(ielem,n,1:3)

            call bie_called(h,p,p0,gxf)

            np = prefix*dsamb(ielem,n,1:3) 
            dgn = dot_product(gxf(2:4),np)

            do   j=1,   ncne
                b_res(j)=b_res(j)+gxf(1)*samb(ielem,n,j)
                a_res(j)=a_res(j)+dgn*samb(ielem,n,j)
            enddo

        enddo
    end subroutine
    !
    !
    ! *************************************************************
    ! *                                                           *
    ! *  The source point is in the mesh, at the node NODJ        *
    ! *                                                           *
    ! *************************************************************
    !
    !    function reorder_node_value(t_kind) result(h_kind)
    !implicit none
    !real(8),intent(in) :: t_kind(:,:) 
    !integer,dimension(2) :: s = shape(t_kind)
    !real(8) :: h_kind(s(1),s(2))!:,:)



    !!integer :: s(2) = shape(t_kind)
    !if (s(2).ne.8) then
    !print *, " try to reorder a non-8-node nodal values"
    !pause
    !end if

    !h_kind(:,1) = t_kind(:,ncon(ielem,1))
    !h_kind(:,2) = t_kind(:,ncon(ielem,3))
    !h_kind(:,3) = t_kind(:,ncon(ielem,5))
    !h_kind(:,4) = t_kind(:,ncon(ielem,7))
    !h_kind(:,5) = t_kind(:,ncon(ielem,2))
    !h_kind(:,6) = t_kind(:,ncon(ielem,4))
    !h_kind(:,7) = t_kind(:,ncon(ielem,6))
    !h_kind(:,8) = t_kind(:,ncon(ielem,8))
    !end function



    !   @params p0 input src point
    subroutine sing_int1(is,ielem,nodj,p0,amatrix,bmatrix) 

        use mvar_mod
        use mfunc_mod
        use hi_intg
        use green_funcs,only:Gcombo1_2
        use proj_cnst,only: ex,ey,xiqsi,xiqet

        implicit none

        integer,intent(in):: is,ielem,nodj
        real(8),intent(in):: p0(3) 
        real(8),intent(inout):: bmatrix(4,8),amatrix(4,8)

        integer ::inode,inodd,j,pwr_g

        real(8) :: src_lcl(2),src_glb(3),origin_offset(3)
        real(8) :: cnr_glb_mtx(3,8)
        real(8) :: passed_nrml(3,8)

        real(8) :: result0(8),result1(8),result2(8),prefix(3)
        real(8) ::  x0,y0,z0,si,eta,xp,yp,zp



        inode=ncon(ielem,nodj) ! corresponding node id
        inodd=ncond(ielem,nodj)!                 normal id
        amatrix(is,:)=0.0d0
        bmatrix(is,:)=0.0d0

        if(ncn(ielem).eq.8)  then 

            si =xiqsi(nodj) !get local coordinate for the src
            eta=xiqet(nodj)

        else if(ncn(ielem).eq.6)  then
            !si =xitsi(nodj)
            !eta=xitet(nodj)
        endif

        ! get src point in sysmetric mesh
        prefix=(/ex(is),ey(is),1.0d0/)
        xp=p0(1)
        yp=p0(2)
        zp=p0(3)

        x0=ex(is)*xp
        y0=ey(is)*yp      
        z0=zp

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

        passed_nrml(:,1) = dxyz(:,ncond(ielem,1))
        passed_nrml(:,2) = dxyz(:,ncond(ielem,3))
        passed_nrml(:,3) = dxyz(:,ncond(ielem,5))
        passed_nrml(:,4) = dxyz(:,ncond(ielem,7))
        passed_nrml(:,5) = dxyz(:,ncond(ielem,2))
        passed_nrml(:,6) = dxyz(:,ncond(ielem,4))
        passed_nrml(:,7) = dxyz(:,ncond(ielem,6))
        passed_nrml(:,8) = dxyz(:,ncond(ielem,8))
        call preset_src(si,eta,xyz(1:3,ncon(ielem,nodj)),origin_offset)

        ! add mirrored sink
        call norm_int(is,ielem,8,xp,yp,zp,amatrix(is,:),bmatrix(is,:),Gcombo1_2)
        !write(9013,'(2i6,8f12.6)') ielem,nodj,amatrix(is,:)

        call eval_singular_elem(cnr_glb_mtx,passed_nrml,result0,result1,result2)

        do j=1, ncn(ielem)
            amatrix(is,j) = amatrix(is,j)+result0(j)
            bmatrix(is,j) = bmatrix(is,j)+result1(j)
        end do
        write(9010,'(2i6,8f12.6)') ielem,nodj,result2
        !bmatrix=0.0d0


    end subroutine 
               

