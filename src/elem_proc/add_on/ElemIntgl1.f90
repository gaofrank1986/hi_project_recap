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
    subroutine norm_ele1(ielem,xp,yp,zp,amatrix,bmatrix)

        use mvar_mod,only:ncn,nsys
        use green_funcs,only:gcombo1

        implicit   none 
        integer,intent(in) :: ielem
        real(8),intent(in) ::  xp,yp,zp 
        real(8),intent(out) ::  bmatrix(4,8),amatrix(4,8)

        bmatrix=0.0d0
        amatrix=0.0d0

        do is=1,   nsys  
            call norm_int(is,ielem,ncn(ielem),xp,yp,zp,amatrix,bmatrix,gcombo1) 
        end do
        end subroutine

!
! ======================================================
!   Integration on an element with source in itself
!   or its mirror ones about any symmetrical axis
!
! ======================================================
!
       SUBROUTINE SING_ELE1(INODE,IELEM,NUMQUA,XP,YP,ZP,AMATRIX,BMATRIX)
       USE MVAR_MOD
       USE MFUNC_mod
!
      IMPLICIT   NONE  
!
      INTEGER I,J,IS,IELEM,INODE,NODNUM,ND,NP,NSAMB,NUMQUA
        REAL*8  XP,YP,ZP,XYZT(3,8),DXYZT(3,8)
        REAL*8 BMATRIX(4,8),AMATRIX(4,8)
!
        DO 5     I=1,  NCN(IELEM)
         
         IF(INODE.EQ.NCON(IELEM,I)) NODNUM=I ! get node num of inode

5       CONTINUE
!
          BMATRIX= 0.0d0     
          AMATRIX= 0.0d0  
!
        IF(NUMQUA.EQ.0)       THEN
            do IS=1,  NSYS
                IF(IS.EQ.1) THEN 
                    CALL SING_INT1(IS,IELEM,NODNUM,XP,YP,ZP,AMATRIX,BMATRIX)
                ELSE IF( IS.NE.1 ) THEN   
                    CALL NORM_INT(IS,IELEM,NCN(IELEM),XP,YP,ZP,AMATRIX,BMATRIX)
                    ! write(*,*) 'After Subroutine SGWP0_1'
                END IF
             end do
!
        ELSE IF(NUMQUA.EQ.2) THEN
        DO 200 IS=1,NSYS     
            IF(IS.EQ.1.OR.IS.EQ.2) THEN  
                CALL SING_INT1(IS,IELEM,NODNUM,XP,YP,ZP,AMATRIX,BMATRIX)
            ELSE IF(IS.EQ.3.OR.IS.EQ.4) THEN  
                CALL NORM_INT(IS,IELEM,NCN(IELEM),XP,YP,ZP,AMATRIX,BMATRIX) 
            END IF
            !
        200      CONTINUE  
!
        ELSE IF(NUMQUA.EQ.4) THEN
         DO 300  IS=1,  NSYS
          IF(IS.EQ.1.OR.IS.EQ.4) THEN   
            CALL SING_INT1(IS,IELEM,NODNUM,XP,YP,ZP,AMATRIX,BMATRIX)
          ELSE  IF(IS.EQ.2.OR.IS.EQ.3) THEN  
            CALL NORM_INT(IS,IELEM,NCN(IELEM),XP,YP,ZP,AMATRIX,BMATRIX) 
          END IF
300      CONTINUE
!
        ELSE IF(NUMQUA.EQ.5) THEN
         DO 400 IS=1, NSYS  
            CALL SING_INT1(IS,IELEM,NODNUM,XP,YP,ZP,AMATRIX,BMATRIX)
400      CONTINUE
        ENDIF
!
        RETURN
        END
                           
! ======================================================
!
! Integration on an element without source point
! 
! ======================================================
!                      
    subroutine norm_int(is,ielem,ncne,xp,yp,zp,amatrix,bmatrix,bie_called)
        use mvar_mod
        use proj_cnst,only : ex,ey
        implicit   none  

        integer,intent(in) :: is,ielem,ncne
        real(8),intent(in) ::  xp,yp,zp
        real(8),intent(inout) ::bmatrix(4,8),amatrix(4,8)
        
        interface
            subroutine bie_called(h,p0,p,gxf)
                real(8),intent(in) :: h,p0(3),p(3)
                real(8),intent(out) :: gxf(4)
            end subroutine
        end interface
        
        real(8) :: dgn,gxf(4),p0(3),p(3),prefix(3),np(3)
        integer :: n,nsamb,j

        amatrix(is,:)=0.0d0
        bmatrix(is,:)=0.0d0

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
                bmatrix(is,j)=bmatrix(is,j)+gxf(1)*samb(ielem,n,j)
                amatrix(is,j)=amatrix(is,j)+dgn*samb(ielem,n,j)
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
    subroutine sing_int1(is,ielem,nodj,xp,yp,zp,amatrix,bmatrix) 

        use mvar_mod
        use trvar_mod    
        use mfunc_mod
        use hi_intg
        use green_funcs,only:Gcombo1_2

        implicit none

        integer,intent(in):: is,ielem,nodj
        real(8),intent(in)::  xp,yp,zp
        real(8),intent(out):: bmatrix(4,8),amatrix(4,8)

        integer ::inode,inodd,j,pwr_g

        real(8) :: src_lcl(2),src_glb(3),origin_offset(3)
        real(8) :: cnr_glb_mtx(3,8)
        real(8) :: passed_nrml(3,8)

        real(8) :: result0(8),result1(8),result2(8)
        real(8) ::  x0,y0,z0,si,eta

        real(8) ::  xiqsi(8),xiqet(8)
        real(8) ::  ex(4),ey(4)


        DATA EX/  1.0d0,  1.0d0, -1.0d0, -1.0d0/                                                  
        DATA EY/  1.0d0, -1.0d0, -1.0d0,  1.0d0/

        ! ** XIQSI AND XIQET: NODE COORDINATES FOR QUADRILATERAL ELEMENTS
        data xiqsi/-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0,-1.0d0,-1.0d0/
        data xiqet/-1.0d0,-1.0d0,-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0/    


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
        call norm_int(is,ielem,8,xp,yp,zp,amatrix,bmatrix,Gcombo1_2)
        !write(9013,'(2i6,8f12.6)') ielem,nodj,amatrix(is,:)

        call eval_singular_elem(cnr_glb_mtx,passed_nrml,result0,result1,result2)

        do j=1, ncn(ielem)
            amatrix(is,j) = amatrix(is,j)+result0(j)
            bmatrix(is,j) = bmatrix(is,j)+result1(j)
        end do
        write(9010,'(2i6,8f12.6)') ielem,nodj,result2
        !bmatrix=0.0d0


    end subroutine 
               
        SUBROUTINE NORM_INT2(IS,IELEM,NCNE,XP,YP,ZP,AMATRIX,BMATRIX)
        USE MVAR_MOD
        IMPLICIT   NONE  
 
        INTEGER IS,IELEM,N,NSAMB,NCNE,J,IP
        REAL*8  XP,YP,ZP,EX(4),EY(4)
        REAL*8  X,X0,Y,Y0,Z,Z0 
        REAL*8  NX,NY,NZ,DGN
        REAL*8  DUM,WKX,PHi
        REAL*8  BMATRIX(4,8),AMATRIX(4,8),GXF(4)
!   
        DATA EX/  1.0d0,  1.0d0, -1.0d0, -1.0d0/                                                  
        DATA EY/  1.0d0, -1.0d0, -1.0d0,  1.0d0/
    write(9013,'(2i6,3f14.8)') ielem,ncne,xp,yp,zp
        write(9013,'(2i6,8f12.6)') ielem,8,amatrix(is,:)
!
!    PRINT *,' IN  NSWP0'
!
        NSAMB=16
        IF(NCNE.EQ.6)   NSAMB=4

        X0=EX(IS)*XP
        Y0=EY(IS)*YP
        Z0= ZP


        DO 100    N=1,   NSAMB     

         X =SAMBXY(IELEM,N,1)! guassian point info
         Y =SAMBXY(IELEM,N,2)
         Z =SAMBXY(IELEM,N,3)

       
       !mirrored sink
       !GXF=-GXF

        CALL DTGRN1(H,X,X0,Y,Y0,Z,Z0,GXF) 
!                      
          NX=EX(IS)*DSAMB(IELEM,N,1)
          NY=EY(IS)*DSAMB(IELEM,N,2)
          NZ=          DSAMB(IELEM,N,3)
          DGN=GXF(2)*Nx+GXF(3)*Ny+GXF(4)*Nz
        !write(9013,'(a5,3f16.6)') 'xyz',x,y,z
        !write(9013,'(a5,3f16.6)') 'xyz0',x0,y0,z0
        !write(9013,'(a5,3f16.6)') 'nxyz',nx,ny,nz
        !write(9013,'(a5,3f16.6)') 'gxf',gxf(2:4)
          
                         
        DO   J=1,   NCNE
          BMATRIX(IS,J)=BMATRIX(IS,J)+GXF(1)*SAMB(IELEM,N,J)!line integration
          AMATRIX(IS,J)=AMATRIX(IS,J)+DGN*SAMB(IELEM,N,J)!line integration
        ENDDO

100     CONTINUE
!
        RETURN
        END
