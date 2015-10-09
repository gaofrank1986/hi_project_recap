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
        SUBROUTINE NORM_ELE1(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)
        USE MVAR_MOD
        IMPLICIT   NONE 
      
        INTEGER IS,IP,N,ND,J,NP,IELEM
      
        REAL*8  XP,YP,ZP 
!       REAL*8  X,Y,Z,X0,Y0,Z0,XSB,YSB,ZSB       

        REAL*8  BMATRIX(4,8),AMATRIX(4,8)
                 
          BMATRIX=0.0D0
          AMATRIX=0.0D0

        do IS=1,   NSYS  
          CALL NORM_INT1(IS,IELEM,NCN(IELEM),XP,YP,ZP,AMATRIX,BMATRIX) 
          ! ncn(?) is element type
        end do
!
        END           

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
            CALL NORM_INT1(IS,IELEM,NCN(IELEM),XP,YP,ZP,AMATRIX,BMATRIX)
           ! write(*,*) 'After Subroutine SGWP0_1'
          END IF
         end do
!
        ELSE IF(NUMQUA.EQ.2) THEN
         DO 200 IS=1,NSYS     
          IF(IS.EQ.1.OR.IS.EQ.2) THEN  
            CALL SING_INT1(IS,IELEM,NODNUM,XP,YP,ZP,AMATRIX,BMATRIX)
          ELSE IF(IS.EQ.3.OR.IS.EQ.4) THEN  
            CALL NORM_INT1(IS,IELEM,NCN(IELEM),XP,YP,ZP,AMATRIX,BMATRIX) 
          END IF
!
200      CONTINUE  
!
        ELSE IF(NUMQUA.EQ.4) THEN
         DO 300  IS=1,  NSYS
          IF(IS.EQ.1.OR.IS.EQ.4) THEN   
            CALL SING_INT1(IS,IELEM,NODNUM,XP,YP,ZP,AMATRIX,BMATRIX)
          ELSE  IF(IS.EQ.2.OR.IS.EQ.3) THEN  
            CALL NORM_INT1(IS,IELEM,NCN(IELEM),XP,YP,ZP,AMATRIX,BMATRIX) 
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
        SUBROUTINE NORM_INT1(IS,IELEM,NCNE,XP,YP,ZP,AMATRIX,BMATRIX)
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
       
        CALL DTGRN(H,X,X0,Y,Y0,Z,Z0,GXF) 
!                      
          NX=EX(IS)*DSAMB(IELEM,N,1)
          NY=EY(IS)*DSAMB(IELEM,N,2)
          NZ=          DSAMB(IELEM,N,3)
          DGN=GXF(2)*Nx+GXF(3)*Ny+GXF(4)*Nz
                         
        DO   J=1,   NCNE
          BMATRIX(IS,J)=BMATRIX(IS,J)+GXF(1)*SAMB(IELEM,N,J)!line integration
          AMATRIX(IS,J)=AMATRIX(IS,J)+DGN*SAMB(IELEM,N,J)!line integration
        ENDDO

100     CONTINUE
!
        RETURN
        END
!
!
! *************************************************************
! *                                                           *
! *  The source point is in the mesh, at the node NODJ        *
! *                                                           *
! *************************************************************
!
        SUBROUTINE SING_INT1(IS,IELEM,NODJ,XP,YP,ZP,AMATRIX,BMATRIX) 
! 
        USE MVAR_MOD
        USE TRVar_mod    
        USE MFUNC_mod
        use hi_intg

        implicit none

        integer,intent(in):: is,ielem,nodj
        real*8,intent(in)::  xp,yp,zp
        real*8,intent(out):: bmatrix(4,8),amatrix(4,8)
        integer :: nf,ndim

        INTEGER N,J,IP   
        Integer Loop1,Loop2,I,NSAMB
        INTEGER LI,LJ,LK,INODE,INODD
!

        real(8) :: src_lcl(2),src_glb(3)
        real(8) :: cnr_glb_mtx(3,8)
        real(8) ::  SF_src(8),DSF_src(2,8),DDSF_src(3,8)
        !shape function ,derivative, double derivatives based on local src point
        real(8) ::  SF_iter(8),DSF_iter(2,8),DDSF_iter(3,8),DET1,DET2,DET3,DET 
        !shape function,..... based on iteration point



        REAL*8  EX(4,4),EY(4,4)
        REAL*8  X,Y,Z,X0,Y0,Z0,DX,DY,DZ,R,GXF(4)!,GXF0(4),GXF1(4)       
        REAL*8  XSB,YSB,ZSB,DUM,WKX,PHI  
        REAL*8  Nx,Ny,Nz,DAREA
        REAL*8  DPOX,DPOY,DPOZ,DGN,DGN0,DGN1
        REAL*8  XIQSI(8),XIQET(8),XITSI(6),XITET(6),SI,ETA,DUMX
        real*8  Xiq(8),Wiq(8)!,Xit(7),EtaT(7),WIT(7)
        REAL*8  SISM,ETASM
        REAL*8  XJ(2,3),XJSM(2,3),XXJ(3,3),XJp(2,3)

        REAL*8  GRN0J(8),GRN1J(8)

        integer :: debug_flag = 1,pwr_g
        integer :: debug_file_id = 108
        real(8) :: result0(8),ctr_glb(3)

        real(8) :: XIT(7),ETAT(7),WIT(9)
        real(8) :: PLO,CSST,SNST,SF,DSF,N0,N1C,N1S,JK0
        real(8) :: JK1C,JK1S,F1,F2,F,TOT,TOTJ,GXF0,GXF1

      DATA EX/  1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
                1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0/

      DATA EY/  1.0d0, -1.0d0, -1.0d0,  1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
                1.0d0, -1.0d0, -1.0d0,  1.0d0/

        ! ** XIQSI AND XIQET: NODE COORDINATES FOR QUADRILATERAL ELEMENTS
        DATA XIQSI/-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0,-1.0d0,-1.0d0/
        DATA XIQET/-1.0d0,-1.0d0,-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0/    

    DATA XIQ/ 0.960289856497536D+00, 0.796666477413626D+00, &
              0.525532409916329D+00, 0.183434642495650D+00, &
             -0.183434642495650D+00,-0.525532409916329D+00, &
             -0.796666477413626D+00,-0.960289856497536D+00/
  
    DATA WIQ/ 0.101228536290376D+00, 0.222381034453374D+00, &
              0.313706645877887D+00, 0.362683783378362D+00, &
              0.362683783378362D+00, 0.313706645877887D+00, &
              0.222381034453374D+00, 0.101228536290376D+00/     
!    ============================================    
      inode=ncon(ielem,nodj) ! corresponding node id
      inodd=ncond(ielem,nodj)!                 normal id


        write(110,*) '   ELEMT ID =',ielem
        write(110,*) '   NODE ID =',INODE

        if(ncn(ielem).eq.8)  then 

            si =xiqsi(nodj) !get local coordinate for the src
            eta=xiqet(nodj)

            call spfunc8_1(si,eta,sf_src,dsf_src,ddsf_src) 
            ! sf_src,dsf_src,ddsf_src is based on si,eta which is one node of the element

        else if(ncn(ielem).eq.6)  then
!           si =xitsi(nodj)
!           eta=xitet(nodj)
!           call spfunc6_1(si,eta,sf_src,dsf_src,ddsf_src) 
        endif

         
        write(debug_file_id,*) '   SI=',SI,'   ETA=',ETA

            X0=EX(IS,1)*XP
            Y0=EY(IS,1)*YP      
            Z0=ZP

            ! if mesh is created by using symmetrical information
            ! basically, we get the same local layout
            ! however, the glb position of the src point is different
            ! this will lead to ...
            ! when the src information is passed to the integration 
            ! both src local pos and global pos is set
            ! they has to match each other, the relationship cannot be reflected here
            ! since both NODJ,and XP-ZP are input information

            src_glb(1) = xp
            src_glb(2) = yp
            src_glb(3) = zp
            src_lcl(1) = si
            src_lcl(2) = eta

            ctr_glb(1) = 0
            ctr_glb(2) = 0
            ctr_glb(3) = 0
        write(110, *)  '  src_ctr_glb    ctr_glb=',ctr_glb

        pwr_g = ncn(ielem)/2+(ncn(ielem)/9)*2 
        !========switch cnr_glb_mtx order for hi module
        cnr_glb_mtx(:,1) = xyz(:,ncon(ielem,1))
        cnr_glb_mtx(:,2) = xyz(:,ncon(ielem,3))
        cnr_glb_mtx(:,3) = xyz(:,ncon(ielem,5))
        cnr_glb_mtx(:,4) = xyz(:,ncon(ielem,7))
        cnr_glb_mtx(:,5) = xyz(:,ncon(ielem,2))
        cnr_glb_mtx(:,6) = xyz(:,ncon(ielem,4))
        cnr_glb_mtx(:,7) = xyz(:,ncon(ielem,6))
        cnr_glb_mtx(:,8) = xyz(:,ncon(ielem,8))

        call set_npwg(pwr_g)
        call set_src_preset(si,eta,xyz(1:3,ncon(ielem,nodj)),ctr_glb)

        nf=8;ndim=3

        call eval_singular_elem(cnr_glb_mtx,result0)

        write (12,*) ielem,NODj,is,result0


        NSAMB=0
!
! **** FOR QUADRILATERIAL ELEMENTS  ****************
!
       write(debug_file_id,*) ' Surface Integration'

        GRN0J(:)=0.0d0  
        GRN1J(:)=0.0d0

    IF(NCN(IELEM).EQ.8)  THEN 

        DO 500 LOOP1=1, 8 
            DO 500 LOOP2=1, 8      
                NSAMB=NSAMB+1     
!         write(debug_file_id,*)  ' NSAMB=',NSAMB                    

            SISM=XIQ(loop1)
            ETASM=XIQ(loop2)  
!             PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
!             CSST=(SISM-SI)/PLO 
!             SNST=(ETASM-ETA)/PLO 
 
            CALL SPFUNC8_1(SISM,ETASM,SF_iter,DSF_iter,DDSF_iter)

            X=0.0D0
            Y=0.0D0
            Z=0.0D0
            NX=0.0D0
            NY=0.0D0
            NZ=0.0D0

            DO  LK=1,  NCN(IELEM)! 8 
                X  = X  + SF_iter(LK)*XYZ(1,NCON(IELEM,LK))  ! change ksi,eta to glb position  
                Y  = Y  + SF_iter(LK)*XYZ(2,NCON(IELEM,LK))
                Z  = Z  + SF_iter(LK)*XYZ(3,NCON(IELEM,LK))
                NX = NX + SF_iter(LK)*DXYZ(1,NCOND(IELEM,LK)) !get normal at eta,ksi  
                NY = NY + SF_iter(LK)*DXYZ(2,NCOND(IELEM,LK))
                NZ = NZ + SF_iter(LK)*DXYZ(3,NCOND(IELEM,LK))
            END DO
           
            DO  LI=1,2
                DO LJ=1,3
                    DUMX=0.0D0
                    DO    LK=1,NCN(IELEM)
                        DUMX=DUMX+DSF_iter(LI,LK)*XYZ(LJ,NCON(IELEM,LK))
                    Enddo
                    XJP(LI,LJ)=DUMX
                end do ! Lj
            end do !Li

            DET1=XJP(1,2)*XJP(2,3)-XJP(1,3)*XJP(2,2)    !  J1
            DET2=XJP(1,3)*XJP(2,1)-XJP(1,1)*XJP(2,3)    !  J2
            DET3=XJP(1,1)*XJP(2,2)-XJP(1,2)*XJP(2,1)    !  J3
            DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)
        !
            DAREA=DET*WIQ(LOOP1)*WIQ(LOOP2) !only for J

!             write(debug_file_id,*)  ' DAREA=',DAREA                    

!  ---------------------------------------------


          
        ! x0-z0, src_ctr_glb for green fucntion


!           CALL DTGRN0(H,X,X0,Y,Y0,Z,Z0,GXF0)
!         DGN0=GXF0(2)*NX+GXF0(3)*NY+GXF0(4)*NZ
!          write(debug_file_id,*)  ' GXF0=',GXF0(2),GXF0(3),GXF0(4)                    
!          write(debug_file_id,*)  ' DGN0=',DGN0                    
              
! !        
!           CALL DTGRN1(H,X,X0,Y,Y0,Z,Z0,GXF1) 
!         DGN1=GXF1(2)*NX+GXF1(3)*NY+GXF1(4)*NZ
! !         write(debug_file_id,*)  ' GXF1=',GXF1                    
! !         write(debug_file_id,*)  ' DGN1=',DGN1                    

        CALL DTGRN (H,X,X0,Y,Y0,Z,Z0,GXF) 
        DGN=GXF(2)*NX+GXF(3)*NY+GXF(4)*NZ
        ! DGN is derivative of Green fucntion over normal vector,
        ! partical G over partial n

        DO J=1, NCN(IELEM)
!             AMATRIX(IS,J)=AMATRIX(IS,J)+DGN*DAREA*SF_iter(J)
            ! difference only exists for using SF_iter(J) and Darea,whch is equivalent ot Jacobian matrix

!             GRN0J(J)=GRN0J(J)+DGN0*DAREA*SF_iter(J)
!             GRN1J(J)=GRN1J(J)+DGN1*DAREA*SF_iter(J) 

            BMATRIX(IS,J)=BMATRIX(IS,J)+GXF(1)*DAREA*SF_iter(J)
            ! B matrix is only G over sufrace

        ENDDO
     
500    CONTINUE
!

!
! **** FOR TRIANGULAR ELEMENTS **********************
! 
        ELSE IF(NCN(IELEM).EQ.6)  THEN
!
!       DO 800 NSAMB=1, 7
      
!        SISM=XIT(NSAMB)
!        ETASM=ETAT(NSAMB) 
!        PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
!        CSST=(SISM-SI)/PLO 
!        SNST=(ETASM-ETA)/PLO 
      
!        CALL SPFUNC6_1(SISM,ETASM,SF,DSF,DDSF)

!       X=0.0D0
!       Y=0.0D0
!       Z=0.0D0
!       NX=0.0D0
!       NY=0.0D0
!       NZ=0.0D0
!       DO  LK=1,  NCN(IELEM)
!         X=X+SF(LK)*XYZ(1,NCON(IELEM,LK))    
!         Y=Y+SF(LK)*XYZ(2,NCON(IELEM,LK))
!         Z=Z+SF(LK)*XYZ(3,NCON(IELEM,LK))
!         NX=NX+SF(LK)*DXYZ(1,NCOND(IELEM,LK))    
!         NY=NY+SF(LK)*DXYZ(2,NCOND(IELEM,LK))
!         NZ=NZ+SF(LK)*DXYZ(3,NCOND(IELEM,LK))
!       END DO
         
!           DO  630  LI=1,2
!           DO  630  LJ=1,3
!           DUMX=0.0D0
!           DO    LK=1,NCN(IELEM)
!          DUMX=DUMX+DSF(LI,LK)*XYZ(LJ,NCON(IELEM,LK))
!         enddo
!            XJP(LI,LJ)=DUMX
! 630      CONTINUE
! !
!           DET1=XJP(1,2)*XJP(2,3)-XJP(1,3)*XJP(2,2)    !  J1
!           DET2=XJP(1,3)*XJP(2,1)-XJP(1,1)*XJP(2,3)    !  J2
!           DET3=XJP(1,1)*XJP(2,2)-XJP(1,2)*XJP(2,1)    !  J3
!           DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)

!         DAREA=DET*WIT(NSAMB)
!  !  ---------------------------------------------
       
!         DO J=1, NCN(IELEM)
!         N0 =SF_src(J)
!         N1C=DSF_src(1,J)
!         N1S=DSF_src(2,J)   
       
!        CALL AREA_COEF(CSST,SNST,JK0,JK1C,JK1S,   &
!                        N0,N1C,N1S,XJ,XXJ,F1,F2)

!         F=F2/PLO/PLO+F1/PLO                       ! (C19)
!         TOT=F/PLO*WIT(NSAMB)
 
!         AMATRIX(IS,J)=AMATRIX(IS,J)-TOT
!         TOTJ(J)=TOTJ(J)-TOT
!        ENDDO

! !  ---------------------------------------------

!         X0=EX(IS,1)*XP
!         Y0=EY(IS,1)*YP        
!         Z0=ZP        
!           CALL DTGRN0(H,X,X0,Y,Y0,Z,Z0,GXF0)       
!         DGN0=GXF0(2)*NX+GXF0(3)*NY+GXF0(4)*NZ
     
!           CALL DTGRN1(H,X,X0,Y,Y0,Z,Z0,GXF1) 
!         DGN1=GXF1(2)*NX+GXF1(3)*NY+GXF1(4)*NZ        

!           CALL DTGRN (H,X,X0,Y,Y0,Z,Z0,GXF) 
!          DGN=GXF(2)*NX+GXF(3)*NY+GXF(4)*NZ
        
!         DO J=1, NCN(IELEM)
!          AMATRIX(IS,J)=AMATRIX(IS,J)+DGN*DAREA*SF(J)
! !        GRN0J(J)=GRN0J(J)+DGN0*DAREA*SF(J)
! !        GRN1J(J)=GRN1J(J)+DGN1*DAREA*SF(J)
!          BMATRIX(IS,J)=BMATRIX(IS,J)+GXF(1)*DAREA*SF(J)
!        ENDDO

700    CONTINUE     

800    CONTINUE 

        ENDIF
!
! ======================================================================
!



!
        DO J=1, NCN(IELEM)
            AMATRIX(IS,J) = result0(J)
        END DO
        write(110,*) '   sieppem result',result0
        write(110,*) 'sum of AMATRIX',SUM(AMATRIX(IS,1:8))
        END
               
