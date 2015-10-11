module mesh
    
    implicit none
    
    real(8),public,allocatable :: XYZB(:,:),DXYZB(:,:)
    !xyzb => (3,node_id)  node data 
    !dxyzb => (3,nrml_id) derivative data
    integer,public,allocatable :: NCONB(:,:),NCONDB(:,:)   
    ! nconb => node list body mesh
    ! ncondb => normal list body
    real(8),allocatable :: XYZE(:,:,:),DXYZE(:,:,:),TXYZE(:,:,:)
    ! xyze =>  before combine, new full mesh for node (3,8,elem_id)
    ! dxyze => before combine, new full mesh for normal (3,8,elem,id)
    real(8),public,allocatable :: XYZTP(:,:),DXYZTP(:,:)
    ! xyztp => convsb, combined node mesh, xyztp(3,node_id)
    ! dxyztp = > combined normal mesh, dxyztp(3,nrml_id)
    integer,allocatable :: NCN(:),NCON(:,:),NCOND(:,:),IETYPE(:)
    ! ncn => elem_type ,namely, the node number in a elem
    ! ncon => combined node list,ncon(8,elem_id)
    ! ncond => combined nrml list,ncond(8,elem_id)
    ! ietype => flag show if a elem is free surface mesh or a body mesh

    integer,public,allocatable :: NNORMN(:),NNORMC(:)
    real(8),allocatable :: xyz(:,:),dxyz(:,:)
    integer :: nsys,nelem,nnode,nnoded,isys
    ! nsys => about symmetr
    ! nelem => elem number in combined mesh
    ! nnode => node number in combined mesh
    ! nnoded => nrml nmber in combined mesh
    integer :: nelemb,nnb,nnbd,nelemf,nnf
    ! nelemb => elem num in body
    ! nnb => node num  in body 
    ! nnbd => nrml num in body
    ! nelemf => elem num in fs
    ! nnf => node num in fs

    real(8),allocatable :: dampe(:,:),dampf(:),damptp(:)
    integer,allocatable :: nodele(:,:),nodnoe(:),nodelj(:,:),nodqua(:) 
contains

    subroutine read_mesh()
        implicit none

        INTEGER IFWKO,IPOL
      !INTEGER NTnum,IFLAG_T,IND,
      
        REAL*8  WL
      !REAL*8  FAMPR(6),FAMPI(6),FORCER(6),FORCEI(6)
      !REAL*8  PL_AMP(6),FORAMP
      !REAL*8  FCD_AMR,  FCD_AMI
    
        OPEN(2, FILE='INPUT/DATBDMS.txt',    STATUS='OLD') 
        OPEN(3, FILE='INPUT/DATWPMS.txt',    STATUS='OLD')

        !==================body mesh================================
        READ(2,*)   ISYS 
        READ(2,*)   NELEMB, NNB, NNBD, IPOL

        
        IF(ISYS.EQ.0) NSYS=1
        IF(ISYS.EQ.1) NSYS=2
        IF(ISYS.EQ.2) NSYS=4

        READ(3,*)   NELEMF
    
        ! =========================================

        NELEM=NELEMB+NELEMF

        ALLOCATE (NCONB(NELEMB,8),NCONDB(NELEMB,8))!body node/normal list
        ALLOCATE (XYZB(3,NNB),DXYZB(3,NNBD))    !body node/normal
        

        ALLOCATE(NCN(NELEM),NCON(NELEM,8),NCOND(NELEM,8),IETYPE(NELEM))
    
        allocate(NNORMN(8*NELEM) )

        ALLOCATE( XYZE(3,8,NELEM),DXYZE(3,8,NELEM))
        !allocate(DAMPE(8,NELEM))
        ALLOCATE( TXYZE(3,8,NELEM))
        ALLOCATE( XYZTP(3,8*NELEM),DXYZTP(3,8*NELEM))!,DAMPTP(8*NELEM))

        !ALLOCATE(SAMB(NELEM,16,0:8),SAMBXY(NELEM,16,3),
        !&         DSAMB(NELEM,16,6))


!        call MESHFS4()! Read in data on free surface mesh
!        call MESHBD(IPOL) ! Read in data on body mesh
!
!        close(2)
!        close(1)
!        close(3)
!
!        call convsb()
!        call prepare_mesh()
!
!        deallocate(xyzb,dxyzb,nconb,ncondb,nnormn,xyztp,dxyztp)

    end subroutine


    subroutine prepare_mesh()
        implicit none
        
        integer :: ind,node_max
        integer :: j,inode,l,ielem
        
        node_max = max(nnode,nnoded)
        
        allocate(xyz(3,node_max),dxyz(3,node_max))
        
        do ind = 1,nnode
            xyz(1:3,ind) = xyztp(1:3,ind)
        end do
       ! didn't assing damptp 
        do ind = 1,nnoded
            dxyz(1:3,ind) = dxyztp(1:3,ind)
        end do
        
        !deallocate(xyztp,dxyztp)

        !==============================================
        allocate(nodele(nnode,64),nodnoe(nnode),nodelj(nnode,64)&
                &,nodqua(nnode))
        
        
        DO 50 INODE=1, NNODE 
        L=0
        DO 40 IELEM=1,  NELEM
        DO 30 J=1,      NCN(IELEM)
        IF(INODE.EQ.NCON(IELEM,J)) THEN
        L=L+1
        NODELE(INODE,L)=IELEM
        NODELJ(INODE,L)=J
        ENDIF
30      CONTINUE
40      CONTINUE
        NODNOE(INODE)=L
!                          
        NODQUA(INODE)=0
        IF( NSYS .GE. 2) THEN
          IF( DABS(XYZ(2,INODE)).LT.1.0E-06 ) THEN
          NODQUA(INODE)=2
          END IF
        END IF
!
        IF( NSYS .EQ. 4) THEN
          IF( DABS(XYZ(1,INODE)).LT.1.0E-06.AND.&
             & DABS(XYZ(2,INODE)).LT.1.0E-06) THEN
           NODQUA(INODE)=5
          ELSE IF( DABS(XYZ(1,INODE)).LT.1.0E-06 ) THEN
           NODQUA(INODE)=4
          ENDIF
        END IF
!
50      CONTINUE
    end subroutine 
!       MESHFS4 + MESHBD 
!
!C *******************************************************************
!C *                                                                 *
!C *  Read in the data file of the mesh on the free surface          *
!C *                                                                 *
!C *                                                                 *
!C *******************************************************************
!C 
    subroutine meshfs42()

        implicit none

        integer :: IE,J,M
        
        XYZE(3,1:8, 1:NELEMF)=0.0d0
        
        DO 100 IE=1, NELEMF
            IETYPE(IE)=2
            READ(3, *)    M, NCN(IE)
            READ(3, *) (XYZE(1,J,IE), J=1, NCN(IE))
            READ(3, *) (XYZE(2,J,IE), J=1, NCN(IE))

      !dampe(:,:)=0.0d0

100    CONTINUE
        ! in dat file ,positive nrml is set as pointing into fulid field
      DXYZE(1, 1:8, 1:NELEMF)= 0.0d0
      DXYZE(2, 1:8, 1:NELEMF)= 0.0d0
      DXYZE(3, 1:8, 1:NELEMF)= 1.0d0 
!     DXYZE(3, 1:8, 1:NELEMF)=-1.0d0 

      RETURN
      END

!C ***************************************************************
!C *                                                             *
!C *  Generate nodal and element data on the body surface of     *
!C *  an arbitrary body                                          *
!C *                                                             *
!C ***************************************************************
!C 
        SUBROUTINE MESHBD_2(NCOR)

!     USE MVAR_MOD
        IMPLICIT   NONE  

!     INTEGER NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNBD,NNF,NNTCH


        INTEGER NCOR,IPL,IPOLAR(50)
        INTEGER I,IE,M,INODE,NCNN,K,KK
        REAL*8  A1,A2,R1,R2,Z1
        REAL*8  XOFSET(50),YOFSET(50),ZOFSET(50)
        real(8),parameter :: pi = 3.141592653589793
        !
        ! --------------------------------------------
        DO 5 I=1, NCOR
        READ(2,*)  M, IPOLAR(I), XOFSET(I),YOFSET(I),ZOFSET(I)
5       CONTINUE
!
        DO 10 INODE=1, NNB
        READ(2,*) M, IPL, R1, A1, Z1
        IF(IPOLAR(IPL).EQ.0) THEN
         XYZB(1,INODE)=R1+XOFSET(IPL)
         XYZB(2,INODE)=A1+YOFSET(IPL)
        ELSE IF(IPOLAR(IPL).EQ.1) THEN
         XYZB(1,INODE)=R1*DCOS(A1*PI/180.0D0)+XOFSET(IPL)
         XYZB(2,INODE)=R1*DSIN(A1*PI/180.0D0)+YOFSET(IPL)
        ENDIF
       XYZB(3,INODE)=Z1+ZOFSET(IPL)
10      CONTINUE
!       
        DO 11 INODE=1, NNBD
        READ(2,*) M,IPL,R2,A2,DXYZB(3,INODE)

        IF(IPOLAR(IPL).EQ.0) THEN
          DXYZB(1,INODE)= R2
          DXYZB(2,INODE)= A2
        ELSE IF(IPOLAR(IPL).NE.0) THEN
          DXYZB(1,INODE)=R2*DCOS( A2*PI/180.0D0 )
          DXYZB(2,INODE)=R2*DSIN( A2*PI/180.0D0 )
        ENDIF
11      CONTINUE

      
        DO 20 I=1, NELEMB
      IE=I+NELEMF
      IETYPE(IE)=1
        READ(2,*) M, NCN(IE)    
        READ(2,*) (NCONB(I,K), K=1, NCN(IE))
20    CONTINUE
!
     DO 30 I=1,  NELEMB
      IE=I+NELEMF
        READ(2,*)  M, NCNN
        READ(2,*) (NCONDB(I,K), K=1, NCN(IE))
30     CONTINUE
!
! ===============================================================
!    ×ª»»ÍøžñÊýŸÝžñÊœ£¬œ«Ë®ÏßÉÏœÚµãŽÓÎïÃæÉÏÏû³ý£¬¹éÈëµœË®ÃæÉÏ
!
       DO 100 I=1, NELEMB
      IE=I+NELEMF
      DO 80 K=1, NCN(IE)
         XYZE(1,K,IE)=XYZB(1, NCONB(I,K))
         XYZE(2,K,IE)=XYZB(2, NCONB(I,K))
         XYZE(3,K,IE)=XYZB(3, NCONB(I,K))
80      CONTINUE
!
100    CONTINUE
!
       DO 200 I=1, NELEMB
      IE=I+NELEMF
      DO 180 K=1, NCN(IE)
         DXYZE(1,K,IE)=DXYZB(1,NCONDB(I,K))
         DXYZE(2,K,IE)=DXYZB(2,NCONDB(I,K))
         DXYZE(3,K,IE)=DXYZB(3,NCONDB(I,K))
180     CONTINUE
!
200    CONTINUE
!
       END
  

    SUBROUTINE CONVSB2

       !USE MVAR_mod
     IMPLICIT NONE

       INTEGER IEB,IEL,IND,I,L,M,NND0,kk
     REAL(8) X,Y,Z,DX,DY,DZ,DR2,tmp1,tmp2,tmp3,tol

        !C ------------------------------------------------------------
        !C  NELEM : total number of elements 
        !C  NELEMF: total number of elements 
        !C  NELEMB: total number of elements 
        !C
        !C  NNODE : total number of nodes 
        !C  NNF   : total number of nodes on the free surface 
        !C  NNB   : total number of nodes on the body surface
        !C ------------------------------------------------------------
        !C    
        !C ** check nodes' situation 
          

        !WRITE(6,*) '  Inside CONVSB'
        !write(11,*) '  Inside CONVSB'
        !write(11,*) '  NELEMF=',NELEMF, '  NELEMB=',NELEMB
     
     
        !DO  IEL=1, NELEMF+NELEMB
         !   write(11,1000) IEL,NCN(IEL)
         !   write(11,1100) (XYZE(1,K,IEL),K=1,NCN(IEL))
         !   write(11,1100) (XYZE(2,K,IEL),K=1,NCN(IEL))
         !   write(11,1100) (XYZE(3,K,IEL),K=1,NCN(IEL))     
        !ENDDO
     
1100    FORMAT(8E15.6)
     
     
      M=0

      IF (NELEMF .EQ. 0) THEN
        NNF=0
      ELSE

        DO 200 IEL=1,  NELEMF
            DO 180 I=1, NCN(IEL)

                X=XYZE(1,I,IEL)
                Y=XYZE(2,I,IEL)
                Z=XYZE(3,I,IEL)
                !========check coincident node
                DO 150 L=1, M
                    tmp1 = DABS( X-XYZTP(1,L) )
                    tmp2 = DABS( Y-XYZTP(2,L) )
                    tmp3 = DABS( Z-XYZTP(3,L) )
                    tol = 1.0E-04
               IF((tmp1.LT.tol).AND.(tmp2.LT.tol).AND.(tmp3.lt.tol)) THEN
                 NCON(IEL, I)=L
                 GOTO 180
                    END IF
 150           CONTINUE

!C ** if the point (XPM, YPM) does not concide with any other
!C    points had been checked,  give a new code
!C
             M=M+1
             NCON(IEL, I)=M
             ! a new node matrix,xyztp 
             XYZTP(1,M)=X
             XYZTP(2,M)=Y
             XYZTP(3,M)=Z

180       CONTINUE
200     CONTINUE

       NNF=M     
       END IF    
!
! -----------------------------------------
!
         DO 300 IEB=1,  NELEMB
          IEL=IEB+NELEMF
           DO 280 I=1, NCN(IEL)


         X=XYZE(1,I,IEL)
         Y=XYZE(2,I,IEL)
         Z=XYZE(3,I,IEL)

             DO 250 L=1, M
               tmp1 = DABS( X-XYZTP(1,L) )
        tmp2 = DABS( Y-XYZTP(2,L) )
        tmp3 = DABS( Z-XYZTP(3,L) )
        tol = 1.0E-04
               IF((tmp1.LT.tol).AND.(tmp2.LT.tol).AND.(tmp3.lt.tol)) THEN
          NCON(IEL, I)=L
                 GOTO 280
               END IF
 250         CONTINUE
!C
!C ** if the point (XPM, YPM) does not concide with any other
!C    points had been checked,  give a new code
!C
             M=M+1
             NCON(IEL, I)=M
             XYZTP(1,M)=X
             XYZTP(2,M)=Y
             XYZTP(3,M)=Z

280       CONTINUE
300     CONTINUE

       NNODE=M
     NNB=NNODE-NNF
    ! write(11,*)  ' NNB=',NNB,' NNF=',NNF,' NNODE=',NNODE
     !write(11,*) 
!
! ===========================================================
!
!
       NNODED=NNODE
!
!  对每个节点先找出一个方向向量
!
       DO 350 IND=1,  NNODE
         DO 320 IEL=1,  NELEM
           DO 320 I=1,  NCN(IEL)

          IF(NCON(IEL, I) .EQ. IND) THEN
               DX=DXYZE(1,I,IEL)
             DY=DXYZE(2,I,IEL)
               DZ=DXYZE(3,I,IEL)
             NCOND(IEL, I)=IND
             DXYZTP(1,IND)=DX
             DXYZTP(2,IND)=DY
             DXYZTP(3,IND)=DZ
             NNORMN(IND)=IND
             GOTO 350
         ENDIF
320   CONTINUE
350  CONTINUE
!
!  对每个节点查找是否还有其他方向向量
!
     DO 500 IND=1,  NNODE

      NND0=NNODED

         DO 480 IEL=1,  NELEM
           DO 450 I=1,  NCN(IEL)

          IF(NCON(IEL, I)  .EQ. IND) THEN      
            DX=DXYZE(1,I,IEL)
            DY=DXYZE(2,I,IEL)
            DZ=DXYZE(3,I,IEL)
          
            DR2=((DX-DXYZTP(1,IND))**2+ &
                 (DY-DXYZTP(2,IND))**2+&
                 (DZ-DXYZTP(3,IND))**2)

           IF( DR2 .LT. 1.0E-6)  THEN 
             NCOND(IEL, I)=IND
             GOTO 450
           ENDIF
          
           DO KK=NND0+1,  NNODED    ! 检查是否与已定义的法向量相同
             DR2=(( DX-DXYZTP(1,KK))**2+&
                  ( DY-DXYZTP(2,KK))**2+&
                  ( DZ-DXYZTP(3,KK))**2)
           
            IF( DR2 .LT. 1.0E-6) THEN 
             NCOND(IEL, I)=KK
             GOTO 450
            ENDIF
              
           ENDDO
!
           
              NNODED=NNODED+1
            NCOND(IEL, I)=NNODED
            DXYZTP(1,NNODED)=DX
            DXYZTP(2,NNODED)=DY
            DXYZTP(3,NNODED)=DZ

            XYZTP(1,NNODED)=XYZTP(1,IND)
            XYZTP(2,NNODED)=XYZTP(2,IND)
            XYZTP(3,NNODED)=XYZTP(3,IND)
              NNORMN(NNODED)=IND

          ENDIF

!
450       CONTINUE
!       

480      CONTINUE
500     CONTINUE

1000   FORMAT(1X,8I6)
  
       END
   end module
