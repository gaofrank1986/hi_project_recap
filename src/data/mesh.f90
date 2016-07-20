!!-------------------------------------------------------------------------------
! DUTWAV
!-------------------------------------------------------------------------------
!  MODULE: mesh
!
!> @brief
!! 
!!
!! @author
!! DUT 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 7.20.2016 : add read damp info  
!
!------------------------------------------------------------------------------- 
module mesh

    use wave 

    implicit none
    ! @var : [xyzb] store body node with dimension(3,node_id) 
    ! @var : [dxyzb] store body nrml with dimension (3,nrml_id)
    real(8),private,allocatable :: xyzb(:,:),dxyzb(:,:)

    integer,private,allocatable :: nconb(:,:),ncondb(:,:)   
    ! nconb => node list body mesh
    ! ncondb => normal list body
    real(8),protected,allocatable :: XYZE(:,:,:),DXYZE(:,:,:),TXYZE(:,:,:)
    ! @var : [xyze]  before combine, new full mesh for node (3,8,elem_id)
    ! @var : [dxyze] before combine, new full mesh for normal (3,8,elem,id)
    real(8),private,allocatable :: XYZTP(:,:),DXYZTP(:,:)
    ! xyztp => convsb, combined node mesh, xyztp(3,node_id)
    ! dxyztp = > combined normal mesh, dxyztp(3,nrml_id)
    integer,protected,allocatable :: NCN(:),NCON(:,:),NCOND(:,:),IETYPE(:)
    ! ncn => elem_type ,namely, the node number in a elem
    ! ncon => combined node list,ncon(8,elem_id)
    ! ncond => combined nrml list,ncond(8,elem_id)
    ! ietype => flag show if a elem is free surface mesh or a body mesh

    integer,private,allocatable :: NNORMN(:),NNORMC(:)
    real(8),protected,allocatable :: xyz(:,:)
    real(8),public,allocatable ::dxyz(:,:)
    integer,protected :: nsys,nelem,nnode,nnoded,isys
    ! nsys => about symmetr
    ! nelem => elem number in combined mesh
    ! nnode => node number in combined mesh
    ! nnoded => nrml nmber in combined mesh
    integer,protected :: nelemb,nnb,nnbd,nelemf,nnf
    ! nelemb => elem num in body
    ! nnb => node num  in body 
    ! nnbd => nrml num in body
    ! nelemf => elem num in fs
    ! nnf => node num in fs

    real(8),protected,allocatable :: dampe(:,:),dampf(:),damptp(:)
    integer,protected,allocatable :: nodele(:,:),nodnoe(:),nodelj(:,:),nodqua(:) 
contains
    subroutine read_mesh()
        implicit none

        INTEGER IFWKO,IPOL
      !INTEGER NTnum,IFLAG_T,IND,
      
!        REAL*8  WL
      !REAL*8  FAMPR(6),FAMPI(6),FORCER(6),FORCEI(6)
      !REAL*8  PL_AMP(6),FORAMP
      !REAL*8  FCD_AMR,  FCD_AMI
    
        open(2, FILE='INPUT/DATBDMS.txt',    status='old') 
        open(3, FILE='INPUT/DATWPMS.txt',    status='old')
        open(11, FILE='output/outmesh.txt',    status='unknown')

        !==================body mesh================================
        READ(2,*)   ISYS 
        READ(2,*)   NELEMB, NNB, NNBD, IPOL

        
        IF(ISYS.EQ.0) NSYS=1
        IF(ISYS.EQ.1) NSYS=2
        IF(ISYS.EQ.2) NSYS=4

        READ(3,*)   NELEMF
    
        ! =========================================

        NELEM=NELEMB+NELEMF

        allocate (nconb(nelemb,8),ncondb(nelemb,8))!body node/normal list
        allocate (xyzb(3,nnb),dxyzb(3,nnbd))    !body node/normal
        allocate(ncn(nelem),ncon(nelem,8),ncond(nelem,8),ietype(nelem))

        allocate(nnormn(8*nelem) )

        allocate( xyze(3,8,nelem),dxyze(3,8,nelem))
        allocate( txyze(3,8,nelem))
        allocate(xyztp(3,8*nelem),dxyztp(3,8*nelem))
        allocate(dampe(8,nelem),damptp(8*nelem))
        !ALLOCATE(SAMB(NELEM,16,0:8),SAMBXY(NELEM,16,3),
        !&         DSAMB(NELEM,16,6))


        call MESHFS4_2(0)! Read in data on free surface mesh
        call MESHBD_2(IPOL) ! Read in data on body mesh

        close(2)
        close(1)
        close(3)

        call convsb_2()
!        call prepare_mesh()
        call pre_mesh_2()
!
!        deallocate(xyzb,dxyzb,nconb,ncondb,nnormn,xyztp,dxyztp)

    end subroutine


subroutine pre_mesh_2()
        implicit none
        integer :: ndmax,ind
          NDMAX=MAX(NNODE,NNODED)
       ALLOCATE(NODELE(NNODE,64),NODNOE(NNODE),NODELJ(NNODE,64),&
     &         NODQUA(NNODE),NNORMC(NNODED),&
     &           XYZ(3,NDMAX),DXYZ(6,NDMAX),DAMPF(NNF))
!
! ---------------------------------------------------
!
!         WRITE(10,*) '  IND       X          Y          Z        DAMPing'
         DO IND=1, NNF
          DAMPF(IND)=DAMPTP(IND)
          XYZ(1,IND)=XYZTP(1,IND)
          XYZ(2,IND)=XYZTP(2,IND)
          XYZ(3,IND)=XYZTP(3,IND)
          !WRITE(10,111) IND, XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),DAMPF(IND)
         END DO
!                
         DO IND=NNF+1, NNODE
          XYZ(1,IND)=XYZTP(1,IND)
          XYZ(2,IND)=XYZTP(2,IND)
          XYZ(3,IND)=XYZTP(3,IND)
          !WRITE(10,111) IND, XYZ(1,IND),XYZ(2,IND),XYZ(3,IND)
         END DO
!
         DO IND=1, NNODED
          DXYZ(1,IND)=DXYZTP(1,IND)
          DXYZ(2,IND)=DXYZTP(2,IND)
          DXYZ(3,IND)=DXYZTP(3,IND)
          NNORMC(IND)=NNORMN(IND)
         ! WRITE(10,111) IND, DXYZ(1,IND),DXYZ(2,IND),DXYZ(3,IND)
         END DO
!
!
         DEALLOCATE(XYZTP,DXYZTP,NNORMN)
!
         DAMPF(:)=W1*DAMPF(:)
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
    ! @func : read in data file of free surface mesh
    ! @param: [flag] 0 if no damping info, 1 if damp info presents

    subroutine meshfs4_2(flag)

        implicit none

        integer :: IE,J,M,flag

        XYZE(3,1:8, 1:NELEMF)=0.0d0

        DO 100 IE=1, NELEMF
            IETYPE(IE)=2
            READ(3, *)    M, NCN(IE)
            READ(3, *) (XYZE(1,J,IE), J=1, NCN(IE))
            READ(3, *) (XYZE(2,J,IE), J=1, NCN(IE))
        if (flag.eq.0) then
            dampe(:,:)=0.0d0
        else if(flag.eq.1) then
            read(3, *) (dampe(j,ie), j=1, ncn(ie))
        end if


        100    CONTINUE
        ! in dat file ,positive nrml is set as pointing into fulid field
        DXYZE(1, 1:8, 1:NELEMF)= 0.0d0
        DXYZE(2, 1:8, 1:NELEMF)= 0.0d0
        DXYZE(3, 1:8, 1:NELEMF)= 1.0d0 
        !     DXYZE(3, 1:8, 1:NELEMF)=-1.0d0 


    end subroutine

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
!        real(8),parameter :: pi = 3.141592653589793
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
  

    SUBROUTINE CONVSB_2

       !USE MVAR_mod
     IMPLICIT NONE

       INTEGER IEB,IEL,IND,I,L,M,NND0,kk,N
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
         write(11,*) Isys
        write(11,*)  NELEM, NNODE,NNODED,1
        write(11,*) ' 1  0  0.0  0.0  0.0'

	 DO N=1, NNODE
	  write(11,1010) N, 1, XYZTP(1,N),XYZTP(2,N),XYZTP(3,N)
	 ENDDO
	 DO N=1, NNODED
	  write(11,1010) N, 1, DXYZTP(1,N),DXYZTP(2,N),DXYZTP(3,N)
	 ENDDO
!
	 DO IEL=1, NELEM
	  write(11,1001) IEL, NCN(IEL)
        write(11,1005) (NCON(IEL, I), I=1, NCN(IEL))
	 ENDDO
!
	 DO IEL=1, NELEM
	  write(11,1001) IEL, NCN(IEL)
        WRITE(11,1005) (NCOND(IEL, I), I=1, NCN(IEL))
	 ENDDO
1001   FORMAT(1X,2I6,3F14.6)
1010   FORMAT(1X,I6,I4,3F14.6)
1005   FORMAT(8(1X,I6))
1020   FORMAT(1X,I6,I4,6F14.6)
       END

    subroutine comp_link(ielem,inode,ii) 
        !use mvar_mod
        ! return how many nodes in this element 
        implicit none
        integer,intent(in) :: inode,ielem
        integer,intent(out) :: ii

        integer :: i

        ii = 0 
        do i=1, nodnoe(inode)!!list of linked element by inode (how manys times inode appear in element)
        if(ielem .eq. nodele(inode,i)) then!
          ii=ii+1
        endif
        enddo
    end subroutine 

    logical function is_connected(ielem,inode) 
        !use mvar_mod
        ! return how many nodes in this element 
        implicit none
        integer,intent(in) :: inode,ielem
        !logical :: is_connected

        integer :: i,ii

        ii = 0 
        do i=1, nodnoe(inode)!!list of linked element by inode (how manys times inode appear in element)
        if(ielem .eq. nodele(inode,i)) then!
          ii=ii+1
        endif
        enddo
        if (ii==0) then
            is_connected=.false.
        else
            is_connected=.true.
        endif
    end function 

    subroutine topology_analysis()
        !use mvar_mod
        implicit none
        integer :: inode,ielem,j,l
        do 50 inode=1, nnode 
        l=0
        do 40 ielem=1,  nelem
        do 30 j=1,      ncn(ielem)
        if(inode.eq.ncon(ielem,j)) then
        l=l+1
        nodele(inode,l)=ielem! elem num linked to inode
        nodelj(inode,l)=j!in node-linked-element, inode appear as j-th node 
        endif
30      continue
40      continue
              nodnoe(inode)=l !total number of links
!        below related to symmetry
        nodqua(inode)=0
        if( nsys .ge. 2) then
          if( abs(xyz(2,inode)).lt.1.0e-06 ) then
          nodqua(inode)=2
          end if
        end if
!
        if( nsys .eq. 4) then
          if( abs(xyz(1,inode)).lt.1.0e-06.and.&
     &        abs(xyz(2,inode)).lt.1.0e-06) then
           nodqua(inode)=5
          else if( abs(xyz(1,inode)).lt.1.0e-06 ) then
           nodqua(inode)=4
          endif
        end if
!
50      continue
        PRINT *,"topology analysis finished"
    end subroutine
   end module
