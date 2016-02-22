module mesh
    use wave 
    implicit none
    
    real(8),public,allocatable :: xyzb(:,:),dxyzb(:,:)
    !xyzb => (3,node_id)  node data 
    !dxyzb => (3,nrml_id) derivative data
    integer,public,allocatable :: nconb(:,:),ncondb(:,:)   
    ! nconb => node list body mesh
    ! ncondb => normal list body
    real(8),allocatable :: xyze(:,:,:),dxyze(:,:,:),txyze(:,:,:)
    ! xyze =>  before combine, new full mesh for node (3,8,elem_id)
    ! dxyze => before combine, new full mesh for normal (3,8,elem,id)
    real(8),public,allocatable :: xyztp(:,:),dxyztp(:,:)
    ! xyztp => convsb, combined node mesh, xyztp(3,node_id)
    ! dxyztp = > combined normal mesh, dxyztp(3,nrml_id)
    integer,allocatable :: ncn(:),ncon(:,:),ncond(:,:),ietype(:)
    ! ncn => elem_type ,namely, the node number in a elem
    ! ncon => combined node list,ncon(8,elem_id)
    ! ncond => combined nrml list,ncond(8,elem_id)
    ! ietype => flag show if a elem is free surface mesh or a body mesh

    integer,public,allocatable :: nnormn(:),nnormc(:)
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
        include './include/prepare_mesh_2.f90'
    subroutine read_mesh()
        implicit none

        integer ifwko,ipol
      !integer ntnum,iflag_t,ind,
      
!        real*8  wl
      !real*8  fampr(6),fampi(6),forcer(6),forcei(6)
      !real*8  pl_amp(6),foramp
      !real*8  fcd_amr,  fcd_ami
    
        open(2, file='input/datbdms.txt',    status='old') 
        open(3, file='input/datwpms.txt',    status='old')
        open(11, file='output/outmesh.txt',    status='unknown')

        !==================body mesh================================
        read(2,*)   isys 
        read(2,*)   nelemb, nnb, nnbd, ipol

        
        if(isys.eq.0) nsys=1
        if(isys.eq.1) nsys=2
        if(isys.eq.2) nsys=4

        read(3,*)   nelemf
    
        ! =========================================

        nelem=nelemb+nelemf

        allocate (nconb(nelemb,8),ncondb(nelemb,8))!body node/normal list
        allocate (xyzb(3,nnb),dxyzb(3,nnbd))    !body node/normal
        

        allocate(ncn(nelem),ncon(nelem,8),ncond(nelem,8),ietype(nelem))
    
        allocate(nnormn(8*nelem) )

        allocate( xyze(3,8,nelem),dxyze(3,8,nelem))
        allocate( txyze(3,8,nelem))
        allocate(xyztp(3,8*nelem),dxyztp(3,8*nelem))
        allocate(dampe(8,nelem),damptp(8*nelem))
        !allocate(samb(nelem,16,0:8),sambxy(nelem,16,3),
        !&         dsamb(nelem,16,6))


        call meshfs4_2()! read in data on free surface mesh
        call meshbd_2(ipol) ! read in data on body mesh

        close(2)
        close(3)

        call convsb_2()
!        call prepare_mesh()
        call pre_mesh_2()
!
!        deallocate(xyzb,dxyzb,nconb,ncondb,nnormn,xyztp,dxyztp)

    end subroutine


    subroutine prepare_mesh()
        implicit none
        
        integer :: ind,node_max
        integer :: j,inode,l,ielem
        
        node_max = max(nnode,nnoded)
        
        allocate(xyz(3,node_max),dxyz(3,node_max),nnormc(nnoded))
        
        do ind = 1,nnode
            xyz(1:3,ind) = xyztp(1:3,ind)
            if (ind <= nnf) then
                dampf(ind) = damptp(ind)
            end if
        end do
       ! didn't assing damptp 
        do ind = 1,nnoded
            dxyz(1:3,ind) = dxyztp(1:3,ind)
                nnormc(ind)=nnormn(ind)

        end do
        dampf(:) = w1*dampf(:) 
        !deallocate(xyztp,dxyztp)

        !==============================================
        allocate(nodele(nnode,64),nodnoe(nnode),nodelj(nnode,64)&
                &,nodqua(nnode))
        
        
!        do 50 inode=1, nnode 
!        l=0
!        do 40 ielem=1,  nelem
!        do 30 j=1,      ncn(ielem)
!        if(inode.eq.ncon(ielem,j)) then
!        l=l+1
!        nodele(inode,l)=ielem
!        nodelj(inode,l)=j
!        endif
!30      continue
!40      continue
!        nodnoe(inode)=l
!!                          
!        nodqua(inode)=0
!        if( nsys .ge. 2) then
!          if( dabs(xyz(2,inode)).lt.1.0e-06 ) then
!          nodqua(inode)=2
!          end if
!        end if
!!
!        if( nsys .eq. 4) then
!          if( dabs(xyz(1,inode)).lt.1.0e-06.and.&
!             & dabs(xyz(2,inode)).lt.1.0e-06) then
!           nodqua(inode)=5
!          else if( dabs(xyz(1,inode)).lt.1.0e-06 ) then
!           nodqua(inode)=4
!          endif
!        end if
!!
!50      continue
    end subroutine 
!       meshfs4 + meshbd 
!
!c *******************************************************************
!c *                                                                 *
!c *  read in the data file of the mesh on the free surface          *
!c *                                                                 *
!c *                                                                 *
!c *******************************************************************
!c 
    subroutine meshfs4_2()

        implicit none

        integer :: ie,j,m
        
        xyze(3,1:8, 1:nelemf)=0.0d0
        
        do 100 ie=1, nelemf
            IETYPE(IE)=2
            READ(3, *)    M, NCN(IE)
            READ(3, *) (XYZE(1,J,IE), J=1, NCN(IE))
            READ(3, *) (XYZE(2,J,IE), J=1, NCN(IE))
            READ(3, *) (DAMPE(J,IE), J=1, NCN(IE))

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
        close(11)
1001   FORMAT(1X,2I6,3F14.6)
1010   FORMAT(1X,I6,I4,3F14.6)
1005   FORMAT(8(1X,I6))
1020   FORMAT(1X,I6,I4,6F14.6)
       END
   end module
