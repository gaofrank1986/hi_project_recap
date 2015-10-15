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
