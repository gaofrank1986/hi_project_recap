!C  TASSB0
!C *********************************************************
!C *                                                       *
!C * Calculate the element contribution and assembly the   *
!C * coefficients of the corresponding system of equationn *
!C *                                                       *
!C *********************************************************
    include './add_on/common_block.f90'
    subroutine comp_link(ielem,inode,ii) 
        use MVAR_MOD
        implicit none
        integer,intent(in) :: inode,ielem
        integer,intent(out) :: ii

        integer :: i

        ii = 0 
        do i=1, nodnoe(inode)
        if(ielem .eq. nodele(inode,i)) then
          ii=ii+1
        endif
        enddo
    end subroutine 

    subroutine tassb0
        use mvar_mod
        use pvar_mod
        use mfunc_mod
        use sebsm_mod

        implicit   none  
        integer  inode,ielem,j,jnode,ind,indd,ip
        integer ::    i,ii,is,l
        real*8  xp,yp,zp,dpox,dpoy,dpoz,phi
        real*8  rsn(4,4)
        real*8  bmatrix(4,8),amatrix(4,8),bmat(4)

        real*8  s_angle
        real*8  dsign
        real(8) :: fterm_coef(0:3,4)


        DATA RSN /1.,  1.,  1.,  1., &
     &            1., -1.,  1., -1., &
     &            1.,  1., -1., -1., &
     &            1., -1., -1.,  1./ 
!
!
!  ----------------------------------------------------
          WRITE(10, *)   ' IN TASSB0 '
          DSDT(:)=0.0
!                 
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
     &        DABS(XYZ(2,INODE)).LT.1.0E-06) THEN
           NODQUA(INODE)=5
          ELSE IF( DABS(XYZ(1,INODE)).LT.1.0E-06 ) THEN
           NODQUA(INODE)=4
          ENDIF
        END IF
!
50      CONTINUE

        amata = 0.0d0
        bmata = 0.0d0

!        open (402,file = 'xyz_matrix.txt',status = 'unknown')
!        open (403,file = 'dxyz_matrix.txt',status = 'unknown')
!        open (404,file = 'dxyze_matrix.txt',status = 'unknown')
!        open (405,file = 'ncn.txt',status = 'unknown')
!        open (406,file = 'ncon.txt',status = 'unknown')
!        
!        do i = 1,nnode
!            write(402,54) xyz(:,i)
!        end do
!        do i = 1,nnoded
!            write(403,54) dxyz(:,i)
!        end do
!        do i = 1,nelem
!            write(404,52) dxyze(1,:,nelem)
!            write(404,52) dxyze(2,:,nelem)
!            write(404,52) dxyze(3,:,nelem)
!        end do
!        do i = 1,nelem
!            write(405,53) i,ncn(i)
!            write(406,55) i,ncon(:,i)
!        end do
!51    FORMAT(6I5)
!52      format(8f10.6)
!53      format(2i4)
!54      format(3f10.6)
!55      format(9i6)
! =======================================================================

        WRITE(9, *) '   INODE      XP       YP       ZP       S_ANGLE'
        WRITE(10,*) '  INODE      ANGLE         A3   ',&
     &               '       C31           C32           C33'

        do  500   inode=1,  nnf

            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode)

            fterm_coef = 0
            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
             &                        h,xyz,dxyze,s_angle)    

            s_angle=1.0d0-s_angle

            write(9,102)  inode, xp, yp, zp, s_angle
            write(*,102)  inode, xp, yp, zp, s_angle
            write(101,102)  inode, xp, yp, zp, s_angle

            102  format(i6,3f12.4,f15.6) 

            angle(inode)=s_angle
            amata(inode,inode,1:nsys)= angle(inode)
            !  ---------------------------
            !  Integration on the free surface

            do   ielem=1,  nelemf

                write(101,*)
                write(101,*) ' ielem=',ielem
                write(101,*)  ' xp,yp,zp=',xp,yp,zp

                call comp_link(ielem,inode,ii)
                if (ii .eq. 0)   then 
                    call norm_ele1(ielem,xp,yp,zp,amatrix,bmatrix)
                else if (ii .ne. 0)   then 
                    call sing_ele1(inode,ielem,nodqua(inode),xp,yp,zp,&
                        &                   amatrix,bmatrix)
                end if 
                call common_block(0,0,ielem,inode,amatrix,bmatrix,fterm_coef)
            end do

            !  Integration on the body surface

            do    ielem=nelemf+1,  nelem
                write(101,*)
                write(101,*)  ' ielem=',ielem
                write(101,*)  ' xp,yp,zp=',xp,yp,zp

                call comp_link(ielem,inode,ii)
                if (ii .eq. 0)   then 
                    call norm_ele1(ielem,xp,yp,zp,amatrix,bmatrix)
                else if (ii .ne. 0)   then 
                    call sing_ele1(inode,ielem,nodqua(inode),xp,yp,zp,&
                     &                   amatrix,bmatrix)
                end if
                call common_block(1,0,ielem,inode,amatrix,bmatrix,fterm_coef)

            end do

            fra3(inode)=fterm_coef(0,1)!
            frc31(inode)=fterm_coef(1,1)-fterm_coef(0,1)*xp!
            frc32(inode)=fterm_coef(2,1)-fterm_coef(0,1)*yp!
            frc33(inode)=fterm_coef(3,1)-fterm_coef(0,1)*zp!

            !         write(10,620) inode,angle(inode),fra3(inode),&
            !     &                 frc31(inode),frc32(inode),frc33(inode)
            write (402,499) fterm_coef(0:3,1)
            499 format(4f10.6)
             phi=poxy(xp,yp,zp) ! this in known since src on free surface
             call dinp(xp,yp,zp,dpox,dpoy,dpoz)! this is unkown, so is this a mistake
             bmata(inode,1)=bmata(inode,1)-fra3(inode)*phi-&
                 &     frc31(inode)*dpox-frc32(inode)*dpoy
        500     continue
!
! =======================================================================
!    Source point is on the body surface
!
        do  1000   inode=nnf+1, nnode   

            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode) 

            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
             &                    h,xyz,dxyze,s_angle) 

            s_angle=1.0d0-s_angle

            write(9,102)  inode, xp, yp, zp, s_angle
            write(*,102)  inode, xp, yp, zp, s_angle

            angle(inode)=s_angle

            amata(inode,inode,1:nsys)= angle(inode)

! ------------------------------
! Intergration on the free surface
! 
        do   ielem=1,  nelemf

            ii=0   
            call norm_ele0(ielem,xp,yp,zp,amatrix,bmatrix)
            call common_block(0,1,ielem,inode,amatrix,bmatrix,fterm_coef)

        end do

! Intergration on the body surface    

        do  ielem=1+nelemf, nelem

            call comp_link(ielem,inode,ii)!
            if (ii .eq. 0)   then 
                call norm_ele0(ielem,xp,yp,zp,amatrix,bmatrix)
            else if (ii .ne. 0)   then 
                call sing_ele0(inode,ielem,nodqua(inode),xp,yp,zp,&
                 &                   amatrix,bmatrix)
            end if

            call common_block(1,1,ielem,inode,amatrix,bmatrix,fterm_coef)
        end do
1000     continue
!
! =============================================

        if( nsys .eq. 2) then
            do inode=1, nnf
            if(nodqua(inode) .eq. 2) then
                amata(inode,inode,2)=1.0e20         
            endif
            enddo
        else if( nsys .eq. 4) then
            do inode=1, nnf
            if(nodqua(inode) .eq. 2) then
                amata(inode,inode,2)=1.0e20
                amata(inode,inode,4)=1.0e20            
            else if(nodqua(inode) .eq. 4) then
                amata(inode,inode,3)=1.0e20
                amata(inode,inode,2)=1.0e20
            else if(nodqua(inode) .eq. 5) then
                amata(inode,inode,2)=1.0e20
                amata(inode,inode,3)=1.0e20            
                amata(inode,inode,4)=1.0e20            
            endif
            enddo
        endif
        !-------output amata,bmata to txt file
        do i = 1,nnode
            do j = 1,nnode
                write(400,*) amata(i,j,1:nsys)
        end do;end do
        do i = 1,nnode
            write(401,*) bmata(i,1:nsys)
        end do
! =============================================
!
            DO IP=1, NSYS
        WRITE(101, *) '  IP=',IP
          WRITE(101, *) '    INODE=',nNODE,'      AMATA' 
        DO INODE=1,  NNODE
         WRITE(101, *) '  INODE=',INODE
         DO IND= 1,  NNODE 
         WRITE(101, 620) IND,XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),&
     &                    AMATA(INODE,IND,IP)
         ENDDO      
        ENDDO
       ENDDO
!
!
        WRITE(102, *) '  =========== Before RLUDCMP =============='
         DO IP=1, NSYS
        WRITE(102, *) '  IP=',IP
          WRITE(102, *) '    INODE          BMATA' 
         DO IND= 1,  NNODE 
         WRITE(102, 620) IND,XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),&
     &                    BMATA(IND,IP) 
         ENDDO      
       ENDDO
!       
!
!
        WRITE(102, *) 
        WRITE(102, *)
        WRITE(102, *) '  =========== After RLUDCMP =============='
           DO IP=1, NSYS
           WRITE(6, *) '  IP=',IP,'    Before RLUDCMP'
             CALL RLUDCMP(IP,AMATA,NNODE,NNODE,NSYS,INDX,DSIGN)  
           ENDDO
!
         DO IS=1, NSYS   
           CALL RLUBKSB(IS,AMATA,NNODE,NNODE,1,NSYS,1,INDX,BMATA)
         ENDDO
!         
!
         DO IP=1, NSYS
        WRITE(102, *) '  IP=',IP
          WRITE(102, *) '    INODE          BMATA' 
         DO IND= 1,  NNODE 
         WRITE(102, 620) IND,XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),&
     &                    BMATA(IND,IP) 
         ENDDO      
       ENDDO
!
! =======================================================                 
! ** output the results
!      
        DO 1500 IP=1, NSYS 
        DO 1360 IND=1, NNODE 
        BMAT(IP)=(0.0D0, 0.0D0)
        DO 1350 IS=1, NSYS
1350     BMAT(IP)=BMAT(IP)+BMATA(IND,IS)*RSN(IP,IS)
        BMAT(IP)=BMAT(IP)/NSYS
        UNKN(IND,IP)=BMAT(IP)
1360     CONTINUE
!
1500     CONTINUE


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!用于输出Ax=B中的x
!       IF(TIMERK .GT. 1.0) THEN
!  
         WRITE(9, *)
           WRITE(9, *)  '  Direvative of potential on the upper surface'
         WRITE(9, *)  '    INODE     XP      YP     ZP     Unkn    DPDZ'
           
        DO  INODE=1, NNF 
        XP=XYZ(1,INODE)
        YP=XYZ(2,INODE)
        ZP=XYZ(3,INODE) 
        CALL DINP(XP,YP,ZP,DPOX,DPOY,DPOZ)       

        WRITE(9, 620)  INODE,XP,YP,ZP,unkn(INODE,1),DPOZ
        ENDDO
        
         WRITE(9, *)
           WRITE(9, *)  '  Potential on the body surface'
         WRITE(9, *)  '    INODE     XP      YP     ZP     Unkn    POXY'

        DO  INODE=1+NNF, NNODE
        XP=XYZ(1,INODE)
        YP=XYZ(2,INODE)
        ZP=XYZ(3,INODE) 
        WRITE(9, 620)  INODE,XP,YP,ZP,unkn(INODE,1),POXY(XP,YP,ZP)
        ENDDO

        !do i = 1,NELEMF
          !write (200,999) i,ncon(i,1:8)
!C           x = 0.5*(XYZ(1,NCON(2))+XYZ(1,NCON(6)))
!C           Y = 0.5*(XYZ(2,NCON(2))+XYZ(2,NCON(6)))
        !end DO
!
         !do i = NELEMF+1,NELEM
          !write (201,999) i,ncon(i,1:8)
!C           x = 0.5*(XYZ(1,NCON(2))+XYZ(1,NCON(6)))
!C           Y = 0.5*(XYZ(2,NCON(2))+XYZ(2,NCON(6)))
        !end DO


!c 610     FORMAT(10x,'NNODE=',I6,/,2X,'IND',10X,'AMATA(IND,IND,1)')
!    
 999     format(9i4)
 610   FORMAT(10x,'NNODE     BMATA(IND,1)     T=',F14.6)
 620   FORMAT(1X,I4,6(1X,F13.6))      
 630   FORMAT(10x,'NNODE     UNKN(IND,1)     T=',F14.6)     
 640   FORMAT(1X,I4,2(2X,F13.6,2X,F13.6))                     
                    
                      
      RETURN
      END


