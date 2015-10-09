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
        integer :: inode,ielem
        integer,intent(out) :: ii
        integer :: i

        ii = 0 
        do i=1, nodnoe(inode)
        if(ielem .eq. nodele(inode,i)) then
          ii=ii+1
        endif
        enddo
    end subroutine 

    SUBROUTINE TASSB0
          USE MVAR_MOD
          USE PVAR_MOD
          USE MFUNC_mod
        USE SEBSM_MOD
        !use elem_intgl
        
        IMPLICIT   NONE  
          INTEGER  INODE,IELEM,J,JNODE,IND,INDD,IP
          integer ::              I,II,IS,JNCON,KNCON,L
          REAL*8  XP,YP,ZP,XSB,YSB,ZSB,R
          REAL*8  DX,DY,DZ,NX,NY,NZ
          REAL*8 RSN(4,4),EX(4),EY(4)
          
        REAL*8  BMATRIX(4,8),AMATRIX(4,8),BMAT(4)

          REAL*8  S_ANGLE
!         REAL*8  POXY 
          REAL*8  DSIGN
          real(8) :: fterm_coef(0:3,4)
        REAL*8  CELE31(4),CELE32(4),CELE33(4),AL1(4)
          REAL*8  DPOX,DPOY,DPOZ,DPDN,PHI


        DATA RSN /1.,  1.,  1.,  1., &
     &            1., -1.,  1., -1.,&
     &            1.,  1., -1., -1., &
     &            1., -1., -1.,  1./ 
!
          DATA EX /  1.0d0,  1.0d0, -1.0d0, -1.0d0/                       
        DATA EY /  1.0d0, -1.0d0, -1.0d0,  1.0d0/
!
!  ----------------------------------------------------
!
          WRITE(10, *)   ' IN TASSB0 '
          DSDT(:)=0.0
!         
!  ----------------------------------------------------
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

        amata(:,:,:)=(0.0d0,0.0d0)
        bmata(:,:)=(0.0d0,0.0d0)

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

        do  500   inode=1,  nnf   ! source point is on the free surface

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
!
!  ---------------------------
!  Integration on the free surface

            do  300   ielem=1,  nelemf

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

300     CONTINUE  

!  Integration on the body surface

        DO  400   IELEM=NELEMF+1,  NELEM
        WRITE(101,*)
        WRITE(101,*)  ' IELEM=',IELEM
        WRITE(101,*)  ' XP,YP,ZP=',XP,YP,ZP
         
        call comp_link(ielem,inode,ii)
!
        if (ii .eq. 0)   then 
            call norm_ele1(ielem,xp,yp,zp,amatrix,bmatrix)
        else if (ii .ne. 0)   then 
            call sing_ele1(inode,ielem,nodqua(inode),xp,yp,zp,&
             &                   amatrix,bmatrix)
        end if                
!       write(101,*) ' BMATRIX=',BMATRIX(1,:)
!       write(101,*) ' AMATRIX=',AMATRIX(1,:)
!!
!        do  320   j=1,  ncn(ielem) 
!          jncon=ncon(ielem,j)
!          kncon=ncond(ielem,j)
!         do  320   ip=1, nsys 
!             xsb=ex(ip)*xyz(1,jncon)
!             ysb=ey(ip)*xyz(2,jncon)
!             zsb=       xyz(3,jncon)
!
!             nx=ex(ip)*dxyz(1,kncon)
!             ny=ey(ip)*dxyz(2,kncon)
!             nz=       dxyz(3,kncon)
!
!             call dinp(xsb,ysb,zsb,dpox,dpoy,dpoz)       
!             dpdn=dpox*nx+dpoy*ny+dpoz*nz 
!
!             do  is=1, nsys    
!             if(jncon .gt. nnf)  then
!                 amata(inode,jncon,ip)=amata(inode,jncon,ip)-&
!         &                          rsn(is,ip)*amatrix(is,j)
!             else
!                 phi=poxy(xsb,ysb,zsb)
!         bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)*amatrix(is,j)*phi   !  * ******
!             endif
!
!         bmata(inode,ip)=bmata(inode,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn  !  * ******
!             enddo
!
!
!        do i = 0,3
!        call dinp0(i,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)       
!        dpdn=dpox*nx+dpoy*ny+dpoz*nz
!        do    is=1, nsys             
!        fterm_coef(i,ip)=fterm_coef(i,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn        
!        fterm_coef(i,ip)=fterm_coef(i,ip)+rsn(is,ip)*amatrix(is,j)*phi
!        enddo
!        end do 


        call common_block(1,0,ielem,inode,amatrix,bmatrix,fterm_coef)

400     CONTINUE

        FrA3(INODE)=fterm_coef(0,1)!AL1(1)
        FrC31(INODE)=fterm_coef(1,1)!CELE31(1)
        FrC32(INODE)=fterm_coef(2,1)!CELE32(1)
        FrC33(INODE)=fterm_coef(3,1)!CELE33(1)
   
         WRITE(10,620) INODE,ANGLE(INODE),FrA3(INODE),&
     &                 FrC31(INODE),FrC32(INODE),FrC33(INODE)

         phi=poxy(xp,yp,zp)
         call dinp(xp,yp,zp,dpox,dpoy,dpoz)       
         bmata(inode,1)=bmata(inode,1)-fra3(inode)*phi-&
     &                 frc31(inode)*dpox-frc32(inode)*dpoy
!        
500     CONTINUE
!
! =======================================================================
!    Source point is on the body surface
!
        DO  1000   INODE=NNF+1, NNODE   
!

        XP=XYZ(1,INODE)
        YP=XYZ(2,INODE)
        ZP=XYZ(3,INODE) 
!               
        CALL SOLIDANGLE(INODE,NNODE,NELEM,NCN,NCON,NODQUA,&
     &                    H,XYZ,DXYZE,S_ANGLE) 

            S_ANGLE=1.0d0-S_ANGLE

            WRITE(9,102)  INODE, XP, YP, ZP, S_ANGLE
            WRITE(*,102)  INODE, XP, YP, ZP, S_ANGLE

            ANGLE(INODE)=S_ANGLE

        DO    IP=1,  NSYS 
            AMATA(INODE,INODE,IP)= ANGLE(INODE)
        ENDDO
! ------------------------------
! Intergration on the free surface
! 
        DO  800   IELEM=1,  NELEMF
        II=0   
!C
!C Using TNORWP0 to integrate 
!C
         CALL NORM_ELE0(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)
        call common_block(0,1,ielem,inode,amatrix,bmatrix,fterm_coef)
!  ------
!
!        DO  710   J=1,  NCN(IELEM) 
!         JNCON=NCON(IELEM,J)
!         KNCON=NCOND(IELEM,J)
!        DO  710   IP=1, NSYS  
!          XSB=EX(IP)*XYZ(1,JNCON)
!            YSB=EY(IP)*XYZ(2,JNCON)
!            ZSB=       XYZ(3,JNCON)
!   
!          NX=EX(IP)*DXYZ(1,KNCON)
!            NY=EY(IP)*DXYZ(2,KNCON)
!            NZ=       DXYZ(3,KNCON)
!          CALL DINP(XSB,YSB,ZSB,DPOX,DPOY,DPOZ)       
!          DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ 
!         
!         DO  710   IS=1, NSYS    
!           AMATA(INODE,JNCON,IP)=AMATA(INODE,JNCON,IP)+&
!     &                           RSN(IS,IP)*BMATRIX(IS,J)
!             BMATA(INODE,IP)=BMATA(INODE,IP)+RSN(IS,IP)*AMATRIX(IS,J)*&
!     &               POXY(XSB,YSB,ZSB)           
!710     CONTINUE  
!

800     CONTINUE
          
!  --------------------------
! Intergration on the body surface    
! 
        DO  900   IELEM=1+NELEMF, NELEM
        II=0   
!C
!C Using TSING if the source point is in the element or its mirror
!C     elements about any symmetrical axis, otherwise using TINBOD
!C
        call comp_link(ielem,inode,ii)!
        IF (II .EQ. 0)   THEN 
         CALL NORM_ELE0(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)
        ELSE IF (II .NE. 0)   THEN 
         CALL SING_ELE0(INODE,IELEM,NODQUA(INODE),XP,YP,ZP,&
     &                   AMATRIX,BMATRIX)
        END IF                
!

         DO 820 J=1,  NCN(IELEM) 
          JNCON=NCON(IELEM,J)         
          KNCON=NCOND(IELEM,J)
         DO 820 IP=1, NSYS                   
           XSB=EX(IP)*XYZ(1,JNCON)
             YSB=EY(IP)*XYZ(2,JNCON)
             ZSB=       XYZ(3,JNCON)
           NX=EX(IP)*DXYZ(1,KNCON)
             NY=EY(IP)*DXYZ(2,KNCON)
             NZ=       DXYZ(3,KNCON)
          CALL DINP(XSB,YSB,ZSB,DPOX,DPOY,DPOZ)       
          DPDN=DPOX*NX+DPOY*NY+DPOZ*NZ 
          
           DO  IS=1, NSYS    
             IF(JNCON .GT. NNF)  THEN
           AMATA(INODE,JNCON,IP)=AMATA(INODE,JNCON,IP)-&
     &                          RSN(IS,IP)*AMATRIX(IS,J)
             ELSE       
              PHI=POXY(XSB,YSB,ZSB)                 
           BMATA(INODE,IP)=BMATA(INODE,IP)+RSN(IS,IP)*AMATRIX(IS,J)*PHI 
             ENDIF
            BMATA(INODE,IP)=BMATA(INODE,IP)-&
     &                      RSN(IS,IP)*BMATRIX(IS,J)*DPDN       
           ENDDO
820     CONTINUE
!
900     CONTINUE
!
1000     CONTINUE
!
! =============================================

           IF( NSYS .EQ. 2) THEN
!
                   DO INODE=1, NNF
                   IF(NODQUA(INODE) .EQ. 2) THEN
                           AMATA(INODE,INODE,2)=1.0E20         
                   ENDIF
                   ENDDO
!
           ELSE IF( NSYS .EQ. 4) THEN
!
                   DO INODE=1, NNF
                   IF(NODQUA(INODE) .EQ. 2) THEN
                           AMATA(INODE,INODE,2)=1.0E20
                           AMATA(INODE,INODE,4)=1.0E20            
                   ELSE IF(NODQUA(INODE) .EQ. 4) THEN
                           AMATA(INODE,INODE,3)=1.0E20
                           AMATA(INODE,INODE,2)=1.0E20
                   ELSE IF(NODQUA(INODE) .EQ. 5) THEN
                           AMATA(INODE,INODE,2)=1.0E20
                           AMATA(INODE,INODE,3)=1.0E20            
                           AMATA(INODE,INODE,4)=1.0E20            
                   ENDIF
                   ENDDO
           ENDIF
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


