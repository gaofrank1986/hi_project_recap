! 
! ********************************************** 
!        
!     4th order Runge_Kutta method 
! 
! ********************************************** 
! 
         SUBROUTINE Time_intg_RK4 
         USE MVAR_MOD 
         USE PVAR_MOD 
         use time_mod
         use wave_func,only:eti2
! 
         IMPLICIT NONE 
!         
         Integer INODE,JNODE,N0,NELE,NLOC,MOD,IP,IPLOT 
         Integer K 
         real(8)  VECT(3) 
         real(8)  AMPF(6),RAL 
 
! 
!  ============================================ 
!  RK1 
 
         RAL=180.0d0/PI 
 
         TimeRK=TIME 
!  Computing velocity on the body surface, and normal velocity of free surface by BEM 
!       
         CALL Runge_Kutta(1) 
 
        WRITE(*,*) 'RK1 COMPLETED' 
! 
!  ============================================ 
!  RK2 
! 
         TimeRK=TIME+Tstep/2.0d0 
!  
!  Compute wave force, and body response 
! 
 
         call runge_kutta(2) 
 
        WRITE(*,*) 'RK2 COMPLETED'  
! 
!  ============================================= 
!  R-K3 
 
          TimeRK=TIME+Tstep/2.0d0 
  
! 
          CALL Runge_Kutta(3) 
 
         WRITE(*,*) 'RK3 COMPLETED'  
! 
!  ============================================= 
!  R-K4 
! 
          TimeRK=TIME+Tstep 
!  
         CALL Runge_Kutta(4) 
 
        WRITE(*,*) 'RK4 COMPLETED' 
! 
! 
!C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
! 
        ! bkn is boudnary value
       do ip=1, nsys 
         do inode=1, nnf 
         !wave height
            et(inode,ip)=et_o(inode,ip)+tstep/6.0* &
     &                       (dh(1,inode,ip)+2.0*dh(2,inode,ip)+ &
     &               2.0*dh(3,inode,ip)+dh(4,inode,ip) ) 
        !wave potential on surface 
            bkn(inode,ip)=bkn_o(inode,ip)+tstep/6.0d0* &
     &                    (dp(1,inode,ip)+2.0d0*dp(2,inode,ip)+ &
     &                2.0d0*dp(3,inode,ip)+dp(4,inode,ip) ) 
 
       enddo 
       enddo 
       write(*,5001) dh(1:4,1,1)
       write(*,5001) dp(1:4,1,1)
       print *,"et=",et(1,1)
       print *,"bkn=",bkn(1,1)
       5001 format(4f14.8)

        if (itime.eq.endtime) then
        write (3000,*) et
        write (3001,*) bkn
        end if
 
         !IF (NPLOUT.EQ.1) THEN 
           !IPLOT=MOD(ITIME, 1) 
           !IF (IPLOT .EQ. 0  )  CALL PLOTOUT8 
         !ELSE IF (NPLOUT.EQ.2) THEN 
           !IPLOT=MOD(ITIME, 2) 
           !IF (IPLOT .EQ. 0  )  CALL PLOTOUT8 
         !ELSE IF (NPLOUT.EQ.3) THEN 
           !IPLOT=MOD(ITIME, 4) 
           !IF (IPLOT .EQ. 0  )  CALL PLOTOUT8 
         !ELSE IF (NPLOUT.EQ.4) THEN 
           !IPLOT=MOD(ITIME, 8) 
           !IF (IPLOT .EQ. 0  )  CALL PLOTOUT8 
         !ELSE IF (NPLOUT.EQ.5) THEN 
           !IPLOT=MOD(ITIME, 16) 
           !IF (IPLOT .EQ. 0  )  CALL PLOTOUT8 
         !ENDIF 
 
          FORCE=FORCEW 
          WRITE(21, 1010)  TIME, FORCE(1), FORCE(2), FORCE(3),& 
                              &        FORCE(4), FORCE(5), FORCE(6) 
 
!C 
! SG FIXME: try to make body fixed 
!C Displacement and d_phi/dt of the body 
!C 
         !DO 120 K=1, 6 
           !DSDT(K)=DSDT_O(K)+(Dposi(1,K)+2.0D0*Dposi(2,K)+ 
        !1                       2.0D0*Dposi(3,K)+Dposi(4,K))/6.0D0 
!120      CONTINUE 
 
         !DO 140 K=1, 6 
         !DISP(K)=DISP_O(K)+ 
        !1                Tstep*DSDT_O(K)+Tstep*(Dposi(1,K)+Dposi(2,K) 
     !2               +Dposi(3,K))/6.0D0 
!140      CONTINUE 
 
! 
 
!         IF(TIME .GT. 2.0d0*TPER) THEN 
 
          !WRITE(22, 1010)  TIME, DISP(1), DISP(2), DISP(3), &
        !&             RAL*DISP(4), RAL*DISP(5), RAL*DISP(6) 
!! 
          !WRITE(23, 1010)  TIME, DSDT(1),     DSDT(2), &
        !&               DSDT(3),     RAL*DSDT(4), &
        !&               RAL*DSDT(5), RAL*DSDT(6) 
!!         ENDIF 
         DO 150 K=1, 6 
        IF( DISP(K).GT. 1000.0d0) DISP(K)= 1000.0d0 
        IF( DISP(K).LT.-1000.0d0) DISP(K)=-1000.0d0 
150      CONTINUE 
 
 
! 
1010    FORMAT(F7.3,1x,7E12.4)  
! 
!C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* 
! 
10    FORMAT(1X,I6,2X,3(E16.9,3X)) 
21    FORMAT(1X,5(F10.5,3X)) 
22    FORMAT(1X,9(I7,2X)) 
        RETURN 
        END 
 
 
 
! 
!  ========================================== 
! 
! 
!  ========================================== 
! 
         SUBROUTINE Runge_Kutta(N) 
         USE MVAR_MOD 
         USE PVar_mod 
         use time_mod
         use mfunc_mod,only:rlubksb
 
! 
         IMPLICIT real(8) (A-H,O-Z) 
!  
         INTEGER INODE,N,J 
         real(8) BMAT(4),RSN(4,4) 
         real(8) COMP(6),AMPF(6)!,TRXYZ(3,3) 
 
 
       DATA RSN /1.,  1.,  1.,  1.,  &
     &           1., -1.,  1., -1., &
     &           1.,  1., -1., -1.,  &
     &           1., -1., -1.,  1./ 
! 
! 
! RAMPF: ramp function for incient potential 
! RAMPV: ramp function for damping  
! 
 
        WRITE(11,*)  
        WRITE(11,*) ' Inside  Runge_Kutta        N=',N 
 
        IF(TimeRK .LT. 2.0d0*TPER) THEN 
                RAMPF=0.5d0*(1.0d0-COS(PI*TimeRK/2.0d0/TPER)) 
        ELSE 
                RAMPF=1.0D0 
        END IF 
        !rampf=1.0d0
 
        IF(TimeRk .LT. 6.0d0*TPER) THEN 
                RAMPV=0.5d0*(3.0D0+Cos(PI*TimeRK/6.0d0/TPER)) 
        ELSE 
                RAMPV=1.0D0 
        END IF 
! 
!  ================================================ at RK1 step 
! 
          IF(N  .EQ. 1)   THEN 
! 
!  Computing wave height, potential at water surface  
! 
!         
          DISP(:)= DISP_O(:) 
          DSDT(:)= DSDT_O(:) 
!         prin 
          CALL TASSBT 
          !solve for unkn

          DO IP=1,    NSYS 
          DO INODE=1, NNF 
               DH(1,INODE,IP)=  UNKN(INODE,IP) -DAMPF(INODE)*ET(INODE,IP) 
               !d\eti / dt
               DP(1,INODE,IP)= -G*ET(INODE,IP) -DAMPF(INODE)*BKN(INODE,IP) 
               !{d \phi}/dt
          ENDDO   
          ENDDO 
          !print *,"after dh,dp"
!        print 
!  Compute wave force, and body response 
!  dpdt is d_phi/dt potential deritive over time??  
          !dsdt is velocity
           DPDT(1:NNF,:)=DP(1,1:NNF,:) 
           DPDT(NNF+1:NNODE,:)=(UNKN(NNF+1:NNODE,:) &
     &                      -UNKN_O(NNF+1:NNODE,:))/Tstep 

          CALL TFORCE                      
           !processing in time
            UNKN_O(:,:)=UNKN(:,:) 
 
! 
         ! CALL TRMAX(DISP,TRXYZ) 
          !CALL ACCEL(1,TRXYZ,COMP) 
         ! DO K=1, 6 
         !dposi is acceleration?? SG
          !Dposi(1,K)=COMP(K) ! can set to 0 if no force
          !END DO 
          Dposi(1,1:6)=0.0d0
!        Print *,' Rk1  DPosi(1)=',DPosi(1,1) 
 
! 
!  ================================================ at RK2 step 
! 
! 
          ELSE IF(N .EQ. 2)  THEN 
 
                  ET(:,:)      = ET_O(:,:)     +DH(1,:,:)*Tstep/2.0d0 
                  BKN(1:NNF,:) = BKN_O(1:NNF,:)+DP(1,:,:)*Tstep/2.0d0 
                  !boundary potential changed on surface

                  DISP(:)=DISP_O(:) +Tstep*DSDT_O(:)/2.0d0  !gives displacement
                  DSDT(:)= DSDT_O(:)+Dposi(1,:)/2.0d0!  gives velocity 
                  ! tassbt require dsdt and bkn updated before running
                  CALL TASSBT 
! 
                  DO IP=1,    NSYS 
                  DO INODE=1, NNF 
                  DH(2,INODE,IP)=  UNKN(INODE,IP) -DAMPF(INODE)*ET(INODE,IP) 
                  DP(2,INODE,IP)= -G*ET(INODE,IP) -DAMPF(INODE)*BKN(INODE,IP) 
                  ENDDO   
                  ENDDO   
 
 
!  
!  Compute wave force, and body response 
! 
                  DPDT(1:NNF,:)=DP(2,1:NNF,:) 
                  DPDT(NNF+1:NNODE,:)=(UNKN(NNF+1:NNODE,:) &
     &                            -UNKN_O(NNF+1:NNODE,:))/(0.5*Tstep) 


           
 
          CALL TFORCE 
! 
         ! CALL TRMAX(DISP,TRXYZ) 
!         CALL MESHT(TRXYZ) 
          !CALL ACCEL(2,TRXYZ,COMP) 
          !DO K=1, 6 
           ! Dposi(2,K)=COMP(K) 
          !END DO 
        Dposi(2,1:6) = 0.0d0 
!        Print *,' Rk2  DPosi(1)=',DPosi(2,1) 
! 
!  ================================================ at RK3 step 
! 
          ELSE IF(N .EQ. 3)  THEN 
 
           ET(:,:)      = ET_O(:,:)     +DH(2,:,:)*Tstep/2.0d0 
           BKN(1:NNF,:) = BKN_O(1:NNF,:)+DP(2,:,:)*Tstep/2.0d0 
           
            DISP(:)=DISP_O(:)   + &
        &                Tstep*DSDT_O(:)/2.0d0+Tstep*Dposi(1,:)/4.0d0 
            DSDT(:)= DSDT_O(:)+Dposi(2,:)/2.0d0 
 
          CALL TASSBT 
! 
          DO IP=1,    NSYS 
          DO INODE=1, NNF 
            DH(3,INODE,IP)=  UNKN(INODE,IP) -DAMPF(INODE)*ET(INODE,IP) 
          DP(3,INODE,IP)= -G*ET(INODE,IP) -DAMPF(INODE)*BKN(INODE,IP) 
          ENDDO   
          ENDDO   
! 
!  
!  Compute wave force, and body response 
! 
              DPDT(1:NNF,:)=DP(3,1:NNF,:) 
              DPDT(NNF+1:NNODE,:)=(UNKN(NNF+1:NNODE,:) &
     &                           -UNKN_O(NNF+1:NNODE,:))/(0.5*Tstep) 
 
 
!         write(10,*) 'DTDP   -RK3' 
!         DO J=1, NNODE 
!          write(10,*) J,DPDT(J,:) 
!         ENDDO 
           
 
          CALL TFORCE      
          !CALL TRMAX(DISP,TRXYZ) 
!         CALL MESHT(TRXYZ) 
          !CALL ACCEL(3,TRXYZ,COMP) 
          !DO K=1, 6 
          !  Dposi(3,K)=COMP(K) 
          !END DO 

        Dposi(3,1:6) = 0.0d0 
!        Print *,' Rk3  DPosi(1)=',DPosi(3,1) 
! 
!  ================================================ at RK4 step 
! 
          ELSE IF(N .EQ. 4)  THEN 
 
           ET(:,:)      = ET_O(:,:)     +DH(3,:,:)*Tstep 
           BKN(1:NNF,:) = BKN_O(1:NNF,:)+DP(3,:,:)*Tstep 
 
            DISP(:)=DISP_O(:)   + &
        &                        Tstep*DSDT_O(:)+Tstep*Dposi(2,:)/2.0d0 
            DSDT(:)= DSDT_O(:)+Dposi(3,:) 
! 
!          HEIGHT(4,:,:)=HEIGHT(4,:,:)+DH(3,:,:)*Tstep 
!          PFREEN(4,:,:)=PFREEN(4,:,:)+DP(3,:,:)*Tstep 
!          BKN(:,:) = BKN_O(:,:) 
 
          CALL TASSBT 
! 
          DO IP=1,    NSYS 
          DO INODE=1, NNF 
            DH(4,INODE,IP)=  UNKN(INODE,IP) -DAMPF(INODE)*ET(INODE,IP) 
          DP(4,INODE,IP)= -G*ET(INODE,IP) -DAMPF(INODE)*BKN(INODE,IP) 
          ENDDO   
          ENDDO   
! 
!  
!  Compute wave force, and body response 
! 
              DPDT(1:NNF,:)=DP(4,1:NNF,:) 
              DPDT(NNF+1:NNODE,:)=(UNKN(NNF+1:NNODE,:) &
     &                           -UNKN_O(NNF+1:NNODE,:))/Tstep 
           
 
          CALL TFORCE 
 
                           
          !CALL TRMAX(DISP,TRXYZ) 
!         CALL MESHT(TRXYZ) 
 
          !CALL ACCEL(4,TRXYZ,COMP) 
          !DO K=1, 6 
           ! Dposi(4,K)=COMP(K) 
          !END DO 
!        
        Dposi(4,1:6) = 0.0D0 
!        Print *,' Rk4  DPosi(1)=',DPosi(4,1) 
 
          ENDIF 
! 
!    pause 
!         
        
         END 
! 
! ==================================================== 
! 
! 
! 
! 
         SUBROUTINE PLOTOUT8 
          USE MVAR_MOD 
          use time_mod
         use wave_func,only:eti2
        IMPLICIT   NONE   
 
         CHARACTER*16 NAME 
         CHARACTER*6  FIRST 
 
         INTEGER  I,INODE,NE 
       real(8) XP,YP,ETin 
        !real(8),EXTERNAL::  ETI2 
 
         CALL OPENFILE(ITIME,FIRST) 
! 
!         Print *,' Itime=',itime,' First=',first 
! 
         NAME='OUTtime\WAVE'//FIRST//'.DAT' 
       OPEN(102,FILE=NAME,STATUS='UNKNOWN') 
! 
!        OPEN(9,  FILE='OUTPUT\OUTPUT1.txt',    STATUS='UNKNOWN') 
!               
         WRITE(102,*) 'TITLE = "3D Mesh Grid Data for Element Boundary"' 
         WRITE(102,*) 'VARIABLES = "X", "Y", "Z"' 
         WRITE(102,*) 'ZONE N=',NNF,',','E=',NELEMF,',',& 
     &  'F=FEPOINT,ET=QUADRILATERAL' 
! 
       DO INODE=1, NNF 
        XP=XYZ(1,INODE) 
        YP=XYZ(2,INODE) 
!        ETin=ETI1(Amp,BETA,WK,W1,TimeRk,RampF,XP,YP) 
        ETin=ETI2(H,G,Ampn,Phi_w,BETA,WKN,FREQ,Time,RampF, &
     &            XP,YP,NFreq,Nwave,IOrder) 
          WRITE(102,21)  XP,YP, ET(INODE,1),ETin 
!         WRITE(102,21)  XP,YP, ET(INODE,1)+ETin 
       ENDDO 
 
         DO NE=1,  NELEMF 
        WRITE(102,22) NCON(NE,1),NCON(NE,3),NCON(NE,5),NCON(NE,7) 
         ENDDO 
! 
         Close(102) 
! 
10     FORMAT(1X,I6,2X,3(E16.9,3X)) 
21     FORMAT(1X,5(F10.5,3X)) 
22     FORMAT(1X,9(I7,2X)) 
! 
         END 
 
! 
! ==================================================== 
! 
         SUBROUTINE OPENFILE(MTIME,FIRST) 
        !INTEGER MTIME 
         CHARACTER*6 FIRST 
          
         IF(MTIME.LT.10) THEN 
          MN0=MTIME 
          FIRST=CHAR(MN0+48) 
         ELSE IF(MTIME.LT.100) THEN 
          MN0=MOD(MTIME,10)      
          MN1=(MTIME-MN0)/10 
          FIRST=CHAR(MN1+48)//CHAR(MN0+48) 
         ELSE IF(MTIME.LT.1000) THEN 
          MN01=MOD(MTIME,100) 
          MN0=MOD(MN01,10)       
          MN1=(MN01-MN0)/10 
          MN2=(MTIME-10*MN1-MN0)/100 
          FIRST=CHAR(MN2+48)//CHAR(MN1+48)//CHAR(MN0+48) 
         ELSE IF(MTIME.LT.10000) THEN  
          MN012=MOD(MTIME,1000) 
          MN01=MOD(MN012,100) 
          MN0=MOD(MN01,10) 
          MN1=(MN01-MN0)/10 
          MN2=(MN012-10*MN1-MN0)/100 
          MN3=(MTIME-MN012)/1000 
          FIRST=CHAR(MN3+48)//CHAR(MN2+48)//CHAR(MN1+48)//CHAR(MN0+48) 
         ENDIF 
! 
! *-*-*--* 
! 
         RETURN 
         END 
 
