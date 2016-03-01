module wave2! 
! =================================================== 
!    
!  Incident wave profile  
!  The first and the second order incident wave profile 
!                                                
!       June 23, 2015           by   Bin Teng            
! =================================================== 
! 
! 
contains
        function eti2(h,g,ampn,phi,beta,wkn,freq,time,ramp,x,y,nn,nwave,iorder) 
  
 !DEC$ATTRIBUTES DLLEXPORT::ETI2 
 
        IMPLICIT  NONE 
! 
!            
        real*8  ETI2 
        Integer,INTENT(IN):: NN,Nwave,IOrder 
            real*8,INTENT(IN):: H,G,Ampn(NN),Phi(NN),BETA,WKN(NN),FREQ(NN),Time,Ramp,X,Y 
!         
        Integer I,J 
        real*8  RXY,WK_S,WK_D,WKX,FREQ_S,FREQ_D,PHY_S,PHY_D 
        real*8  H2_S,H2_D,GD_S,GD_D 
        real*8  DUM1,DUM2,DUM3,DUM4,DUM5,SH3 
! 
           ETI2=0.0d0 
          
           RXY=X*COS(BETA)+Y*SIN(BETA) 
!   
        IF (IORDER==1)  THEN     
          DO I=1, Nwave 
            WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
            ETI2=ETI2+AmpN(I)*Cos(WKX) 
              ENDDO           
!             
        ELSE IF (IORDER==2)  THEN 
              
             IF(H .LE. 0.0D0) THEN 
                DO I=1, Nwave 
                  WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                   ETI2=ETI2+0.5*AMPN(I)*AMPN(I)*WKN(I)*COS(2*WKX) 
                END DO 
                  
               IF (Nwave .ge. 2) then     
                 DO J=1, Nwave-1 
                   DO I=J+1, Nwave 
                     WK_S=WKN(I)+WKN(J) 
                     WK_D=WKN(I)-WKN(J) 
                     FREQ_S=FREQ(I)+FREQ(J) 
                    FREQ_D=FREQ(I)-FREQ(J) 
                     PHY_S=Phi(I)+Phi(J) 
                     PHY_D=Phi(I)-Phi(J) 
                     DUM1=WK_S*COS(WK_S*RXY-FREQ_S*TIME+PHY_S) 
                     DUM2=WK_D*COS(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                     ETI2=ETI2+0.5*AMPN(I)*AMPN(J)*(DUM1-DUM2) 
                     END DO 
                  END DO         
                END IF                    
                   
            ELSE IF(H .GT. 0.0D0) THEN                     
                  
               ETI2=0.0    
                DO I=1, Nwave 
                  SH3=( SINH(WKN(I)*H) )**3 
                 WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                  ETI2=ETI2+AMPN(I)*AMPN(I)*WKN(I)/4*COSH(WKN(I)*H)*(2.0+   & 
                          COSH(2.0*WKN(I)*H))*COS(2*WKX)/SH3 
                END DO 
 
               IF (Nwave .ge. 2) THEN 
                 DO J=1,  Nwave-1 
                  DO I=J+1,  Nwave 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
 
                    DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_S              & 
                         *(1-TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2+    & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2             
                    DUM5=G*WK_S*TANH(WK_S*H)-FREQ_S**2 
         
                    GD_S=-G*G*(DUM3+DUM4)/DUM5 
 
                    DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_D             & 
                           *(1+TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2-   & 
                           0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2  
                   DUM5=G*WK_D*TANH(WK_D*H)-FREQ_D**2           
         
                    GD_D=-G*G*(DUM3+DUM4)/DUM5 
 
                    H2_S=FREQ_S/G*GD_S-0.5*G*WKN(I)*WKN(J)/FREQ(I)/FREQ(J)    & 
                         *COSH(WK_D*H)/COSH(WKN(I)*H)/COSH(WKN(J)*H)+0.5*     & 
                         (WKN(I)*TANH(WKN(I)*H)+WKN(J)*TANH(WKN(J)*H)) 
 
                    H2_D=FREQ_D/G*GD_D-0.5*G*WKN(I)*WKN(J)/FREQ(I)/FREQ(J)    & 
                          *COSH(WK_S*H)/COSH(WKN(I)*H)/COSH(WKN(J)*H)+0.5*    & 
                         (WKN(I)*TANH(WKN(I)*H)+WKN(J)*TANH(WKN(J)*H)) 
 
                   DUM1=H2_S*COS(WK_S*RXY-FREQ_S*TIME+PHY_S) 
  
                   DUM2=H2_D*COS(WK_D*RXY-FREQ_D*TIME+PHY_D) 
 
                    ETI2=ETI2+AMPN(I)*AMPN(J)*(DUM1+DUM2) 
                  END DO 
                END DO 
               END IF  
              END IF 
              
         ENDIF 
              
         ETI2=RAMP*ETI2 
          
        RETURN 
        END FUNCTION ETI2 
! 
! ======================================================== 
! X and Y derivatives of first and second order incident wave profile  
! 
! ======================================================== 
! 
        SUBROUTINE DETW2DXY(H,G,Ampn,Phi,BETA,WKN,FREQ,Time,Ramp,X,Y,NN,Nwave,IOrder, & 
                             DETDX,DETDY) 
          
   !DEC$ATTRIBUTES DLLEXPORT::DETW2DXY 
       
        IMPLICIT  NONE 
! 
        INTEGER, INTENT(IN)::  NN,Nwave,IOrder 
         real*8, INTENT(IN):: H,G,Ampn(NN),Phi(NN),BETA,WKN(NN),FREQ(NN),Time,Ramp,X,Y 
         real*8, INTENT(OUT)::  DETDX,DETDY 
 
        INTEGER I,J 
        real*8  RXY,WKX,WK_S,WK_D,FREQ_S,FREQ_D,PHY_S,PHY_D 
        real*8  H2_S,H2_D,GD_S,GD_D 
        real*8  DUM1,DUM2,DUM3,DUM4,DUM5,SH3 
 
          RXY=X*COS(BETA)+Y*SIN(BETA) 
 
           DETDX=0.0D0 
          DETDY=0.0D0 
 
      IF(IORDER.EQ.1) THEN 
 
           DO I=1,  Nwave 
              WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I)            
              DETDX=DETDX-AMPN(I)*WKN(I)*COS(BETA)*SIN(WKX) 
              DETDY=DETDY-AMPN(I)*WKN(I)*SIN(BETA)*SIN(WKX) 
           END DO 
 
      ELSE IF(IORDER.EQ.2) THEN 
           
            IF(H .LE. 0.0D0) THEN 
                DO I=1, Nwave 
                 WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                  DETDX=DETDX-AMPN(I)*AMPN(I)*WKN(I)**2*2*COS(BETA)*SIN(2*WKX) 
                 DETDY=DETDY-AMPN(I)*AMPN(I)*WKN(I)**2*2*SIN(BETA)*SIN(2*WKX) 
                END DO 
                  
               IF (Nwave .ge.2) THEN    
                 DO J=1, Nwave-1 
                    DO I=J+1, Nwave 
                        WK_S=WKN(I)+WKN(J) 
                        WK_D=WKN(I)-WKN(J) 
                        FREQ_S=FREQ(I)+FREQ(J) 
                        FREQ_D=FREQ(I)-FREQ(J) 
                        PHY_S=Phi(I)+Phi(J) 
                        PHY_D=Phi(I)-Phi(J) 
                        DUM1=-WK_S**2*COS(BETA)*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S) 
                       DUM2=-WK_D**2*COS(BETA)*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                       DUM3=-WK_S**2*SIN(BETA)*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S) 
                       DUM4=-WK_D**2*SIN(BETA)*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                       DETDX=DETDX+0.5*AMPN(I)*AMPN(J)*(DUM1-DUM2) 
                      DETDY=DETDY+0.5*AMPN(I)*AMPN(J)*(DUM3-DUM4) 
                     END DO 
                  END DO                                  
                 END IF  
               
              ELSE      IF(H .GT. 0.0D0) THEN                      
                     
                 DO I=1, Nwave 
                    SH3=( SINH(WKN(I)*H) )**3 
                   WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                    DETDX=DETDX-AMPN(I)*AMPN(I)*WKN(I)**2*COS(BETA)/2*COSH(WKN(I)*H)*(2.0+   & 
                           COSH(2.0*WKN(I)*H))*SIN(2*WKX)/SH3 
                   DETDY=DETDY-AMPN(I)*AMPN(I)*WKN(I)**2*SIN(BETA)/2*COSH(WKN(I)*H)*(2.0+   & 
                           COSH(2.0*WKN(I)*H))*SIN(2*WKX)/SH3 
                END DO 
 
               IF (Nwave .ge. 2) THEN 
                 DO J=1,  Nwave-1 
                   DO I=J+1,  Nwave 
                      WK_S=WKN(I)+WKN(J) 
                      WK_D=WKN(I)-WKN(J) 
                      FREQ_S=FREQ(I)+FREQ(J) 
                     FREQ_D=FREQ(I)-FREQ(J) 
                      PHY_S=Phi(I)+Phi(J) 
                      PHY_D=Phi(I)-Phi(J) 
 
                    DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_S              & 
                         *(1-TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2+    & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2             
                    DUM5=G*WK_S*TANH(WK_S*H)-FREQ_S**2 
         
                    GD_S=-G*G*(DUM3+DUM4)/DUM5 
 
                    DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_D             & 
                           *(1+TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2-   & 
                           0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2  
                   DUM5=G*WK_D*TANH(WK_D*H)-FREQ_D**2           
         
                    GD_D=-G*G*(DUM3+DUM4)/DUM5 
 
                    H2_S=FREQ_S/G*GD_S-0.5*G*WKN(I)*WKN(J)/FREQ(I)/FREQ(J)    & 
                         *COSH(WK_D*H)/COSH(WKN(I)*H)/COSH(WKN(J)*H)+0.5*     & 
                         (WKN(I)*TANH(WKN(I)*H)+WKN(J)*TANH(WKN(J)*H)) 
 
                    H2_D=FREQ_D/G*GD_D-0.5*G*WKN(I)*WKN(J)/FREQ(I)/FREQ(J)    & 
                          *COSH(WK_S*H)/COSH(WKN(I)*H)/COSH(WKN(J)*H)+0.5*    & 
                         (WKN(I)*TANH(WKN(I)*H)+WKN(J)*TANH(WKN(J)*H)) 
 
                   DUM1=-H2_S*WK_S*COS(BETA)*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S) 
  
                   DUM2=-H2_D*WK_D*COS(BETA)*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
 
                   DUM3=-H2_S*WK_S*SIN(BETA)*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S) 
  
                   DUM4=-H2_D*WK_D*SIN(BETA)*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
 
                    DETDX=DETDX+AMPN(I)*AMPN(J)*(DUM1+DUM2) 
                   DETDY=DETDY+AMPN(I)*AMPN(J)*(DUM3+DUM4) 
                  END DO 
                END DO 
               END IF  
             
             END IF 
 
       ENDIF 
 
        DETDX=RAMP*DETDX 
        DETDY=RAMP*DETDY 
 
        RETURN 
        END SUBROUTINE DETW2DXY 
! 
! ====================================================================== 
! 
! The time derivatives of second order incident wave profile  
! 
! ====================================================================== 
! 
! 
        FUNCTION DETW2DT(H,G,Ampn,Phi,BETA,WKN,FREQ,Time,Ramp,X,Y,NN,Nwave,IOrder) 
   
   !DEC$ATTRIBUTES DLLEXPORT:: DETW2DT 
 
       IMPLICIT  NONE 
! 
       INTEGER, INTENT(IN)::  NN,Nwave,IOrder 
           real*8, INTENT(IN):: H,G,Ampn(NN),Phi(NN),BETA,WKN(NN),FREQ(NN),Time,Ramp,X,Y 
 
           INTEGER  I,J 
 
       real*8  DETW2DT 
       real*8  RXY,WKX 
       real*8  WK_S,WK_D,FREQ_S,FREQ_D,PHY_S,PHY_D 
           real*8  SH3,DUM1,DUM2,DUM3,DUM4,DUM5 
           real*8  H2_S,H2_D,GD_S,GD_D                
 
    
         RXY=X*COS(BETA)+Y*SIN(BETA) 
 
    
      IF(IORDER.EQ.1) THEN 
           
           DETW2DT=0.0D0 
           DO I=1, Nwave         
              WKX=WKN(I)*RXY-FREQ(I)*TIME+Phi(I)     
              DETW2DT=DETW2DT+AMPN(I)*FREQ(I)*SIN(WKX) 
           END DO 
           
      ELSE IF(IORDER.EQ.2) THEN 
   
            IF(H .LE. 0.0D0) THEN 
              DETW2DT=0.0D0 
              DO I=1, Nwave      
              WKX=WKN(I)*RXY-FREQ(I)*TIME+Phi(I)     
              DETW2DT=DETW2DT+AMPN(I)*AMPN(I)*WKN(I)*FREQ(I)*SIN(2*WKX) 
              END DO 
 
             IF (Nwave .GE. 2) THEN 
              DO J=1,Nwave-1 
                DO I=J+1,Nwave 
                  WK_S=WKN(I)+WKN(J) 
                  WK_D=WKN(I)-WKN(J) 
                  FREQ_S=FREQ(I)+FREQ(J) 
                FREQ_D=FREQ(I)-FREQ(J) 
                  PHY_S=Phi(I)+Phi(J) 
                  PHY_D=Phi(I)-Phi(J) 
                  DUM1=WK_S*FREQ_S*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S)  
                DUM2=WK_D*FREQ_D*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                DETW2DT=DETW2DT+0.5*AMPN(I)*AMPN(J)*(DUM1-DUM2) 
                END DO 
               END DO 
              END IF 
 
            ELSE 
                DETW2DT=0.0D0 
               DO I=1,Nwave 
                  SH3=( SINH(WKN(I)*H) )**3 
                  DETW2DT=DETW2DT+AMPn(I)*AMPn(I)*WKN(I)*FREQ(I)/2*    & 
                          COSH(WKN(I)*H)*(2.0+COSH(2.0*WKN(I)*H))*SIN(2*WKX)/SH3 
               END DO 
 
              IF (Nwave .GE. 2) THEN 
               DO J=1,Nwave-1 
                 DO I=J+1,Nwave 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
 
                    DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_S   & 
                         *(1-TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2+   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2             
                    DUM5=G*WK_S*TANH(WK_S*H)-FREQ_S**2 
         
                    GD_S=-G*G*(DUM3+DUM4)/DUM5 
 
                    DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_D*(1+TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2-   & 
                              0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2  
                   DUM5=G*WK_D*TANH(WK_D*H)-FREQ_D**2           
         
                   GD_D=-G*G*(DUM3+DUM4)/DUM5 
 
                    H2_S=FREQ_S/G*GD_S-0.5*G*WKN(I)*WKN(J)/FREQ(I)/FREQ(J)     & 
                         *COSH(WK_D*H)/COSH(WKN(I)*H)/COSH(WKN(J)*H)          & 
                         +0.5*(WKN(I)*TANH(WKN(I)*H)+WKN(J)*TANH(WKN(J)*H)) 
 
                    H2_D=FREQ_D/G*GD_D-0.5*G*WKN(I)*WKN(J)/FREQ(I)/FREQ(J)    & 
                         *COSH(WK_S*H)/COSH(WKN(I)*H)/COSH(WKN(J)*H)         & 
                         +0.5*(WKN(I)*TANH(WKN(I)*H)+WKN(J)*TANH(WKN(J)*H)) 
 
                   DUM1=H2_S*FREQ_S*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S) 
  
                   DUM2=H2_D*FREQ_D*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
 
                   DETW2DT=DETW2DT+AMPN(I)*AMPN(J)*(DUM1+DUM2) 
                END DO 
              END DO 
             END IF 
 
            END IF 
 
         END IF 
 
         DETW2DT=RAMP*DETW2DT 
 
        RETURN 
        END FUNCTION DETW2DT 
         
!  ============================================== 
! Incident potential  
!  accurate to the first order or the second order) 
! 
!                                                
!       June 23, 2015           by   Bin Teng            
!  ============================================== 
! 
 
        FUNCTION POXY2(H,G,Ampn,Phi,BETA,WKN,Freq,Time,Ramp,X,Y,Z,NN,Nwave,IOrder) 
!DEC$ATTRIBUTES DLLEXPORT::POXY2 
      
        IMPLICIT  NONE 
 
        INTEGER, INTENT(IN)::  NN,Nwave,IOrder 
            real*8,INTENT(IN):: H,G,Ampn(NN),Phi(NN),BETA,WKN(NN),FREQ(NN),Time,Ramp,X,Y,Z 
         
        INTEGER I,J 
         
        real*8  POXY2 
         real*8  RXY,WKX,CSHKH1 
! 
        real*8  WK_S,WK_D,FREQ_S,FREQ_D,PHY_S,PHY_D 
         real*8  GD_S,GD_D 
         real*8  SH4,DUM,DUM1,DUM2,DUM3,DUM4,DUM5 
 
          POXY2=0.0 
         RXY=X*COS(BETA)+Y*SIN(BETA) 
! 
         IF(IORDER.EQ.1) THEN 
 
             IF(H .LE. 0.0D0) THEN 
                  DO I=1, NWAVE 
                    DUM=AMPN(I)*G/FREQ(I)*EXP(WKN(I)*Z) 
                    WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                    POXY2=POXY2+DUM*SIN(WKX) 
                END DO 
              ELSE 
                  DO I=1,  NWAVE 
                    DUM=AMPN(I)*G/FREQ(I)*COSH(WKN(I)*(Z+H))/COSH(WKN(I)*H) 
                    WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                    POXY2=POXY2+DUM*SIN(WKX) 
                END DO 
            END IF 
               
        ELSE IF(IORDER.EQ.2) THEN 
 
              IF(H .LE. 0.0D0) THEN 
                  POXY2=0.0 
                    
                   IF (Nwave .ge. 2) THEN 
                    DO J=1,NWAVE-1 
                      DO I=J+1, NWAVE 
                         WK_D=WKN(I)-WKN(J) 
                         FREQ_D=FREQ(I)-FREQ(J) 
                         PHY_D=Phi(I)-Phi(J) 
                         POXY2=POXY2-FREQ(I)*AMPN(I)*AMPN(J)*EXP(WK_D*Z)*     & 
                                   SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                      END DO 
                    END DO 
                  END IF 
              ELSE 
                   POXY2=0.0 
                    DO I=1, NWAVE 
                       SH4=(SINH(WKN(I)*H))**4 
                        DUM=3.0*FREQ(I)*AMPN(I)*AMPN(I)/8.0*          & 
                                    COSH(2.0*WKN(I)*(Z+H))/SH4 
                       WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                        POXY2=POXY2+DUM*SIN(2.0*WKX) 
                    END DO 
                     
                  IF (Nwave .ge. 2) THEN 
                    DO J=1, NWAVE-1 
                      DO I=J+1, NWAVE 
                         WK_S=WKN(I)+WKN(J) 
                         WK_D=WKN(I)-WKN(J) 
                         FREQ_S=FREQ(I)+FREQ(J) 
                        FREQ_D=FREQ(I)-FREQ(J) 
                         PHY_S=Phi(I)+Phi(J) 
                         PHY_D=Phi(I)-Phi(J) 
 
                        DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_S           & 
                                *(1-TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                        DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2+   & 
                                0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2             
                        DUM5=G*WK_S*TANH(WK_S*H)-FREQ_S**2 
         
                        GD_S=-G*G*(DUM3+DUM4)/DUM5 
 
                        DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_D  & 
                                *(1+TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                        DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2-   & 
                                0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2  
                        DUM5=G*WK_D*TANH(WK_D*H)-FREQ_D**2           
         
                        GD_D=-G*G*(DUM3+DUM4)/DUM5 
                       DUM1=GD_S*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S)*COSH(WK_S*(Z+H))/COSH(WK_S*H) 
                       DUM2=GD_D*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D)*COSH(WK_D*(Z+H))/COSH(WK_D*H) 
                   
                         POXY2=POXY2+AMPN(I)*AMPN(J)*(DUM1+DUM2) 
                      END DO 
                    END DO 
                  END IF 
                 END IF                               
           
             ENDIF 
! 
         POXY2=RAMP*POXY2 
 
 
         RETURN 
         END FUNCTION POXY2 
! 
! =================================================== 
! Time Derivatives of incident wave potential 
!    
! 
        FUNCTION DPOT2(H,G,Ampn,Phi,BETA,WKN,Freq,Time,Ramp,X,Y,Z,NN,Nwave,IOrder) 
!DEC$ATTRIBUTES DLLEXPORT::DPOT2 
! 
        IMPLICIT  NONE 
! 
        real*8  DPOT2 
        Integer,INTENT(IN):: NN,Nwave,IOrder 
! 
            real*8,INTENT(IN):: H,G,Ampn(NN),Phi(NN),BETA,WKN(NN),FREQ(NN),Time,Ramp,X,Y,Z 
 
        INTEGER  I,J 
         
            real*8  WKX,DUM,CSHKH1 
 
        real*8  RXY,WK_S,WK_D,FREQ_S,FREQ_D,PHY_S,PHY_D 
        real*8  H2_S,H2_D,GD_S,GD_D 
        real*8  DUM1,DUM2,DUM3,DUM4,DUM5,SH3,SH4 
 
! 
         RXY=X*COS(BETA)+Y*SIN(BETA) 
 
      IF(IORDER.EQ.1) THEN 
 
             IF(H .LE. 0.0D0) THEN 
               DPOT2=0.0 
                 DO I=1, NWAVE 
                DUM=-AMPN(I)*G*EXP(WKN(I)*Z) 
                WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                DPOT2=DPOT2+DUM*COS(WKX) 
                END DO 
              ELSE 
               DPOT2=0.0 
                DO I=1,  NWAVE 
                DUM=-AMPN(I)*G*COSH(WKN(I)*(Z+H))/COSH(WKN(I)*H) 
                WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                DPOT2=DPOT2+DUM*COS(WKX) 
                END DO 
            END IF 
               
      ELSE IF(IORDER.EQ.2) THEN 
 
              IF(H .LE. 0.0D0) THEN 
                 DPOT2=0.0 
                   
                  IF (Nwave .ge. 2) THEN 
                    DO J=1,NWAVE-1 
                      DO I=J+1, NWAVE 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_D=Phi(I)-Phi(J) 
                    DPOT2=DPOT2+FREQ(I)*FREQ_D*AMPN(I)*AMPN(J)*EXP(WK_D*Z)*     & 
                                   COS(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                      END DO 
                    END DO 
                   END IF 
              ELSE 
               DPOT2=0.0 
                DO I=1, NWAVE 
               SH4=(SINH(WKN(I)*H))**4 
                  DUM=-3.0*FREQ(I)**2*AMPN(I)**2/4.0*          & 
                        COSH(2.0*WKN(I)*(Z+H))/SH4 
                 WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                  DPOT2=DPOT2+DUM*COS(2.0*WKX) 
                    END DO 
 
                  IF (Nwave .ge. 2) THEN 
                    DO J=1, NWAVE-1 
                      DO I=J+1, NWAVE 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
 
                    DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_S           & 
                         *(1-TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2+   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2             
                    DUM5=G*WK_S*TANH(WK_S*H)-FREQ_S**2 
         
                    GD_S=-G*G*(DUM3+DUM4)/DUM5 
 
                    DUM3=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_D  & 
                         *(1+TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM4=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2-   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2  
                   DUM5=G*WK_D*TANH(WK_D*H)-FREQ_D**2           
         
                    GD_D=-G*G*(DUM3+DUM4)/DUM5 
                    DUM1=-GD_S*FREQ_S*COS(WK_S*RXY-FREQ_S*TIME+PHY_S)*COSH(WK_S*(Z+H))/COSH(WK_S*H) 
                    DUM2=-GD_D*FREQ_D*COS(WK_D*RXY-FREQ_D*TIME+PHY_D)*COSH(WK_D*(Z+H))/COSH(WK_D*H) 
                   
                   DPOT2=DPOT2+AMPN(I)*AMPN(J)*(DUM1+DUM2) 
                    END DO 
                  END DO 
               END IF 
              END IF                                          
           
         ENDIF 
         DPOT2=RAMP*DPOT2 
 
        RETURN 
        END FUNCTION DPOT2 
 
 
! 
! 
!============================================================================ 
!  
!    Spacial Derivatives of the first and the second order components  
!    of incident wave potential 
!       IORDER=1: for the first order potential 
!       IORDER=2: for the second order potential  
! 
!============================================================================ 
 
        SUBROUTINE  DPOXYZ(H,G,Ampn,Phi,BETA,WKN,Freq,Time,Ramp,X,Y,Z,NN,Nwave,IOrder, & 
                     DPOX,DPOY,DPOZ) 
 
!DEC$ATTRIBUTES DLLEXPORT::DPOXYZ 
 
          IMPLICIT    NONE 
            
         INTEGER, INTENT(IN)::  NN,Nwave,IOrder 
          real*8,INTENT(IN)::   H,G,Ampn(NN),Phi(NN),BETA,WKN(NN),FREQ(NN),Time,Ramp,X,Y,Z 
         real*8,INTENT(OUT)::  DPOX,DPOY,DPOZ 
 
          INTEGER I,J 
 
           real*8  RXY,WKX,SH4 
          real*8  WK_S,WK_D,FREQ_S,FREQ_D,PHY_S,PHY_D 
           real*8  DUM,DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7 
           real*8  GD_S,GD_D 
 
! 
         RXY=X*COS(BETA)+Y*SIN(BETA) 
 
       IF(IORDER.EQ.1) THEN 
           
              IF(H .LE. 0.0D0) THEN 
                DPOX=0.D0 
                DPOY=0.D0 
                DPOZ=0.D0 
                DO I=1,Nwave 
                    DUM=AMPN(I)*G/FREQ(I) 
                   WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                   DPOX=DPOX+DUM*WKN(I)*COS(BETA)*EXP(WKN(I)*Z)*COS(WKX) 
                   DPOY=DPOY+DUM*WKN(I)*SIN(BETA)*EXP(WKN(I)*Z)*COS(WKX) 
                   DPOZ=DPOZ+DUM*WKN(I)*EXP(WKN(I)*Z)*SIN(WKX) 
                END DO 
              ELSE 
                DPOX=0.D0 
                DPOY=0.D0 
                DPOZ=0.D0 
                DO I=1,Nwave 
                    DUM=AMPN(I)*G/FREQ(I) 
                   WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I)             
                   DPOX=DPOX+DUM*WKN(I)*COS(BETA)*COSH(WKN(I)*(Z+H))/    & 
                               COSH(WKN(I)*H)*COS(WKX) 
                   DPOY=DPOY+DUM*WKN(I)*SIN(BETA)*COSH(WKN(I)*(Z+H))/    & 
                               COSH(WKN(I)*H)*COS(WKX) 
                   DPOZ=DPOZ+DUM*WKN(I)*SINH(WKN(I)*(Z+H))/     & 
                               COSH(WKN(I)*H)*SIN(WKX) 
                END DO 
              END IF             
! 
        ELSEIF(IORDER.EQ.2) THEN 
 
              IF(H .LE. 0.0D0) THEN 
                   DPOX=0.D0 
                  DPOY=0.D0 
                  DPOZ=0.D0 
                 
                IF (Nwave .ge. 2) THEN 
                   DO J=1,Nwave-1 
                  DO I=J+1,Nwave 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
                    DPOX=DPOX-FREQ(I)*AMPN(I)*AMPN(J)*WK_D*EXP(WK_D*Z)    & 
                                     *COS(BETA)*COS(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                    DPOY=DPOY-FREQ(I)*AMPN(I)*AMPN(J)*WK_D*EXP(WK_D*Z)    & 
                                     *SIN(BETA)*COS(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                    DPOZ=DPOZ-FREQ(I)*AMPN(I)*AMPN(J)*WK_D*EXP(WK_D*Z)    & 
                                    *SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                   END DO 
                   END DO 
                 END IF 
              ELSE 
                   DPOX=0.D0 
                   DPOY=0.D0 
                   DPOZ=0.D0 
                  DO I=1,Nwave 
                     DUM=3.0*FREQ(I)*AMPN(I)*AMPN(I)/8.0 
                    WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I)             
                    SH4=(SINH(WKN(I)*H))**4 
 
                     DPOX=DPOX+DUM*2.0*WKN(I)*COS(BETA)*COSH(2.0*WKN(I)   & 
                          *(Z+H))/SH4*COS(2.0*WKX) 
                     DPOY=DPOY+DUM*2.0*WKN(I)*SIN(BETA)*COSH(2.0*WKN(I)   & 
                           *(Z+H))/SH4*COS(2.0*WKX) 
                     DPOZ=DPOZ+DUM*2.0*WKN(I)*SINH(2.0*WKN(I)*(Z+H))/SH4*SIN(2.0*WKX) 
                END DO 
 
            if (nwave .ge. 2) then 
                 DO J=1,Nwave-1 
                  DO I=J+1,Nwave 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
 
                    DUM5=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_S    & 
                         *(1-TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM6=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2+   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2      
                    DUM7=G*WK_S*TANH(WK_S*H)-FREQ_S**2 
         
                    GD_S=-G*G*(DUM5+DUM6)/DUM7 
 
                    DUM5=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_D     & 
                         *(1+TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM6=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2-   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2  
                   DUM7=G*WK_D*TANH(WK_D*H)-FREQ_D**2           
         
                    GD_D=-G*G*(DUM5+DUM6)/DUM7 
 
                   DUM1=GD_S*COS(WK_S*RXY-FREQ_S*TIME+PHY_S)*COSH(WK_S*(Z+H))/COSH(WK_S*H) 
  
                   DUM2=GD_D*COS(WK_D*RXY-FREQ_D*TIME+PHY_D)*COSH(WK_D*(Z+H))/COSH(WK_D*H) 
 
                    DUM3=GD_S*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S)*SINH(WK_S*(Z+H))/COSH(WK_S*H) 
  
                    DUM4=GD_D*SIN(WK_D* RXY-FREQ_D*TIME+PHY_D)*SINH(WK_D*(Z+H))/COSH(WK_D*H) 
                           
                    DPOX=DPOX+AMPN(I)*AMPN(J)*COS(BETA)*(WK_S*DUM1+WK_D*DUM2) 
                    DPOY=DPOY+AMPN(I)*AMPN(J)*SIN(BETA)*(WK_S*DUM1+WK_D*DUM2) 
                    DPOZ=DPOZ+AMPN(I)*AMPN(J)*(WK_S*DUM3+WK_D*DUM4) 
                      END DO 
                    END DO 
                end if 
                END IF 
 
          ENDIF 
             DPOX=RAMP*DPOX 
             DPOY=RAMP*DPOY 
             DPOZ=RAMP*DPOZ 
 
        RETURN 
        END SUBROUTINE  DPOXYZ 
 
 
 
 
! 
!============================================================================ 
!  
!    Second spacial and time Derivatives of the first order component  
!    of incident wave potential 
!   d^2 phi/dt/dx,  d^2 phi/dt/dy,  d^2 phi/dt/dz  
!============================================================================ 
! 
! 
        SUBROUTINE DPDXYZT(H,G,Ampn,Phi,BETA,WKN,Freq,Time,Ramp,X,Y,Z,NN,Nwave,IOrder,DPDXDT,DPDYDT,DPDZDT) 
 
!DEC$ATTRIBUTES DLLEXPORT::DPDXYZT 
 
           IMPLICIT    NONE 
            
        INTEGER, INTENT(IN)::  NN,Nwave,IOrder 
            real*8,INTENT(IN)::   H,G,Ampn(NN),Phi(NN),BETA,WKN(NN),FREQ(NN),Time,Ramp,X,Y,Z 
        real*8,INTENT(OUT)::  DPDXDT,DPDYDT,DPDZDT 
 
            INTEGER I,J 
         
           real*8  RXY,WKX,SH4 
       real*8  WK_S,WK_D,FREQ_S,FREQ_D,PHY_S,PHY_D 
           real*8  DUM,DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7 
           real*8  GD_S,GD_D 
 
! 
! 
         RXY=X*COS(BETA)+Y*SIN(BETA) 
 
     IF(IORDER.EQ. 1) THEN 
            IF(H .LE. 0.0D0) THEN 
              DPDXDT=0.D0 
              DPDYDT=0.D0 
              DPDZDT=0.D0 
              DO I=1,Nwave 
                   DUM=AMPN(I)*G*WKN(I) 
                  WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                   DPDXDT=DPDXDT+DUM*COS(BETA)*EXP(WKN(I)*Z)*SIN(WKX) 
                   DPDYDT=DPDYDT+DUM*SIN(BETA)*EXP(WKN(I)*Z)*SIN(WKX)   
                  DPDZDT=DPDZDT-DUM*EXP(WKN(I)*Z)*COS(WKX) 
             END DO 
            ELSE 
              DPDXDT=0.D0 
              DPDYDT=0.D0 
              DPDZDT=0.D0 
              DO I=1,Nwave 
                    DUM=AMPN(I)*G*WKN(I)        
                   WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
                    DPDXDT=DPDXDT+DUM*COS(BETA)*COSH(WKN(I)*(Z+H))/  & 
                                  COSH(WKN(I)*H)*SIN(WKX) 
                    DPDYDT=DPDYDT+DUM*SIN(BETA)*COSH(WKN(I)*(Z+H))/   & 
                                  COSH(WKN(I)*H)*SIN(WKX)   
                   DPDZDT=DPDZDT-DUM*SINH(WKN(I)*(Z+H))/   & 
                                 COSH(WKN(I)*H)*COS(WKX) 
              END DO 
            END IF 
           
        ELSE IF(IORDER.EQ.2) THEN 
               IF(H .LE. 0.0D0) THEN 
                   DPDXDT=0.D0 
                  DPDYDT=0.D0 
                  DPDZDT=0.D0 
                 IF (Nwave .ge. 2) THEN 
                   DO J=1,Nwave-1 
                  DO I=J+1,Nwave 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
                    DPDXDT=DPDXDT-FREQ(I)*AMPN(I)*AMPN(J)*WK_D*FREQ_D*EXP(WK_D*Z)    & 
                                       *COS(BETA)*sin(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                    DPDYDT=DPDYDT-FREQ(I)*AMPN(I)*AMPN(J)*WK_D*FREQ_D*EXP(WK_D*Z)    & 
                                       *SIN(BETA)*sin(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                    DPDZDT=DPDZDT+FREQ(I)*AMPN(I)*AMPN(J)*WK_D*FREQ_D*EXP(WK_D*Z)    & 
                                    *cos(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                   END DO 
                   END DO 
                 END IF 
              ELSE 
                    DPDXDT=0.D0 
                   DPDYDT=0.D0 
                   DPDZDT=0.D0 
                  DO I=1,Nwave 
                     DUM=3.0*FREQ(I)*AMPN(I)*AMPN(I)/8.0 
                    WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I)             
                    SH4=(SINH(WKN(I)*H))**4 
 
                     DPDXDT=DPDXDT+DUM*4.0*WKN(I)*COS(BETA)*FREQ(I)*COSH(2.0*WKN(I)   & 
                          *(Z+H))/SH4*SIN(2.0*WKX) 
                     DPDYDT=DPDYDT+DUM*4.0*WKN(I)*SIN(BETA)*FREQ(I)*COSH(2.0*WKN(I)   & 
                           *(Z+H))/SH4*SIN(2.0*WKX) 
                     DPDZDT=DPDZDT+DUM*4.0*WKN(I)*FREQ(I)*SINH(2.0*WKN(I)*(Z+H))/SH4*COS(2.0*WKX) 
                END DO 
 
             if (nwave .ge. 2) then 
                 DO J=1,Nwave-1 
                  DO I=J+1,Nwave 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
 
                    DUM5=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_S    & 
                         *(1-TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM6=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2+   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2      
                    DUM7=G*WK_S*TANH(WK_S*H)-FREQ_S**2 
         
                    GD_S=-G*G*(DUM5+DUM6)/DUM7 
 
                    DUM5=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_D     & 
                         *(1+TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM6=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2-   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2  
                   DUM7=G*WK_D*TANH(WK_D*H)-FREQ_D**2           
         
                    GD_D=-G*G*(DUM5+DUM6)/DUM7 
 
                   DUM1=GD_S*FREQ_S*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S)*COSH(WK_S*(Z+H))/COSH(WK_S*H) 
  
                   DUM2=GD_D*FREQ_D*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D)*COSH(WK_D*(Z+H))/COSH(WK_D*H) 
 
                    DUM3=-GD_S*FREQ_S*COS(WK_S*RXY-FREQ_S*TIME+PHY_S)*SINH(WK_S*(Z+H))/COSH(WK_S*H) 
  
                    DUM4=-GD_D*FREQ_D*SIN(WK_D* RXY-FREQ_D*TIME+PHY_D)*SINH(WK_D*(Z+H))/COSH(WK_D*H) 
                           
                    DPDXDT=DPDXDT+AMPN(I)*AMPN(J)*COS(BETA)*(WK_S*DUM1+WK_D*DUM2) 
                    DPDYDT=DPDYDT+AMPN(I)*AMPN(J)*SIN(BETA)*(WK_S*DUM1+WK_D*DUM2) 
                    DPDZDT=DPDZDT+AMPN(I)*AMPN(J)*(WK_S*DUM3+WK_D*DUM4) 
                      END DO 
                    END DO 
                end if 
                END IF         
 
        END IF 
 
         DPDXDT=RAMP*DPDXDT 
         DPDYDT=RAMP*DPDYDT 
         DPDZDT=RAMP*DPDZDT 
 
         RETURN 
         END SUBROUTINE DPDXYZT        
 
! 
!============================================================================ 
!  
!   SECOND ORDER Spacial Derivatives of FIRST ORDER incident wave potential 
! 
!   d^2 phi/dx/dx,  d^2 phi/dy/dy,  d^2 phi/dz/dz  
!   d^2 phi/dx/dy,  d^2 phi/dx/dz,  d^2 phi/dy/dz  
! 
!============================================================================ 
 
        SUBROUTINE DDINP(H,G,Ampn,Phi,BETA,WKN,Freq,Time,Ramp,X,Y,Z,NN,Nwave,IOrder,   & 
                          DPOXX,DPOXY,DPOXZ,DPOYY,DPOYZ,DPOZZ) 
 
!DEC$ATTRIBUTES DLLEXPORT::DDINP 
 
           IMPLICIT    NONE 
            
        INTEGER,INTENT(IN)::  NN,Nwave,IOrder 
            real*8,INTENT(IN)::   H,G,Ampn(NN),Phi(NN),BETA,WKN(NN),FREQ(NN),Time,Ramp,X,Y,Z 
        real*8,INTENT(OUT)::  DPOXX,DPOXY,DPOXZ,DPOYY,DPOYZ,DPOZZ 
 
            INTEGER I,J 
         
 
            real*8  RXY,WKX,SH4 
        real*8  WK_S,WK_D,FREQ_S,FREQ_D,PHY_S,PHY_D 
            real*8  DUM,DUM1,DUM2,DUM3,DUM4,DUM5,DUM6,DUM7,GD_S,GD_D 
 
! 
         RXY=X*COS(BETA)+Y*SIN(BETA) 
 
              DPOXX=0.D0 
              DPOXY=0.D0 
              DPOXZ=0.D0 
              DPOYY=0.D0 
              DPOYZ=0.D0 
              DPOZZ=0.D0 
 
          IF(IORDER.EQ. 1) THEN 
 
             IF(H .LE. 0.0D0) THEN 
 
              DO I=1,Nwave 
                DUM=G*AMPN(I)*WKN(I)*WKN(I)/FREQ(I) 
                WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I) 
 
                DPOXX=DPOXX-DUM*COS(BETA)*COS(BETA)*EXP(WKN(I)*Z)*SIN(WKX) 
                DPOXY=DPOXY-DUM*SIN(BETA)*COS(BETA)*EXP(WKN(I)*Z)*SIN(WKX) 
               DPOXZ=DPOXZ+DUM*COS(BETA)*EXP(WKN(I)*Z)*COS(WKX) 
               DPOYY=DPOYY-DUM*SIN(BETA)*SIN(BETA)*EXP(WKN(I)*Z)*SIN(WKX) 
               DPOYZ=DPOYZ+DUM*SIN(BETA)*EXP(WKN(I)*Z)*COS(WKX) 
               DPOZZ=DPOZZ+DUM*EXP(WKN(I)*Z)*SIN(WKX) 
              END DO 
             
             ELSE 
 
              DO I=1,Nwave 
                DUM=G*AMPN(I)*WKN(I)*WKN(I)/FREQ(I) 
                WKX=WKN(I)*RXY -FREQ(I)*Time+Phi(I) 
 
                DPOXX=DPOXX-DUM*COS(BETA)*COS(BETA)*COSH(WKN(I)*(Z+H))/    & 
                      COSH(WKN(I)*H)*SIN(WKX) 
                DPOXY=DPOXY-DUM*SIN(BETA)*COS(BETA)*COSH(WKN(I)*(Z+H))/    & 
                      COSH(WKN(I)*H)*SIN(WKX) 
               DPOXZ=DPOXZ+DUM*COS(BETA)*SINH(WKN(I)*(Z+H))/       & 
                      COSH(WKN(I)*H)*COS(WKX) 
               DPOYY=DPOYY-DUM*SIN(BETA)*SIN(BETA)*COSH(WKN(I)*(Z+H))/     & 
                      COSH(WKN(I)*H)*SIN(WKX) 
               DPOYZ=DPOYZ+DUM*SIN(BETA)*SINH(WKN(I)*(Z+H))/     & 
                      COSH(WKN(I)*H)*COS(WKX) 
               DPOZZ=DPOZZ+DUM*COSH(WKN(I)*(Z+H))/COSH(WKN(I)*H)*SIN(WKX) 
              END DO 
            END IF 
  
            ELSE IF(IORDER.EQ.2) THEN 
 
                 IF(H .LE. 0.0D0) THEN 
                   DPOXX=0.D0 
                  DPOXY=0.D0 
                   DPOXZ=0.D0 
                  DPOYY=0.D0 
                   DPOYZ=0.D0 
                  DPOZZ=0.D0 
                  
                 IF (Nwave .ge. 2) THEN 
                   DO J=1,Nwave-1 
                     DO I=J+1,Nwave 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
                    DPOXX=DPOXX+FREQ(I)*AMPN(I)*AMPN(J)*WK_D**2*EXP(WK_D*Z)    & 
                                     *COS(BETA)**2*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                   DPOXY=DPOXY+FREQ(I)*AMPN(I)*AMPN(J)*WK_D**2*EXP(WK_D*Z)    & 
                                     *COS(BETA)*SIN(BETA)*COS(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                   DPOXZ=DPOXZ-FREQ(I)*AMPN(I)*AMPN(J)*WK_D**2*EXP(WK_D*Z)    & 
                                    *COS(BETA)*COS(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                    DPOYY=DPOYY+FREQ(I)*AMPN(I)*AMPN(J)*WK_D**2*EXP(WK_D*Z)    & 
                                     *SIN(BETA)**2*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                    DPOYZ=DPOYZ-FREQ(I)*AMPN(I)*AMPN(J)*WK_D**2*EXP(WK_D*Z)    & 
                                     *SIN(BETA)*COS(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                    DPOZZ=DPOZZ-FREQ(I)*AMPN(I)*AMPN(J)*WK_D**2*EXP(WK_D*Z)    & 
                                    *SIN(WK_D*RXY-FREQ_D*TIME+PHY_D) 
                   END DO 
                   END DO 
                 END IF 
              ELSE 
                    DPOXX=0.D0 
                   DPOXY=0.D0 
                   DPOXZ=0.D0 
                    DPOYY=0.D0 
                   DPOYZ=0.D0 
                   DPOZZ=0.D0 
 
                  DO I=1,Nwave 
                     DUM=3.0*FREQ(I)*AMPN(I)*AMPN(I)/8.0 
                    WKX=WKN(I)*RXY-FREQ(I)*Time+Phi(I)             
                    SH4=(SINH(WKN(I)*H))**4 
 
                     DPOXX=DPOXX-DUM*4.0*WKN(I)**2*COS(BETA)**2*COSH(2.0*WKN(I)   & 
                          *(Z+H))/SH4*SIN(2.0*WKX) 
                    DPOXY=DPOXZ-DUM*4.0*WKN(I)**2*COS(BETA)*SIN(BETA)*COSH(2.0*WKN(I)   & 
                          *(Z+H))/SH4*SIN(2.0*WKX) 
                    DPOXZ=DPOXZ+DUM*4.0*WKN(I)**2*COS(BETA)*SINH(2.0*WKN(I)   & 
                          *(Z+H))/SH4*COS(2.0*WKX) 
                     DPOYY=DPOYY-DUM*4.0*WKN(I)**2*SIN(BETA)**2*COSH(2.0*WKN(I)   & 
                           *(Z+H))/SH4*SIN(2.0*WKX) 
                    DPOYZ=DPOYZ+DUM*4.0*WKN(I)**2*SIN(BETA)*SINH(2.0*WKN(I)   & 
                           *(Z+H))/SH4*COS(2.0*WKX) 
                     DPOZZ=DPOZZ+DUM*4.0*WKN(I)**2*COSH(2.0*WKN(I)*(Z+H))/SH4*SIN(2.0*WKX) 
                END DO 
               
              IF (Nwave .ge. 2) THEN 
                 DO J=1,Nwave-1 
                  DO I=J+1,Nwave 
                    WK_S=WKN(I)+WKN(J) 
                    WK_D=WKN(I)-WKN(J) 
                    FREQ_S=FREQ(I)+FREQ(J) 
                   FREQ_D=FREQ(I)-FREQ(J) 
                    PHY_S=Phi(I)+Phi(J) 
                    PHY_D=Phi(I)-Phi(J) 
 
                    DUM5=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_S    & 
                         *(1-TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM6=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2+   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2      
                    DUM7=G*WK_S*TANH(WK_S*H)-FREQ_S**2 
         
                    GD_S=-G*G*(DUM5+DUM6)/DUM7 
 
                    DUM5=WKN(I)*WKN(J)/FREQ(I)/FREQ(J)*FREQ_D     & 
                         *(1+TANH(WKN(I)*H)*TANH(WKN(J)*H)) 
                    DUM6=0.5*WKN(I)*WKN(I)/FREQ(I)/(COSH(WKN(I)*H))**2-   & 
                         0.5*WKN(J)*WKN(J)/FREQ(J)/(COSH(WKN(J)*H))**2  
                   DUM7=G*WK_D*TANH(WK_D*H)-FREQ_D**2           
         
                    GD_D=-G*G*(DUM5+DUM6)/DUM7 
 
                   DUM1=-GD_S*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S)*COSH(WK_S*(Z+H))/COSH(WK_S*H) 
  
                   DUM2=-GD_D*SIN(WK_D*RXY-FREQ_D*TIME+PHY_D)*COSH(WK_D*(Z+H))/COSH(WK_D*H) 
 
                   DUM3=GD_S*SIN(WK_S*RXY-FREQ_S*TIME+PHY_S)*COSH(WK_S*(Z+H))/COSH(WK_S*H) 
  
                    DUM4=GD_D*SIN(WK_D* RXY-FREQ_D*TIME+PHY_D)*COSH(WK_D*(Z+H))/COSH(WK_D*H) 
 
                   DUM5=GD_S*COS(WK_S*RXY-FREQ_S*TIME+PHY_S)*COSH(WK_S*(Z+H))/COSH(WK_S*H) 
  
                   DUM6=GD_D*COS(WK_D*RXY-FREQ_D*TIME+PHY_D)*COSH(WK_D*(Z+H))/COSH(WK_D*H) 
 
 
 
                   DPOXX=DPOXX+AMPN(I)*AMPN(J)*COS(BETA)**2*(WK_S**2*DUM1+WK_D**2*DUM2) 
                   DPOXY=DPOXY+AMPN(I)*AMPN(J)*COS(BETA)*SIN(BETA)*(WK_S**2*DUM1+WK_D**2*DUM2) 
                   DPOXZ=DPOXZ+AMPN(I)*AMPN(J)*COS(BETA)*(WK_S**2*DUM5+WK_D**2*DUM6) 
                    DPOYY=DPOYY+AMPN(I)*AMPN(J)*SIN(BETA)**2*(WK_S**2*DUM1+WK_D**2*DUM2) 
                    DPOYZ=DPOYZ+AMPN(I)*AMPN(J)*SIN(BETA)*(WK_S**2*DUM5+WK_D**2*DUM6) 
                    DPOZZ=DPOZZ+AMPN(I)*AMPN(J)*(WK_S**2*DUM3+WK_D**2*DUM4)             
                     
                      END DO 
                    END DO 
                 END IF 
 
                END IF 
 
            END IF 
 
            DPOXX=RAMP*DPOXX 
            DPOXY=RAMP*DPOXY 
            DPOXZ=RAMP*DPOXZ 
            DPOYY=RAMP*DPOYY 
            DPOYZ=RAMP*DPOYZ 
            DPOZZ=RAMP*DPOZZ 
 
        RETURN 
        END SUBROUTINE DDINP
end module
