C 
C ****************************************************** 
C *                                                    * 
C *  Evaluate the wave force on a 3-D body             * 
C *                                                    * 
C ****************************************************** 
C 
        SUBROUTINE TFORCE 
 
          USE MVAR_MOD 
          USE PVAR_MOD 
          use time_mod
          use wave_func,only:dpot2
! 
        IMPLICIT   NONE 
! 
          INTEGER IP,IELEM,NSAMB,N,K 
        real(8)  EXY(4,2),XP,YP,ZP     
        real(8)  DDUM 
 
        !real(8),EXTERNAL::  DPOT2 

C 
        DATA EXY /1., 1.,-1.,-1.,  1.,-1.,-1., 1./ 
C               
        FORCEW=0.0D0 
C 
        DO 100  IELEM=NELEMF+1,  NELEM  
        NSAMB=16 
        IF(NCN(IELEM) .EQ.6 ) NSAMB=4 
C                     
        DO  100  IP=1,  NSYS    
 
        DO  80  N =1, NSAMB   
C 
        DDUM=0.0D0 
        DO 40   K=1,  NCN(IELEM)  
          DDUM=DDUM+ DPDT(NCON(IELEM,K),IP)*SAMB(IELEM,N,K)  
!         Print *,'  IELEM=',IELEM,'  K=',K,'  NCON=',NCON(IELEM,K) 
!         Print *,'  DPDT=',DPDT(NCON(IELEM,K),IP),' SAMB=',SAMB(IELEM,N,K) 
40      CONTINUE  
C 
!       pause 
 
        XP=EXY(IP,1)*SAMBXY(IELEM,N,1) 
        YP=EXY(IP,2)*SAMBXY(IELEM,N,2)   
        ZP=          SAMBXY(IELEM,N,3)       
! 
        DDUM=DDUM+DPOT2(H,G,Ampn,Phi_w,BETA,WKN,Freq,TimeRK,RampF, 
     1          XP,YP,ZP,NFreq,Nwave,IOrder)*SAMB(IELEM,N,0)  
!          
        
        FORCEW(1)=FORCEW(1)+DDUM* 
     1                        EXY(IP,1)*DSAMB(IELEM,N,1) 
        FORCEW(2)=FORCEW(2)+ DDUM* 
     1                        EXY(IP,2)*DSAMB(IELEM,N,2) 
 
        FORCEW(3)=FORCEW(3) + DDUM* DSAMB(IELEM,N,3)        
C 
        FORCEW(4)=FORCEW(4)+DDUM* 
     2               EXY(IP,2)* DSAMB(IELEM,N,4)  
        FORCEW(5)=FORCEW(5)+DDUM* 
     2               EXY(IP,1)* DSAMB(IELEM,N,5) 
        FORCEW(6)=FORCEW(6)+DDUM* 
     2               EXY(IP,2)* EXY(IP,1)*DSAMB(IELEM,N,6) 
 
!         Print *,'  DDUM=',DDUM,' DX=',DSAMB(IELEM,N,1),' EX=',EXY(IP,1) 
 
80      CONTINUE 
! 
!         Print *,' IELEM=',IELEM,' F1=',FORCEW(1),' F3=',FORCEW(3) 
 
100     CONTINUE 
! 
C 
          FORCEW(:)=RHO*FORCEW(:) 
           
!         Print *,  ' F1=',FORCEW(1),' F3=',FORCEW(3) 
           
!         pause 
 
        RETURN 
        END                       
 
 
