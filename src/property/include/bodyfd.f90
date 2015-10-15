!C
!C *********************************************************************
!C *                                                                   *
!C *  Identify the Gaussian sampling points and evaluate the           *
!C *  corresponding values for the shape functions and Jacobian        *
!C *  matrices                                                         *
!C *                                                                   *
!C *********************************************************************
!C
        subroutine get_gaussian_data(xc,yc,zc)
        use MFUNC_mod
        IMPLICIT   NONE  
        real(8),intent(in) :: xc,yc,zc
	  INTEGER IE,K,J,LK,LI,LJ,ISI,IETA,NSAMB,IND 
        REAL*8 XITSI(4),XITETA(4),WIT(4),XIQ(4),WIQ(4)
        REAL*8 SF(8),DSF(2,8),XJ(3,3) 

        REAL*8 DET,DET1,DET2,DET3,DUM,SI,ETA 

!C                      
!C
!C ** matrix XITSI store the Gauss-Legendre sampling points(3) for the 
!C    triangular element in SI coordinate
!C
        DATA XITSI/0.333333333333333D0, 0.6D0, 0.2D0, 0.2D0/
!C
!C ** matrix XITETA store the Gauss-Legendre sampling points(3) for the 
!C    triangular element in ETA coordinate                         
!C
        DATA XITETA/0.333333333333333D0, 0.2D0, 0.6D0, 0.2D0/

!C
!C ** matrix WIT store the Gauss-Legendre weighting factors for the 
!C    triangular element in SI coordinate
!C
        DATA WIT/-0.28125D0,.260416666666667D0,.260416666666667D0,&
     &          .260416666666667D0/
!C
!C
!C ** matrix XIQ store the Gauss-Legendre sampling points(3) for the 
!C    quadrilateral element 
!c   
        DATA XIQ/-0.861136311594053D0,-0.339981043584856D0,&
     &          0.339981043584856D0, 0.861136311594053D0/ 
!C
!C ** matrix WIQ store the Gauss-Legendre weighting factors for the 
!C    quadrilateral element 
!c
        DATA WIQ/0.347854845137454D0,0.652145154862546D0,&
     &         0.652145154862546D0,0.347854845137454D0/
!c

        DO 500 IE=1, NELEM
!c	Print *,' I=',I
!C
        NSAMB=0
!c
!c       NSAMB:  codes of sampling points inside an element
!C
        IF(NCN(IE).EQ.8) THEN
!C
!C ** Quadrilateral element **
!C
        DO 110 ISI=1,4
        DO 120 IETA=1,4
        NSAMB=NSAMB+1
!c	Print *,' NSAMB=',NSAMB
!C
!C **  calculate the shape function at the sampling points
!C
        SI =XIQ(ISI)
        ETA=XIQ(IETA)

        CALL SPFUNC8(SI,ETA,SF,DSF)
!C
!C ** evaluate the Jacobian matrix at (SI,ETA),  XJ(2,3)
!C        
!c       LI: 1--SI,  2--ETA
!c       LJ: 1--X,   2--Y,   3--Z
!c
      DO 130 LI=1,2
      DO 130 LJ=1,3
      DUM=0.0D0
      DO 140 LK=1, NCN(IE)
      DUM=DUM+DSF(LI,LK)*XYZE(LJ,LK,IE)
!      DUM=DUM+DSF(LI,LK)*XYZ(LJ,NCON(IE,LK))  XYZE(3,J,IE)
140   CONTINUE
130   XJ(LI,LJ)=DUM
!C
!C ** compute the determinant of the Jacobian maxtix at (SI,ETA), DET
!C
      DET1=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2) 
      DET2=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1) 
      DET3=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1) 
      DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)
!C
!C ** transform the local coordinates of the sampling points to 
!C    global coordinates
!                                   
!c       SAMBXY=??
!c       I: code of the element,  LI:1--X, 2--Y, 3--Z  
!c       NCON:code of the node in the global mesh
!c       NSAMB:   
!c       SF: shape function              
!c
      DO 160 LI=1,3
      SAMBXY(IE,NSAMB,LI)=0.0D0
      DO 170 LK=1,NCN(IE)
      SAMBXY(IE,NSAMB,LI)=SAMBXY(IE,NSAMB,LI)+SF(LK)*XYZE(LI,LK,IE)    
170   CONTINUE
160   CONTINUE


      DO 180 LI=1,3
      DSAMB(IE,NSAMB,LI)=0.0D0
      DO 190 LK=1,NCN(IE)
      DSAMB(IE,NSAMB,LI)=DSAMB(IE,NSAMB,LI)+SF(LK)*DXYZE(LI,LK,IE)
190   CONTINUE
180   CONTINUE  

      DSAMB(IE,NSAMB,4)=(SAMBXY(IE,NSAMB,2)-YC)*DSAMB(IE,NSAMB,3)- &
     &      (SAMBXY(IE,NSAMB,3)-ZC)*DSAMB(IE,NSAMB,2)
      DSAMB(IE,NSAMB,5)=(SAMBXY(IE,NSAMB,3)-ZC)*DSAMB(IE,NSAMB,1)- &
     &      (SAMBXY(IE,NSAMB,1)-XC)*DSAMB(IE,NSAMB,3)
      DSAMB(IE,NSAMB,6)=(SAMBXY(IE,NSAMB,1)-XC)*DSAMB(IE,NSAMB,2)- &
     &      (SAMBXY(IE,NSAMB,2)-YC)*DSAMB(IE,NSAMB,1)
!C
!C ** calculate the free surface boundary condition*WIT*DET
!C
      SAMB(IE,NSAMB,0)=WIQ(ISI)*WIQ(IETA)*DET
      DO 20 J=1,NCN(IE)
      SAMB(IE,NSAMB,J)=SF(J)*WIQ(ISI)*WIQ(IETA)*DET
20    CONTINUE


120   CONTINUE
110   CONTINUE

      ELSE
!
!C ** for triangular element **
!C
      DO 210 J=1, 4
      NSAMB=NSAMB+1
!C
!C **  calculate the shape function at the sampling points
!C
      SI =XITSI(J)
      ETA=XITETA(J)
      CALL SPFUNC6(SI,ETA,SF,DSF)
!C
!C ** evaluate the Jacobian matrix at (SI,ETA)
!C
      DO 230 LI=1,2
      DO 230 LJ=1,3
      DUM=0.0D0
      DO 240  LK=1, NCN(IE)
      DUM=DUM+DSF(LI,LK)*XYZE(LJ,LK,IE)
240   CONTINUE
230   XJ(LI,LJ)=DUM
!C
!C ** compute the determinant of the Jacobian maxtix at (SI,ETA)
!C
      DET1=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2) 
      DET2=XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1) 
      DET3=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1) 
      DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)
!
!C ** transform the local coordinates of the sampling points to 
!C    global coordinates
!C
      DO 260 LI=1,3
      SAMBXY(IE,NSAMB,LI)=0.0D0
      DO 270 LK=1,NCN(IE)
      SAMBXY(IE,NSAMB,LI)=SAMBXY(IE,NSAMB,LI)+SF(LK)*XYZE(LI,LK,IE)
270   CONTINUE
260   CONTINUE

      DO 280 LI=1,3
      DSAMB(IE,NSAMB,LI)=0.0D0
      DO 290 LK=1,NCN(IE)
      DSAMB(IE,NSAMB,LI)=DSAMB(IE,NSAMB,LI)+SF(LK)*DXYZE(LI,LK,IE)
290   CONTINUE  
280   CONTINUE 
     
	 DSAMB(IE,NSAMB,4)=(SAMBXY(IE,NSAMB,2)-YC)*DSAMB(IE,NSAMB,3)- &
     &      (SAMBXY(IE,NSAMB,3)-ZC)*DSAMB(IE,NSAMB,2)
       DSAMB(IE,NSAMB,5)=(SAMBXY(IE,NSAMB,3)-ZC)*DSAMB(IE,NSAMB,1)- &
     &      (SAMBXY(IE,NSAMB,1)-XC)*DSAMB(IE,NSAMB,3)
       DSAMB(IE,NSAMB,6)=(SAMBXY(IE,NSAMB,1)-XC)*DSAMB(IE,NSAMB,2)- &
     &      (SAMBXY(IE,NSAMB,2)-YC)*DSAMB(IE,NSAMB,1)
!C                                                                
!C ** calculate the free surface boundary condition*WIT*DET   
!C
      SAMB(IE,NSAMB,0)=WIT(J)*DET
      DO 30 K=1,NCN(IE)
      SAMB(IE,NSAMB,K)=SF(K)*WIT(J)*DET
30    CONTINUE   



220   CONTINUE
210   CONTINUE

      ENDIF

500    CONTINUE


	 DO IND=1, NNODED
	  
	 DXYZ(4,IND)=(XYZ(2,IND)-YC)*DXYZ(3,IND)- &
     &             (XYZ(3,IND)-ZC)*DXYZ(2,IND)
       DXYZ(5,IND)=(XYZ(3,IND)-ZC)*DXYZ(1,IND)- &
     &             (XYZ(1,IND)-XC)*DXYZ(3,IND)
       DXYZ(6,IND)=(XYZ(1,IND)-XC)*DXYZ(2,IND)- &
     &             (XYZ(2,IND)-YC)*DXYZ(1,IND)
	ENDDO


       
      RETURN
      END
