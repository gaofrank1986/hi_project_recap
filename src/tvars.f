
C
C =========================================
C  Varibles used in the Tri-pole transform 
C                                                         
       MODULE TRVar_mod

         INTEGER NOSAMP
         REAL*8 XYNOD(3,50),DXYNOD(6,50),SAMNOD(50,0:8)

       END MODULE TRVar_mod

C
C  **************************************************
C  *                                                                                              *
C  *  This is a module for declaring variables      *
C  *      in B-spline expansion                     *
C  *                                                *
C  *      Aug. 10, 2002          by   Bin Teng      *
C  **************************************************                           
C

       MODULE  BVAR_mod
C
       INTEGER KBSPL,NB,NT,NA,NAA,NBU,NA1
       INTEGER MUSTA,MUEND       
C
C --------------------------------------------------------------------
C KBSPL : DEGREE (order-1) of the B-splines
C
C NB    : Number of given points 
C NT    : Maximum number of given points
C
C NBU   : Number of intervals in the u-direction 
C NA1     : Number of control factors    (NBU+KBSPL)
C NAA     : Number of control factors ??   (NBU+2*KBSPL+1)
C NA    : Maximum number of total control factors  
C -------------------------------------------------------------------
C
         PARAMETER (NT=500, NA=100, KBSPL=3)
         REAL*8 BXYZ(2,NT),BXYZ1(2,NT),UXYZ(NT)
C
C XYZ     : Catersian coordinates of given points
C UXYZ    : U coordinates of given points on the body surface
C
       REAL*8 BJU(NA),B(NA,4),DBJU(NA)
         REAL*8 XYZCON(NA,2),UTJ(NA),UXYZNE(NA)
C
C UTJ    : U coordinates of nodes and co-located points
C XYZCON : factors for body geometry expansion
C UXYZNE : U coordinates of controlling points 
C
       END MODULE BVAR_mod
C
