    !
    ! *************************************************************
    ! *                                                           *
    ! *  The source point is in the mesh                          *
    ! *                                                           *
    ! *************************************************************
    !
    SUBROUTINE SGBD0_1(IS,IELEM,NODN,XP,YP,ZP,VALG,VALDG) 
        USE MVAR_MOD
        USE TRVar_mod    
        USE  MFUNC_mod   
        use proj_cnst,only:ex,ey,xiqet,xiqsi
        use green_funcs,only:gcombo1

        IMPLICIT NONE

        INTEGER IS,IELEM,N,J,IP,NODN   
        integer loop1,loop2,i,NSAMB
        INTEGER INODD,LI,LJ,LK,INODE

        REAL*8  XP,YP,ZP,X,Y,Z,r       
        REAL*8  DPOX,DPOY,DPOZ,DGN
        REAL*8  VALG(8),VALDG(8),GXF(4)
        REAL*8 XITSI(6),XITET(6),SI,ETA,DUMX
        real*8  Xiq(8),Wiq(8),Xit(7),EtaT(7),WIT(7)
        REAL*8  SF0(8),DSF0(2,8),DDSF0(3,8) 
        REAL*8  SF(8),DSF(2,8),DDSF(3,8),DET1,DET2,DET3,DET 
        REAL*8  SISM,ETASM,FITG,SI1,ETA1
        REAL*8  XJ(2,3),XJp(2,3),XJSM(2,3),XXJ(3,3)
        REAL*8  jk0(3),jk1c(3),jk1s(3),jk1(3),n0,n1c,n1s
        REAL*8  JKt(3)
        REAL*8  DX,DY,DZ,nx,ny,nz
        real*8  a1,a2,a3,b1,b2,b3,ast,bst,tot,totF,totf_1,totf_2
        real*8  n1,b30,a30,b31,a31,g31,f,f1,f2,gv,xff(13)    
        REAL*8  CSST,SNST,PLO,BETAs,GAMMA
        REAL*8  AS_3,AS_2
        REAL*8  Line_Ele,Area0_Ele,Area1_Ele,Area2_Ele
        REAL*8  GXF0(4),DGN0

        DATA XIQ/ 0.960289856497536D+00, 0.796666477413626D+00, &
            0.525532409916329D+00, 0.183434642495650D+00, &
            -0.183434642495650D+00,-0.525532409916329D+00, &
            -0.796666477413626D+00,-0.960289856497536D+00/

        DATA WIQ/ 0.101228536290376D+00, 0.222381034453374D+00, &
            0.313706645877887D+00, 0.362683783378362D+00, &
            0.362683783378362D+00, 0.313706645877887D+00, &
            0.222381034453374D+00, 0.101228536290376D+00/     
        !       
        !    ============================================    
        valg=0.0d0 
        valdg=0.0d0
        line_ele=0.0d0
        Area0_Ele=0.0d0
        Area1_Ele=0.0d0
        Area2_Ele=0.0d0
        INODD=NCOND(IELEM,NODN)

        IF(NCN(IELEM).EQ.8)  THEN 
            !               
            SI =XIQSI(NODN)
            ETA=XIQET(NODN)
            CALL SPFUNC8_1(SI,ETA,SF0,DSF0,DDSF0) 

        ENDIF
        !       
        ! ----------------------------------------------  (C2)          
        !
        DO    LI=1,2
            DO    LJ=1,3
                DUMX=0.0D0
                DO   LK=1, NCN(IELEM)
                    DUMX=DUMX+DSF0(LI,LK)*XYZ(LJ,NCON(IELEM,LK))
                ENDDO
                XJ(LI,LJ)=DUMX
            enddo
        enddo!

        !
        ! ** compute the determinant of the Jacobian matrix at (SI,ETA), DET
        !
        ! **计算Jk0                         ------ (C12)

        JK0(1)=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2)    !  J1
        JK0(2)=XJ(1,3)*XJ(2,1)-XJ(1,1)*XJ(2,3)    !  J2
        JK0(3)=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)    !  J3
        !  


        !  --------------------------------------------  (C3)
        !
        DO   LI=1,3
            DO   LJ=1,3
                DUMX=0.0D0
                DO   LK=1,  NCN(IELEM)
                    DUMX=DUMX+DDSF0(LI,LK)*XYZ(LJ,NCON(IELEM,LK))   
                ENDDO
                XXJ(LI,LJ)=DUMX
            end do
        enddo

        ! **计算  JK1=Jk1C*cos（theta）+Jk1S*sin（theta） ------ (C12)

        JK1C(1)=(XXJ(1,2)*XJ(2,3)+XJ(1,2)*XXJ(3,3)) -   &
            (XXJ(1,3)*XJ(2,2)+XJ(1,3)*XXJ(3,2))
        JK1C(2)=(XXJ(1,3)*XJ(2,1)+XJ(1,3)*XXJ(3,1)) -   &
            (XXJ(1,1)*XJ(2,3)+XJ(1,1)*XXJ(3,3))
        JK1C(3)=(XXJ(1,1)*XJ(2,2)+XJ(1,1)*XXJ(3,2)) -  &
            (XXJ(1,2)*XJ(2,1)+XJ(1,2)*XXJ(3,1))

        JK1S(1)=(XXJ(3,2)*XJ(2,3)+XJ(1,2)*XXJ(2,3)) -   &
            (XXJ(3,3)*XJ(2,2)+XJ(1,3)*XXJ(2,2))
        JK1S(2)=(XXJ(3,3)*XJ(2,1)+XJ(1,3)*XXJ(2,1)) -   &
            (XXJ(3,1)*XJ(2,3)+XJ(1,1)*XXJ(2,3))
        JK1S(3)=(XXJ(3,1)*XJ(2,2)+XJ(1,1)*XXJ(2,2)) -  &
            (XXJ(3,2)*XJ(2,1)+XJ(1,2)*XXJ(2,1))    

        ! =========================================================
        !

        DO  J=1, NCN(IELEM)

            ! **计算N0，N1=N1C*cos（theta）+N1S*sin（theta）   (C13)

            N0 =SF0(J)
            N1C=DSF0(1,J)
            N1S=DSF0(2,J)
            !

            CALL CIRBOD_2(NODN,NCN(IELEM),SI,ETA,FITG,   &
                JK0,JK1C,JK1S,N0,N1C,N1S,XJ,XXJ)
            
            !line_ele=line_ele+fitg
            valdg(j) = fitg


        enddo

        !


        NSAMB=0

        IF(NCN(IELEM).EQ.8)  THEN 

            DO  LOOP1=1, 8 
                DO  LOOP2=1, 8      
                    NSAMB=NSAMB+1  

                    SISM=XIQ(loop1)
                    ETASM=XIQ(loop2)  
                    PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
                    CSST=(SISM-SI)/PLO 
                    SNST=(ETASM-ETA)/PLO        

                    CALL SPFUNC8_1(SISM,ETASM,SF,DSF,DDSF)
                    !
                    X=0.0D0
                    Y=0.0D0
                    Z=0.0D0
                    NX=0.0D0
                    NY=0.0D0
                    NZ=0.0D0
                    DO  LK=1,  NCN(IELEM)
                        X=X+SF(LK)*XYZ(1,NCON(IELEM,LK))    
                        Y=Y+SF(LK)*XYZ(2,NCON(IELEM,LK))
                        Z=Z+SF(LK)*XYZ(3,NCON(IELEM,LK))
                        NX=NX+SF(LK)*DXYZ(1,NCOND(IELEM,LK))    
                        NY=NY+SF(LK)*DXYZ(2,NCOND(IELEM,LK))
                        NZ=NZ+SF(LK)*DXYZ(3,NCOND(IELEM,LK))
                    END DO

                    !CALL DTGRN0(-1.0d0,X,XP,Y,YP,Z,ZP,GXF0)
                    call gcombo1(-1.0d0,(/x,y,z/),(/xp,yp,zp/),gxf0)
                    dgn0=gxf0(2)*nx+gxf0(3)*ny+gxf0(4)*nz

                    DO    LI=1,2
                        DO    LJ=1,3
                            DUMX=0.0D0
                            DO    LK=1,NCN(IELEM)
                                DUMX=DUMX+DSF(LI,LK)*XYZ(LJ,NCON(IELEM,LK))
                            enddo
                            XJP(LI,LJ)=DUMX
                        enddo
                    enddo
                    !
                    DET1=XJP(1,2)*XJP(2,3)-XJP(1,3)*XJP(2,2)    !  J1
                    DET2=XJP(1,3)*XJP(2,1)-XJP(1,1)*XJP(2,3)    !  J2
                    DET3=XJP(1,1)*XJP(2,2)-XJP(1,2)*XJP(2,1)    !  J3
                    DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)


                    JKt(1)=Jk0(1)+PLO*(JK1C(1)*CSST+JK1S(1)*SNST)
                    JKt(2)=Jk0(2)+PLO*(JK1C(2)*CSST+JK1S(2)*SNST)
                    JKt(3)=Jk0(3)+PLO*(JK1C(3)*CSST+JK1S(3)*SNST)




                    DO  J=1, NCN(IELEM)

                        ! **计算N0，N1=N1C*cos（theta）+N1S*sin（theta）   (C13)

                        N0 =SF0(J)
                        N1C=DSF0(1,J)
                        N1S=DSF0(2,J)   

                        CALL AREA_COEF(CSST,SNST,JK0,JK1C,JK1S,   &
                            N0,N1C,N1S,XJ,XXJ,F1,F2)

                        F=F2/PLO/PLO+F1/PLO                       ! (C19)

                        TOT=F/PLO*WIQ(LOOP1)*WIQ(LOOP2)
                        TOTf_2=F2/PLO/PLO/PLO*WIQ(LOOP1)*WIQ(LOOP2)
                        TOTf_1=F1/PLO/PLO*WIQ(LOOP1)*WIQ(LOOP2)         

                        TOTF=DGN0*DET*WIQ(LOOP1)*WIQ(LOOP2)*SF(J)

                        Area0_Ele=Area0_Ele+TOTF
                        Area1_Ele=Area1_Ele+TOTf_1
                        Area2_Ele=Area2_Ele+TOTf_2
                        valdg(j)=valdg(j)+totf-tot
                        !valdg(j)=totf

                    enddo

                enddo
            enddo

            !write(*,'(8f14.8)') valdg
            !pause
            !
            ! ======================================================================
            !
            ! **** FOR TRIANGULAR ELEMENTS **********************
            ! 
        endif



        !Area1_sum=Area1_sum-Area1_Ele

        !Area2_sum=Area2_sum-Area2_Ele

        !Area0_sum=Area0_sum+Area0_Ele

        !line_sum=line_sum+line_ele

        END

        !
        !C
        !C **************************************************************
        !C *                                                            *
        !C *  Calculating the one-dimension integral part of the        *
        !C *  singular integral                                         * 
        !C *  The singularity is on the body surface                    *
        !C *     by M. Guiggiani and A. Gigante's Method (Dec. 1990)    *
        !C *                                                            *
        !C **************************************************************
        !C
        SUBROUTINE CIRBOD_2(NODN,NCNE,SI,ETA,FITG,   &
                JK0,JK1C,JK1S,N0,N1C,N1S,XJ,XXJ)
            !        USE MVAR_MOD
            IMPLICIT NONE 
            !
            INTEGER NEWNUM,NODN,IN,I,NCNE
            INTEGER NEW6(6),NEW8(8)
            !
            real*8 XIQSI(32),XIQET(32),Wiq1(17),WIQ2(25)
            REAL*8 XITSI(24),XITET(24),WIT1(9),WIT2(17)
            REAL*8 WIT11(9),WIT21(9),WIT31(9),WIT42(17),WIT52(17),WIT62(17)
            REAL*8 FITG,FITG1,FITG2
            REAL*8 SI,ETA,SISM,ETASM,SISM1,ETASM1
            REAL*8 CSST,SNST,PLO,BETAs,GAMMA
            REAL*8 A1,A2,A3,AST,B1,B2,B3,BST
            REAL*8 A30,B30,jk0(3),jk1c(3),jk1s(3),jk1(3),n0,n1,n1c,n1s
            REAL*8 A31,B31,AS_3,AS_2,PI4
            REAL*8 F1,F2,SUM1,SUM2,G31
            REAL*8 XJ(2,3),XXJ(3,3)
            REAL*8 det1,det2,det3,dumx,sf(8),dsf(2,8),ddsf(3,8)           
            !    
            ! ** SAMPLING POINTS AND WEIGHTING FACTORS FOR QUADRILATERAL ELEMENT
            !
            !   8 nodes at each side
            !       
            DATA NEW8/1, 5, 9, 13, 17, 21, 25, 29/
            DATA XIQSI/-1.d0,-.75d0,-.5d0,-.25d0,0.d0,.25d0,.5d0,.75d0,1.d0,   &
                1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,    &
                .75d0,.5d0,.25d0,0.d0,-.25d0,-.5d0,-.75d0,-1.d0,   &
                -1.d0, -1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0/
            DATA XIQET/-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,-1.d0,  &
                -.75d0,-.5d0,-.25d0,0.d0,.25d0,.5d0,.75d0,1.d0,    &
                1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,   &
                .75d0,.5d0,.25d0,0.d0,-.25d0,-.5d0,-.75d0/ 

            DATA WIQ1/0.062177497273381,0.122489331563432,0.117207837861905,  &
                0.109334472936971,0.099914322536495,0.089926749896239,  &
                0.080115342139031,0.070948527302082,0.066568163775824,  &
                0.070948527302082,0.080115342139031,0.089926749896239,  &
                0.099914322536495,0.109334472936971,0.117207837861905, &
                0.122489331563432,0.062177497273381/


            DATA WIQ2/0.122489331563432,0.231823804500403,0.199261222833210,  &
                0.160875277198321,0.126277137889030,0.098697779924940,  &
                0.077797413988515,0.062177497273381,0.080187721987975,  &
                0.109334472936971,0.117207837861905,0.122489331563432,  &
                0.124354994546761,0.122489331563432,0.117207837861905,  &
                0.109334472936971,0.080187721987975,0.062177497273381,  &
                0.077797413988515,0.098697779924940,0.126277137889030,  &
                0.160875277198321,0.199261222833210,0.231823804500403,  &
                0.122489331563432/
            !    
            ! ** SAMPLING POINTS AND WEIGHTING FACTORS FOR TRIANGULAR ELEMENT
            !   8 nodes at each side
            !
            DATA NEW6/1, 9, 17, 5, 13, 21/     
            DATA XITSI/0.0d0, 0.125d0, 0.25d0, 0.375d0, 0.5d0, 0.615d0, 0.75d0, 0.875d0,   & 
                1.0d0, 0.875d0, 0.75d0, 0.615d0, 0.5d0, 0.375d0, 0.25d0, 0.125d0,   &
                0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0  /
            DATA XITET/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,       &
                0.125d0, 0.25d0, 0.375d0, 0.5d0, 0.615d0, 0.75d0, 0.875d0, 1.0d0,   &
                0.875d0, 0.75d0, 0.615d0, 0.5d0, 0.375d0, 0.25d0, 0.125d0/


            DATA WIT11/0.070948527302082, 0.160875277198321, 0.199261222833210,  &
                0.231823804500403, 0.244978663126864, 0.231823804500403,  &
                0.199261222833210, 0.160875277198321, 0.070948527302082  /          
            DATA WIT21/0.033284081887912, 0.070948527302082, 0.080115342139031,&
                0.089926749896239, 0.099914322536495, 0.109334472936971,  &
                0.117207837861905, 0.122489331563432, 0.062177497273381 /          
            DATA WIT31/0.062177497273381, 0.122489331563432, 0.117207837861905,&
                0.109334472936971, 0.099914322536495, 0.089926749896239,&
                0.080115342139031, 0.070948527302082, 0.033284081887912  /          

            DATA WIT42/0.160875277198321, 0.392699081698724, 0.463647609000806,&
                0.392699081698724, 0.259573057123262, 0.160875277198321, &
                0.103748113217601, 0.070948527302082, 0.057127163980720, &
                0.062177497273381, 0.077797413988515, 0.098697779924940,&
                0.126277137889030, 0.160875277198321, 0.199261222833210,&
                0.231823804500403, 0.122489331563432  /
            DATA WIT52/0.070948527302082, 0.160875277198321, 0.199261222833210, &
                0.231823804500403, 0.244978750549644, 0.231823891923183, &
                0.199261222833210, 0.160875277198321, 0.141897054604164, &
                0.160875277198321, 0.199261222833210, 0.231823804500403, &
                0.244978663126864, 0.231823804500403, 0.199261222833210,&
                0.160875277198321, 0.070948527302082 /
            DATA WIT62/0.122489331563432, 0.231823804500403, 0.199261222833210, &
                0.160875277198321, 0.126277137889030, 0.098697779924940, &
                0.077797413988515, 0.062177497273381, 0.057127163980720, &
                0.070948527302082, 0.103748113217601, 0.160875277198321, &
                0.259573057123261, 0.392699081698724, 0.463647609000806, &
                0.392699081698724, 0.160875277198321 /

            DATA PI4/12.56637061435917295385d0/ 

            !
            FITG=0.0D0
            FITG1=0.0d0
            FITG2=0.0d0
            SUM1=0.0d0
            SUM2=0.0d0          

            WRITE(10,*)
            WRITE(10,*) ' Inside CIRBOD_2   NCNE=',NCNE
            WRITE(10,*) '   NODN=',NODN,' SI=',SI,' ETA=',ETA


            !  =====================================================================         
            ! **** FOR QUADRILATERIAL ELEMENTS  **************************
            !       
            IF (NCNE .EQ. 8) THEN 
                NEWNUM=NEW8(NODN) 
                !                                
                !
                IF(NODN==1 .OR. NODN==3 .OR. NODN==5 .OR. NODN==7)  THEN  
                    !            
                    DO  200   i=1, 17!two edges
                        !    WRITE(10,*) '  I=',I
                        IN=MOD(NEWNUM+I+6, 32) + 1 
                        SISM  = XIQSI(IN) 
                        ETASM = XIQET(IN)   
                        !
                        PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
                        CSST=(SISM-SI)/PLO 
                        SNST=(ETASM-ETA)/PLO 

                        JK1(1)=JK1C(1)*CSST+ JK1S(1)*SNST    !(C12)
                        JK1(2)=JK1C(2)*CSST+ JK1S(2)*SNST
                        JK1(3)=JK1C(3)*CSST+ JK1S(3)*SNST     

                        !       WRITE(10,*) ' JK1=',JK1

                        N1= N1C*CSST+ N1S*SNST               !(C13)

                        A1=XJ(1,1)*CSST+XJ(2,1)*SNST
                        A2=XJ(1,2)*CSST+XJ(2,2)*SNST
                        A3=XJ(1,3)*CSST+XJ(2,3)*SNST

                        B1=XXJ(1,1)*CSST*CSST*0.5+XXJ(3,1)*CSST*SNST +  &
                            XXJ(2,1)*SNST*SNST*0.5
                        B2=XXJ(1,2)*CSST*CSST*0.5+XXJ(3,2)*CSST*SNST +  &
                            XXJ(2,2)*SNST*SNST*0.5
                        B3=XXJ(1,3)*CSST*CSST*0.5+XXJ(3,3)*CSST*SNST +  &
                            XXJ(2,3)*SNST*SNST*0.5

                        AST=DSQRT(A1*A1 + A2*A2 + A3*A3)   
                        BST=DSQRT(B1*B1 + B2*B2 + B3*B3)        !(C7)

                        AS_3=1.0d0/AST**3
                        AS_2=-3.0d0*(A1*B1+A2*B2+A3*B3)/AST**5


                        G31=B1*JK0(1)+B2*JK0(2)+B3*JK0(3)+A1*JK1(1)+A2*JK1(2)+A3*JK1(3)
                        G31=(A3/AST/AST)*G31

                        B30=-JK0(3)
                        A30=B30*N0

                        B31=3*G31-JK1(3)
                        A31=B31*N0+B30*N1

                        !WRITE(*,*) 'jk0(3)=',jk0(3)

                        !       F2=AS_3*JK0(3)
                        !       F1=AS_2*JK0(3)+AS_3*JK1(3)     

                        F1=AS_2*A30+AS_3*A31     
                        F2=AS_3*A30

                        GAMMA=-(A1*B1+A2*B2+A3*B3)/AST**4
                        BETAs=1.0d0/AST

                        !        WRITE(10,*) ' F1=',F1,' GAMMA=',GAMMA
                        !        WRITE(10,*) ' PLO=',PLO,'  BETAs=',BETAs
                        !        WRITE(10,*) ' DLOG(PLO/BETAs)=',DLOG(PLO/BETAs)

                        FITG2=FITG2-F2*(GAMMA/BETAs**2+1.0/PLO)*WIQ1(I) 
                        FITG1=FITG1+F1*DLOG(PLO/BETAs)*WIQ1(I) 

                        !       WRITE(10,*) ' FITG1=',FITG1
                        !       WRITE(10,*)

                        ! --------------------------------------------
                        !
                        !**  线积分中sum1应为0， sum2应为2*Pi    

                        SUM1=SUM1+F1*WIQ1(I)                      !(25)
                        SUM2=SUM2+F2/BETAs*WIQ1(I)                !(26)
                        !       SUM1=SUM1+F1*DLOG(1.0/AST)*WIQ11(I)                      !(25)
                        !       SUM2=SUM2+F2*(-1.0*GAMMA/BETAs**2)*WIQ11(I)              !(26)
                        !
                        200  CONTINUE
                        !
                    ELSE IF(NODN==2 .OR. NODN==4 .OR. NODN==6 .OR. NODN==8)   THEN  

                        DO  400   i=1, 25!3edges
                            !      WRITE(10,*) '  I=',I
                            IN=MOD(NEWNUM+I+2, 32) + 1 
                            SISM  = XIQSI(IN) 
                            ETASM = XIQET(IN)   
                            !C
                            PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
                            CSST=(SISM-SI)/PLO 
                            SNST=(ETASM-ETA)/PLO 

                            JK1(1)=JK1C(1)*CSST+ JK1S(1)*SNST    !(C12)
                            JK1(2)=JK1C(2)*CSST+ JK1S(2)*SNST
                            JK1(3)=JK1C(3)*CSST+ JK1S(3)*SNST     

                            !       WRITE(10,*) ' JK1=',JK1

                            N1= N1C*CSST+ N1S*SNST               !(C13)

                            A1=XJ(1,1)*CSST+XJ(2,1)*SNST
                            A2=XJ(1,2)*CSST+XJ(2,2)*SNST
                            A3=XJ(1,3)*CSST+XJ(2,3)*SNST

                            B1=XXJ(1,1)*CSST*CSST*0.5+XXJ(3,1)*CSST*SNST +  &
                                XXJ(2,1)*SNST*SNST*0.5
                            B2=XXJ(1,2)*CSST*CSST*0.5+XXJ(3,2)*CSST*SNST +  &
                                XXJ(2,2)*SNST*SNST*0.5
                            B3=XXJ(1,3)*CSST*CSST*0.5+XXJ(3,3)*CSST*SNST +  &
                                XXJ(2,3)*SNST*SNST*0.5

                            AST=DSQRT(A1*A1 + A2*A2 + A3*A3)   
                            BST=DSQRT(B1*B1 + B2*B2 + B3*B3)        !(C7)

                            AS_3=1.0d0/AST**3
                            AS_2=-3.0d0*(A1*B1+A2*B2+A3*B3)/AST**5


                            G31=B1*JK0(1)+B2*JK0(2)+B3*JK0(3)+A1*JK1(1)+A2*JK1(2)+A3*JK1(3)
                            G31=(A3/AST/AST)*G31

                            B30=-JK0(3)
                            A30=B30*N0

                            B31=3*G31-JK1(3)
                            A31=B31*N0+B30*N1
                            !WRITE(*,*) 'jk0(3)=',jk0(3)


                            !        F2=AS_3*JK0(3)
                            !        F1=AS_2*JK0(3)+AS_3*JK1(3)     

                            F1=AS_2*A30+AS_3*A31     
                            F2=AS_3*A30

                            GAMMA=-(A1*B1+A2*B2+A3*B3)/AST**4
                            BETAs=1.0d0/AST

                            !       WRITE(10,*) ' F1=',F1,' GAMMA=',GAMMA
                            !       WRITE(10,*) ' PLO=',PLO,'  BETAs=',BETAs
                            !       WRITE(10,*) ' DLOG(PLO/BETAs)=',DLOG(PLO/BETAs)
                            FITG2=FITG2-F2*(GAMMA/BETAs**2+1.0/PLO)*WIQ2(I) 
                            FITG1=FITG1+F1*DLOG(PLO/BETAs)*WIQ2(I) 

                            !      WRITE(10,*) ' FITG1=',FITG1
                            !      WRITE(10,*)

                            ! --------------------------------------------
                            !
                            !**  线积分中sum1应为0， sum2应为2*Pi    

                            SUM1=SUM1+F1*WIQ2(I)                      !(25)
                            SUM2=SUM2+F2/BETAs*WIQ2(I)                !(26)
                            !       SUM1=SUM1+F1*DLOG(1.0/AST)*WIQ11(I)                      !(25)
                            !       SUM2=SUM2+F2*(-1.0*GAMMA/BETAs**2)*WIQ11(I)              !(26)
                            !
                            400  CONTINUE 


                        END IF
                        !      

                                !      

                            END IF             
                            !
                            !  =====================================================================         
                            !
                            FITG=(FITG2+FITG1)/PI4

                            END

                            !
                            ! ===============================================================
                            !
                            !  Determine F1 and F2 for area integration
                            !
                            ! ===============================================================
                            !
                            SUBROUTINE AREA_COEF(CSST,SNST,JK0,JK1C,JK1S,N0,N1C,N1S,XJ,XXJ,F1,F2)

                                IMPLICIT NONE 
                                !
                                REAL*8 CSST,SNST,BETAs,GAMMA
                                REAL*8 jk0(3),jk1c(3),jk1s(3),n0,n1c,n1s
                                !
                                REAL*8  jk1(3),n1       
                                REAL*8  XJ(2,3),XXJ(3,3)
                                real*8  a1,a2,a3,b1,b2,b3,ast,bst,tot,totF,totf_1,totf_2
                                real*8  b30,a30,b31,a31,g31,f,f1,f2,xff(13)    
                                REAL*8  AS_3,AS_2,PI4

                                DATA PI4/12.56637061435917295385d0/ 

                                !
                                ! ==================================================================
                                !
                                !       write(10,*) ' Inside AREA_COEF'
                                !           
                                JK1(1)=JK1C(1)*CSST+ JK1S(1)*SNST    !(C12)
                                JK1(2)=JK1C(2)*CSST+ JK1S(2)*SNST
                                JK1(3)=JK1C(3)*CSST+ JK1S(3)*SNST     

                                !        write(10,*) ' JK1=',JK1

                                N1= N1C*CSST+ N1S*SNST               !(C13)
                                !       
                                !   -----------------------------------------(C6)
                                !             
                                A1=XJ(1,1)*CSST+XJ(2,1)*SNST
                                A2=XJ(1,2)*CSST+XJ(2,2)*SNST
                                A3=XJ(1,3)*CSST+XJ(2,3)*SNST

                                B1=XXJ(1,1)*CSST*CSST*0.5+XXJ(3,1)*CSST*SNST +   &
                                    XXJ(2,1)*SNST*SNST*0.5
                                B2=XXJ(1,2)*CSST*CSST*0.5+XXJ(3,2)*CSST*SNST +   &
                                    XXJ(2,2)*SNST*SNST*0.5
                                B3=XXJ(1,3)*CSST*CSST*0.5+XXJ(3,3)*CSST*SNST +   &
                                    XXJ(2,3)*SNST*SNST*0.5

                                AST=DSQRT(A1*A1 + A2*A2 + A3*A3)         !  (C7)  
                                BST=DSQRT(B1*B1 + B2*B2 + B3*B3)

                                AS_3=1.0d0/AST**3                        !  (C10)
                                AS_2=-3.0d0*(A1*B1+A2*B2+A3*B3)/AST**5

                                G31=B1*JK0(1)+B2*JK0(2)+B3*JK0(3)+A1*JK1(1)+A2*JK1(2)+A3*JK1(3)
                                G31=(A3/AST/AST)*G31                    !  (C16)

                                B30=-JK0(3)                             !  (C17)            
                                B31=3*G31-JK1(3)

                                A30=B30*N0                              !  (C18)
                                A31=B31*N0+B30*N1

                                F1=(AS_2*A30+AS_3*A31)/PI4     
                                F2=AS_3*A30/PI4          
                                !                         
                                END
