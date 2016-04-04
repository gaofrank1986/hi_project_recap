    !
    ! *************************************************************
    ! *                                                           *
    ! *  The source point is in the mesh                          *
    ! *                                                           *
    ! *************************************************************
    !
    SUBROUTINE SGBD0_1(IS,IELEM,NODN,XP,YP,ZP,VALG,VALDG) 
        use kinds
        use mesh,only:xyz,dxyz,ncn,ncon,ncond
        use shape_funcs
        use proj_cnst,only:ex,ey,xiqet,xiqsi
        use linalg,only:cross_product
        use green_funcs,only:gcombo1_1

        IMPLICIT NONE

        INTEGER IS,IELEM,N,J,IP,NODN   
        integer loop1,loop2,i,NSAMB
        INTEGER INODD,LI,LJ,LK,INODE

        REAL*8  XP,YP,ZP       
        REAL*8  VALG(8),VALDG(8)
        REAL*8 XITSI(6),XITET(6),SI,ETA,DUMX
        real*8  Xiq(8),Wiq(8),Xit(7),EtaT(7),WIT(7)
        REAL*8  SF0(8),DSF0(2,8),DDSF0(3,8) 
        REAL*8  SF(8),DSF(2,8),DET,ans
     
        REAL*8  XJ(2,3),XJp(2,3),XXJ(3,3)
        REAL*8  jk0(3),jk1c(3),jk1s(3),n0,n1c,n1s
        REAL*8  xxx(3,8),xxd(3,8),p(3),np(3)
        real*8  tot,totF
        real*8  n1,f,f1,f2,gv
        REAL*8 PLO
        REAL*8  Line_Ele,Area0_Ele,Area1_Ele,Area2_Ele
        REAL*8  GXF0(4),DGN0,xi0(2),n01(3),xi(2),p0(3),jk01(3,3),param(2)

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
        p0=[xp,yp,zp]
        inodd=ncond(ielem,nodn)

        if(ncn(ielem).eq.8)  then 
            !               
            xi0=[xiqsi(nodn),xiqet(nodn)]
            si =xiqsi(nodn)
            eta=xiqet(nodn)
            call spfunc8_1(xi0(1),xi0(2),sf0,dsf0,ddsf0) 

            do lk =1,8
                xxx(1:3,lk) = xyz(1:3,ncon(ielem,lk))
                xxd(1:3,lk) = dxyz(1:3,ncond(ielem,lk))
            end do


        endif

        xj(1:2,1:3) = matmul(dsf0(1:2,1:8),transpose(xxx(1:3,1:8)))
        xxj(1:3,1:3) = matmul(ddsf0(1:3,1:8),transpose(xxx))

        jk0(1:3) = cross_product(xj(1,:),xj(2,:))
        jk1c=cross_product(xxj(1,:),xj(2,:))+cross_product(xj(1,:),xxj(3,:))
        jk1s=cross_product(xxj(3,:),xj(2,:))+cross_product(xj(1,:),xxj(2,:))

            jk01(1,:) = jk0
            jk01(2,:) = jk1c
            jk01(3,:) = jk1s


        DO  J=1, NCN(IELEM)
            n01=[sf0(j),dsf0(1,j),dsf0(2,j)]

            call cirbod_2(nodn,ncn(ielem),si,eta,   &
                jk01,n01,xj,xxj,ans)

            valdg(j) = ans
        enddo


        NSAMB=0

        IF(NCN(IELEM).EQ.8)  THEN 

            do  loop1=1, 8 ;do  loop2=1, 8      
                nsamb=nsamb+1  

                xi=[xiq(loop1),xiq(loop2)]

                call spfunc8(xi(1),xi(2),sf,dsf)
                plo=norm2(xi-xi0)
                xi=(xi-xi0)/plo

                p=matmul(sf,transpose(xxx))
                np=matmul(sf,transpose(xxd))


                call gcombo1_1(-1.0d0,p,[xp,yp,zp],gxf0)
                dgn0=dot_product(gxf0(2:4),np)

                xjp(1:2,1:3)=matmul(dsf,transpose(xxx))
                det =norm2(cross_product(xjp(1,:),xjp(2,:)))

                do  j=1, ncn(ielem)
                    n01=[sf0(j),dsf0(1,j),dsf0(2,j)]
                    call comp_coef(xi(1),xi(2),xj,xxj,jk01,n01,f1,f2,param)


                    f=f2/plo/plo+f1/plo                       ! (c19)

                    tot=f/plo*wiq(loop1)*wiq(loop2)
                    totf=dgn0*det*wiq(loop1)*wiq(loop2)*sf(j)

                    valdg(j)=valdg(j)+totf-tot

                enddo

            enddo; enddo

            !write(*,'(8f14.8)') valdg
            !pause
            !
            ! ======================================================================
            !
            ! **** FOR TRIANGULAR ELEMENTS **********************
            ! 
        endif




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
        SUBROUTINE CIRBOD_2(NODN,NCNE,SI,ETA,   &
                JK01,N01,XJ,XXJ,fitg)
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
            REAL*8 XJ(2,3),XXJ(3,3),xxx(3,8)
            REAL*8 det1,det2,det3,dumx,sf(8),dsf(2,8),ddsf(3,8),jk01(3,3),n01(3),param(2)
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

                            call comp_coef(csst,snst,xj,xxj,jk01,n01,f1,f2,param)

                        !fixme *pi4
                        !FITG2=FITG2-F2*(param(1)/param(2)**2+1.0/PLO)*WIQ1(I)
                        !FITG1=FITG1+F1*DLOG(PLO/param(2))*WIQ1(I)
                        FITG2=FITG2-F2*(param(1)/param(2)**2+1.0/PLO)*WIQ1(I)*pi4
                        FITG1=FITG1+F1*DLOG(PLO/param(2))*WIQ1(I)*pi4 

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


                        

                            call comp_coef(csst,snst,xj,xxj,jk01,n01,f1,f2,param)

                            fitg2=fitg2-f2*(param(1)/param(2)**2+1.0/plo)*wiq2(i)*pi4
                            fitg1=fitg1+f1*dlog(plo/param(2))*wiq2(i)*pi4



                            SUM1=SUM1+F1*WIQ2(I)                      !(25)
                            SUM2=SUM2+F2/BETAs*WIQ2(I)                !(26)
                            
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





         subroutine comp_coef(csst,snst,xj,xxj,jk0,n0,f1,f2,param) 
            
            use proj_cnst,only :pi4

            implicit none
            
            real(8),intent(in) :: csst,snst,xj(2,3),xxj(3,3),jk0(3,3),n0(3)
            real(8),intent(out) :: f1,f2,param(2)

            real(8) :: jk1(3),n1,scos(2)
            real(8) :: a(3),b(3),ast,bst,as_3,as_2,g31,b30,a30,b31,a31

            scos=[csst,snst]

            !jk1(1:3) = jk0(2,1:3)*csst+jk0(3,1:3)*snst
            jk1(1:3) = matmul(scos,jk0(2:3,:))

            !n1=n0(2)*csst+n0(3)*snst
            n1=dot_product(n0(2:3),scos)

            a(1:3) = xj(1,:)*csst+xj(2,:)*snst
            b(1:3) = xxj(1,:)*csst**2*0.5+xxj(2,:)*snst**2*0.5+xxj(3,:)*csst*snst

            ast=norm2(a)
            bst=norm2(b)

            as_3=1.0d0/ast**3
            as_2=-3.0d0*(dot_product(a,b))/ast**5
            g31=dot_product(b,jk0(1,:))+dot_product(a,jk1)
            g31=(a(3)/ast**2)*g31

            b30=-jk0(1,3)
            a30=b30*n0(1)

            b31=3*g31-jk1(3)
            a31=b31*n0(1)+b30*n1

            !!todo neg removed
            !f1=-(as_2*a30+as_3*a31)/pi4
            !f2=-(as_3*a30)/pi4
            f1=(as_2*a30+as_3*a31)/pi4
            f2=(as_3*a30)/pi4

            param(1) = -dot_product(a,b)/ast**4
            param(2) = 1.0d0/ast

        end subroutine
