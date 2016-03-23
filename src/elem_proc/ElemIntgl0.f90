      
!  NORM_ELE0+TSING0
!  NORM_INT0+SING_INT0
! ======================================================
!
!   Integration on an element without source in itself
!   nor its symmetrical ones
!
! ======================================================
        include './add_on/ElemIntgl1.f90'

    subroutine norm_ele0(ielem,xp,yp,zp,aval,bval)

        use mvar_mod
        use green_funcs,only:gcombo0
        implicit   none 

        integer  is,nd,j,np,ielem,ncne
        real*8  xp,yp,zp 
        real*8  aval(4,8),bval(4,8)

        ncne=ncn(ielem)

        aval=0.0d0
        bval=0.0d0

        do    is=1,   nsys  
            call norm_int(is,ielem,ncne,xp,yp,zp,aval(is,:),bval(is,:),gcombo0)
        end do

        end           

!
! ======================================================
!   Integration on an element with source in itself
!   or its mirror ones about any symmetrical axis
!
! ======================================================
!
       SUBROUTINE SING_ELE0(INODE,IELEM,NUMQUA,XP,YP,ZP,AVAL,BVAL)
	   USE MVAR_MOD
	   USE MFUNC_mod
       !use tripole_mod
        IMPLICIT   NONE  
!
	  INTEGER I,J,IS,IELEM,INODE,NODNUM,ND,NP,NUMQUA
        REAL*8  XP,YP,ZP,XYZT(3,8),DXYZT(3,8)
        REAL*8 AVAL(4,8),BVAL(4,8)
!
          AVAL= 0.0d0 
          BVAL= 0.0d0 
!

        DO 5     I=1,  NCN(IELEM)
        XYZT(1, I)  =  XYZE(1, I, IELEM)  
        XYZT(2, I)  =  XYZE(2, I, IELEM)  
        XYZT(3, I)  =  XYZE(3, I, IELEM)
	    
        DXYZT(1, I) = DXYZE(1, I, IELEM)  
        DXYZT(2, I) = DXYZE(2, I, IELEM)  
        DXYZT(3, I) = DXYZE(3, I, IELEM)

        IF(INODE.EQ.NCON(IELEM,I)) NODNUM=I
5       CONTINUE
!
        CALL TRIPOL(NODNUM,NCN(IELEM),XYZT,DXYZT)
!	  PRINT *,' AFTER TRIPOL'
!
!
        !write(*,*) 'NUMQUA=',NUMQUA,'nsys=',nsys
        !write(*,*) 'NODNUM=',nodnum
        !pause
        IF(NUMQUA.EQ.0)       THEN
         DO 100 IS=1,  NSYS
          IF(IS.EQ.1) THEN 
            CALL SING_INT0(IS,IELEM,nodnum,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.NE.1 ) THEN   
            CALL NORM_INT(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
100      CONTINUE
!
        ELSE IF(NUMQUA.EQ.2) THEN
         DO 200 IS=1,NSYS     
          IF(IS.EQ.1.OR.IS.EQ.2) THEN
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.EQ.3.OR.IS.EQ.4) THEN  
            !CALL NORM_INT(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
!
200      CONTINUE  
!
        ELSE IF(NUMQUA.EQ.4) THEN
         DO 300  IS=1,  NSYS
          IF(IS.EQ.1.OR.IS.EQ.4) THEN   
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.EQ.2.OR.IS.EQ.3) THEN  
            !CALL NORM_INT(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
300      CONTINUE
!
        ELSE IF(NUMQUA.EQ.5) THEN
         DO 400 IS=1, NSYS  
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
400      CONTINUE
        ENDIF
!
        RETURN
        END
        
    !<  -----------------------------------------------
    !<  element integration with src point on one node
    !!  method only capable of Gfunc desingularizaito
    
    subroutine sing_int0(is,ielem,nodnum,xp,yp,zp,aval,bval)
        use mvar_mod
        use green_funcs,only:gcombo0
        use proj_cnst,only:ex,ey
        use trvar_mod
        implicit   none  
        
        integer is,ielem,n,j,ip,nodnum
        real*8  xp,yp,zp!,ex(4,4),ey(4,4)
        real*8  x,y,z,x0,y0,z0      
        real*8  nx,ny,nz,dgn
        real*8  aval(4,8),bval(4,8),gxf(4),p(3),p0(3),np(3)


        x0=ex(is)*xp
        y0=ey(is)*yp
        z0= zp
        p0=(/x0,y0,z0/)

        do n=1, nosamp

            p= xynod(1:3,n)
            np=dxynod(1:3,n)

            call gcombo0 (h,p,p0,gxf) 
            dgn = dot_product(gxf(2:4),np)

            do j=1, ncn(ielem)
                aval(is,j)=aval(is,j)+dgn*samnod(n,j)
                bval(is,j)=bval(is,j)+gxf(1)*samnod(n,j)
            enddo
        enddo
            !if(ielem>nelemf) then
            !write(9011,'(2i6,8f12.6)') ielem,nodnum,aval(is,:)
            !call sing_int1(is,ielem,nodnum,xp,yp,zp,aval,bval) 
            !end if
            !
    end subroutine



