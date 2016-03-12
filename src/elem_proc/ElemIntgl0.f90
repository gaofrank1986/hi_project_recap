      
!  NORM_ELE0+TSING0
!  NORM_INT0+SING_INT0
! ======================================================
!
!   Integration on an element without source in itself
!   nor its symmetrical ones
!
! ======================================================
        include './add_on/ElemIntgl1.f90'
!
      subroutine norm_ele0(ielem,xp,yp,zp,aval,bval)

      use mvar_mod
      implicit   none 

      integer  is,nd,j,np,ielem,ncne
      real*8  xp,yp,zp 
      real*8  aval(4,8),bval(4,8)
      real*8  atmp(4,8),btmp(4,8)
      !   
      !            
      ncne=ncn(ielem)

      !        print *,' in tinbod   pi=',pi             
      !
      aval=0.0d0
      bval=0.0d0
      atmp=0
      btmp=0
      !
      do 100   is=1,   nsys  
      call norm_int0(is,ielem,ncne,xp,yp,zp,aval,bval,1)
      !call norm_int0(is,ielem,ncne,xp,yp,-(2*h+zp),atmp,btmp,1)
      aval=aval+atmp
      bval=bval+btmp
      100     continue
      !
      return
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
        REAL*8 ATMP(4,8),BTMP(4,8)
!
          AVAL= 0.0d0 
          ATMP=0
          BTMP=0
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
!         PRINT *,' AFTER TRIPOL'
!
        !write(*,*) 'NUMQUA=',NUMQUA,'nsys=',nsys
        !write(*,*) 'NODNUM=',nodnum
        !pause
        !IF(NUMQUA.EQ.0)       THEN
         DO 100 IS=1,  NSYS
          IF(IS.EQ.1) THEN 
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
        !if (inode.eq.942) then
            !CALL sing_INT0(IS,IELEM,XP,YP,-(2*H+ZP),ATMP,BTMP)!,2,1)
        !else
            !CALL sing_INT0(IS,IELEM,XP,YP,-(2*H+ZP),ATMP,BTMP)!,2,0)
        !end if

          aval=aval+atmp
          bval=bval+btmp
          ELSE IF(IS.NE.1 ) THEN   
            !CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
100      CONTINUE
!
!        ELSE IF(NUMQUA.EQ.2) THEN
         !DO 200 IS=1,NSYS     
          !IF(IS.EQ.1.OR.IS.EQ.2) THEN
            !CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          !ELSE IF(IS.EQ.3.OR.IS.EQ.4) THEN  
            !CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          !END IF
!!
!200      CONTINUE  
!!
        !ELSE IF(NUMQUA.EQ.4) THEN
         !DO 300  IS=1,  NSYS
          !IF(IS.EQ.1.OR.IS.EQ.4) THEN   
            !CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          !ELSE IF(IS.EQ.2.OR.IS.EQ.3) THEN  
            !CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          !END IF
!300      CONTINUE
!!
        !ELSE IF(NUMQUA.EQ.5) THEN
         !DO 400 IS=1, NSYS  
            !CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
!400      CONTINUE
        !ENDIF
!
        RETURN
        END
                           
!C ======================================================
!C
!C Integration on an element without source point
!C 
!C ======================================================
!C                      
      subroutine norm_int0(is,ielem,ncne,xp,yp,zp,aval,bval,kind,debug)
      use mvar_mod
      use green_funcs,only : GFunc,DGFunc,mirror
      use proj_cnst,only:ex,ey
      !use mfunc_mod
      implicit   none  

      integer is,ielem,n,nsamb,ncne,j,kind

      real*8  xp,yp,zp
      real(8) :: p(3),p0(3),np(3)
      integer,optional,intent(in) :: debug

      real*8  v1,v2
      real*8  aval(4,8),bval(4,8)
      
      p0 = (/ex(is)*xp,ey(is)*yp,zp/)
      !if (present(debug).and.(debug)) print *,x0,y0,z0
      nsamb=16
      if(ncne.eq.6)   nsamb=4

      do 100    n=1,   nsamb     

      p =sambxy(ielem,n,1:3)
      np = dsamb(ielem,n,1:3)
      v1 = GFunc(p,p0)+GFunc(p,mirror(h,p0))
      v2 = dot_product(np,DGFunc(p,p0)+DGFunc(p,mirror(h,p0)))
      if (kind .eq.1) then
            do  j=1,   ncne
            bval(is,j)=bval(is,j)+v1*samb(ielem,n,j)
            aval(is,j)=aval(is,j)+v2*samb(ielem,n,j)
            enddo
      else

            do  j=1,   ncne
            bval(is,j)=bval(is,j)+v2*samb(ielem,n,j)
            aval(is,j)=aval(is,j)+v1*samb(ielem,n,j)
            enddo
      end if
      100     continue
      end
!
! ===========================================================
! Integration on an element with source point
!
       SUBROUTINE SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
           USE MVAR_MOD
            USE MFUNC_mod
           USE TRVAR_MOD
       IMPLICIT   NONE  
!
           INTEGER IS,IELEM,N,J,IP       
           REAL*8  XP,YP,ZP,EX(4,4),EY(4,4)
           REAL*8  X,Y,Z,X0,Y0,Z0      
           REAL*8  NX,NY,NZ,DGN
       REAL*8  AVAL(4,8),BVAL(4,8),GXF(4)
!
          DATA EX/  1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
                1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0/
!                                                  
          DATA EY/  1.0d0, -1.0d0, -1.0d0,  1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
                1.0d0, -1.0d0, -1.0d0,  1.0d0/
!
!   
201       FORMAT(' In SING_INT0    XP, YP, ZP=',3F12.4)
!
 
          X0=EX(IS,1)*XP
          Y0=EY(IS,1)*YP
          Z0= ZP

      DO 130 N=1, NOSAMP
        
          X =XYNOD(1,N)
          Y =XYNOD(2,N)
          Z =XYNOD(3,N)

          !CALL TGRN (H,X,X0,Y,Y0,Z,Z0,GXF) 
          CALL  TGRN2(H,X,X0,Y,Y0,Z,Z0,GXF) 
!
        DGN=GXF(2)*DXYNOD(1,N)+GXF(3)*DXYNOD(2,N)+   &
              GXF(4)*DXYNOD(3,N)             

!         write(6,*) ' DGN=',DGN

        DO J=1, NCN(IELEM)
         AVAL(IS,J)=AVAL(IS,J)+DGN*SAMNOD(N,J)
         BVAL(IS,J)=BVAL(IS,J)+GXF(1)*SAMNOD(N,J)
       ENDDO
!

130     CONTINUE
!
        RETURN
        END


