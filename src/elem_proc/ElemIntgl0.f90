      
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
        SUBROUTINE NORM_ELE0(IELEM,XP,YP,ZP,AVAL,BVAL)
        
	    USE MVAR_MOD
        IMPLICIT   NONE 
	  
	    INTEGER  IS,ND,J,NP,IELEM,NCNE
	    REAL*8  XP,YP,ZP 
        REAL*8  AVAL(4,8),BVAL(4,8)
!   
!            
	      NCNE=NCN(IELEM)
	                 
!	 Print *,' In TINBOD   PI=',PI             
!
          AVAL=0.0D0
          BVAL=0.0D0
!
        DO 100   IS=1,   NSYS  
          CALL NORM_INT0(IS,IELEM,NCNE,XP,YP,ZP,AVAL,BVAL)
100     CONTINUE
!
        RETURN
        END           

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
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.NE.1 ) THEN   
            CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
100      CONTINUE
!
        ELSE IF(NUMQUA.EQ.2) THEN
         DO 200 IS=1,NSYS     
          IF(IS.EQ.1.OR.IS.EQ.2) THEN
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.EQ.3.OR.IS.EQ.4) THEN  
            CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
          END IF
!
200      CONTINUE  
!
        ELSE IF(NUMQUA.EQ.4) THEN
         DO 300  IS=1,  NSYS
          IF(IS.EQ.1.OR.IS.EQ.4) THEN   
            CALL SING_INT0(IS,IELEM,XP,YP,ZP,AVAL,BVAL)
          ELSE IF(IS.EQ.2.OR.IS.EQ.3) THEN  
            CALL NORM_INT0(IS,IELEM,NCN(IELEM),XP,YP,ZP,AVAL,BVAL)
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
                           
!C ======================================================
!C
!C Integration on an element without source point
!C 
!C ======================================================
!C                      


      subroutine norm_int0(is,ielem,ncne,xp,yp,zp,aval,bval)


      use mvar_mod
      use green_funcs,only : GFunc,DGFunc,mirror
      use proj_cnst,only:ex,ey
      implicit   none  

      integer is,ielem,n,nsamb,ncne,j

      real*8  xp,yp,zp
      real(8) :: p(3),p0(3),np(3)
      !integer,optional,intent(in) :: debug

      real*8  v1,v2
      real*8  aval(4,8),bval(4,8)
      
      aval=0.0d0
      bval=0.0d0
      
      p0 = (/ex(is)*xp,ey(is)*yp,zp/)
      !if (present(debug).and.(debug)) print *,x0,y0,z0
      nsamb=16
      if(ncne.eq.6)   nsamb=4

      do 100    n=1,   nsamb     

      p =sambxy(ielem,n,1:3)
      np = dsamb(ielem,n,1:3)
      v1 = GFunc(p,p0)+GFunc(p,mirror(h,p0))
      v2 = dot_product(np,DGFunc(p,p0)+DGFunc(p,mirror(h,p0)))
      do  j=1,   ncne
      bval(is,j)=bval(is,j)+v1*samb(ielem,n,j)
      aval(is,j)=aval(is,j)+v2*samb(ielem,n,j)
      enddo
      100     continue

      end subroutine

! ===========================================================
! Integration on an element with source point
!
      subroutine sing_int0(is,ielem,xp,yp,zp,aval,bval)

            
        use trvar_mod
        use mvar_mod
        use green_funcs,only : GFunc,DGFunc,mirror
        use proj_cnst,only:ex,ey
        implicit   none  

      real(8) :: p(3),p0(3),np(3)
      integer is,ielem,n,j,ip       
      real*8  xp,yp,zp!,ex(4,4),ey(4,4)
      real*8  v1,v2
      real*8  aval(4,8),bval(4,8)

      p0 = (/ex(is)*xp,ey(is)*yp,zp/)

      do 130 n=1, nosamp

      p =xynod(1:3,n)
      np = dxynod(1:3,n)
      
      v1 = gfunc(p,p0)+gfunc(p,mirror(h,p0))
      v2 = dot_product(np,dgfunc(p,p0)+dgfunc(p,mirror(h,p0)))

      do j=1, ncn(ielem)
      aval(is,j)=aval(is,j)+v2*samnod(n,j)
      bval(is,j)=bval(is,j)+v1*samnod(n,j)
      enddo

      130     continue
      return
      end

