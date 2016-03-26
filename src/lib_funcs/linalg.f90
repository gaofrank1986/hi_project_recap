module linalg
contains
! ********************************************* 
! * SOLUTION OF LINEAR EQUATIONS  [A]{X}= {B} * 
! *           LU DECOMPOSITION                * 
! * FROM 'NUMERICAL RECIPES'     pp. 35-37    * 
! ********************************************* 
!
!
! ---------------------------------------- 
! NSYS: ¾ØÕóA[:,:]µÄ¸öÊý£¬ÓÃÓÚ¿ªÊý×é 
! IP  : ±¾´Î¼ÆËãµÄ¾ØÕóÐòºÅ
! A   :
! N   :
! NP  :
! LI  :
! NMOD:
! INDX:
! B   :
! ---------------------------------------- 
! 

        SUBROUTINE RLUDCMP(IP,A,N,NP,NSYS,INDX,D)           
        IMPLICIT REAL*8(A-H,O-Z)  
        PARAMETER (NMAX=5000, TINY=1.0E-20) 
        DIMENSION INDX(NP,NSYS) 
        REAL*8 SUM,DUM,AAMAX 
        REAL*8 A(NP,NP,NSYS),VV(NMAX) 

        D=1. 
        DO 12 I=1, N 
        AAMAX=(0., 0.) 
          DO 11 J=1, N 
          IF ( DABS(A(I,J,IP)).GT. DABS(AAMAX) )  &
                             AAMAX=A(I,J,IP) 
11      CONTINUE 
!
        IF (DABS(AAMAX) .EQ. 0.0)  THEN	  
	    Print  *, ' SINGULAR MATRIX   inside RLUDCMP' 
	    Print  *, ' IP=',IP,' I=',I
		PAUSE
        ENDIF
!
	  VV(I)=1./AAMAX 
12      CONTINUE 
        DO 19 J=1, N 
          DO 14 I=1, J-1 
            SUM=A(I,J,IP) 
            DO 13 K=1, I-1 
              SUM=SUM-A(I,K,IP)*A(K,J,IP) 
13            CONTINUE 
            A(I,J,IP)=SUM 
14          CONTINUE 
          AAMAX=(0., 0.) 
          DO 16 I=J, N                          
            SUM=A(I,J,IP) 
            DO 15 K=1, J-1 
              SUM=SUM-A(I,K,IP)*A(K,J,IP) 
15            CONTINUE 
            A(I,J,IP)=SUM 
            DUM=VV(I)*SUM 
            IF (DABS(DUM) .GE. DABS(AAMAX)) THEN 
            IMAX=I 
            AAMAX=DUM 
            END IF 
16          CONTINUE 
         IF(J .NE. IMAX) THEN 
         DO 17 K=1, N 
           DUM=A(IMAX,K,IP) 
           A(IMAX,K,IP)=A(J,K,IP) 
           A(J,K,IP)=DUM 
17         CONTINUE 
          D=-D 
          VV(IMAX)=VV(J) 
          END IF 
          INDX(J,IP)=IMAX 
          IF(A(J,J,IP).EQ.0.) A(J,J,IP)=TINY 
          IF(J.NE.N) THEN 
          DUM=1./A(J,J,IP) 
          DO 18 I=J+1,N 
            A(I,J,IP)=A(I,J,IP)*DUM 
18          CONTINUE 
          END IF 
19        CONTINUE 
        RETURN 
	 END SUBROUTINE RLUDCMP





!    
! ---------------------------------------- 
! A   :
! N   :
! NP  :
! IP  :
! NSYS:
!
! B   :
! LI  :
! NMOD:
!
! INDX:
! ----------------------------------------
!
        SUBROUTINE RLUBKSB(IP,A,N,NP,LI,NSYS,NMOD,INDX,B)  
        IMPLICIT REAL*8(A-H,O-Z)  
        DIMENSION INDX(NP,NSYS) 
        REAL*8 SUM,A(NP,NP,NSYS),B(NP,NMOD,NSYS) 

        II=0 
        DO 12 I=1,N 
         LL=INDX(I,IP) 
        SUM=B(LL,LI,IP) 
        B(LL,LI,IP)=B(I,LI,IP) 
        IF(II.NE.0) THEN 
        DO 11 J=II, I-1 
        SUM=SUM-A(I,J,IP)*B(J,LI,IP) 
11      CONTINUE 
        ELSE IF (DABS(SUM).NE. 0.) THEN 
        II=I 
        END IF 
        B(I,LI,IP)=SUM 
12      CONTINUE 
        DO 14 I=N, 1, -1 
        SUM=B(I,LI,IP) 
        DO 13 J=I+1, N 
        SUM=SUM-A(I,J,IP)*B(J,LI,IP) 
13      CONTINUE 
        B(I,LI,IP)=SUM/A(I,I,IP) 
14      CONTINUE 
        RETURN 
	  END SUBROUTINE RLUBKSB


    !C 
    !C ********************************************* 
    !C * Inverse of a matric  [A] by               * 
    !C *           LU decomposition                * 
    !C * From 'NUMERICAL RECIPES'     pp. 35-37    * 
    !C ********************************************* 
    !C 
       SUBROUTINE INVERSE(a,y,n,np)
       INTEGER n,np,indx(n)
       REAL*8 a(np,np),b(np),y(np,np),d

	 CALL ludcmp(a,n,np,indx,d)
	  y=0.0d0
       do 10 i=1, n
	  y(i,i)=1.0d0
10	 continue

	 do 100 j=1, n
	  do 20 i=1, n
	   b(i)=y(j,i) !i,j
20	  continue

        call  lubksb(a,n,np,indx,b)

	  do 30 i=1, n
	   y(j,i)=b(i) !j,i
30	  continue

100	 continue


      END	 SUBROUTINE INVERSE




    !C ********************************************* 
    !C * SOLUTION OF LINEAR EQUATIONS  [A]{X}= {B} * 
    !C *           LU DECOMPOSITION                * 
    !C * FROM 'NUMERICAL RECIPES'     pp. 35-37    * 
    !C ********************************************* 

      SUBROUTINE ludcmp(a,n,np,indx,d)

      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (DAbs(a(i,j)).gt.aamax) aamax=DAbs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*DAbs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax) then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n) then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END SUBROUTINE ludcmp




      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12     continue
!       
	 do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14     continue
       return
       END	 SUBROUTINE lubksb



!C ***********************************************************
!C         From Chapter  2.6 of Numerical Recipes
!C       
!C
!C ---------------------------------------------------------
!C   Given a matrix A, with logical dimensions M by N and physical 
!C dimensions MP by NP, this routine computes its singular value 
!C decomposition, A=U.W.V^T. The matrix U replaces A on output. 
!C The diagonal matrix of singular values W is output as a vector W. 
!C The matrix V (not the transpose V^T) is output as V. M must be 
!C greater or equal to N; if it is amaller, then A should be filled 
!C up to square with zero rows. 
!C **************************************************************
!      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)  
!      IMPLICIT  REAL*4(A-H, O-Z)
!      REAL*8    a(mp,np),v(np,np),w(np)
!      PARAMETER (NMAX=5000)
!C
!      DIMENSION rv1(NMAX)
!C
!      g=0.0
!      scale=0.0
!      anorm=0.0
!      do 25 i=1,n
!        l=i+1
!        rv1(i)=scale*g
!        g=0.0
!        s=0.0
!        scale=0.0
!        if(i.le.m)then
!          do 11 k=i,m
!            scale=scale+Abs(a(k,i))
!11        continue
!          if(scale.ne.0.0)then
!            do 12 k=i,m
!              a(k,i)=a(k,i)/scale
!              s=s+a(k,i)*a(k,i)
!12          continue
!            f=a(i,i)
!            g=-sign(Sqrt(s),f)
!            h=f*g-s
!            a(i,i)=f-g
!            do 15 j=l,n
!              s=0.0
!              do 13 k=i,m
!                s=s+a(k,i)*a(k,j)
!13            continue
!              f=s/h
!              do 14 k=i,m
!                a(k,j)=a(k,j)+f*a(k,i)
!14            continue
!15          continue
!            do 16 k=i,m
!              a(k,i)=scale*a(k,i)
!16          continue
!          endif
!        endif
!        w(i)=scale *g
!        g=0.0
!        s=0.0
!        scale=0.0
!        if((i.le.m).and.(i.ne.n))then
!          do 17 k=l,n
!            scale=scale+Abs(a(i,k))
!17        continue
!          if(scale.ne.0.0)then
!            do 18 k=l,n
!              a(i,k)=a(i,k)/scale
!              s=s+a(i,k)*a(i,k)
!18          continue
!            f=a(i,l)
!            g=-sign(Sqrt(s),f)
!            h=f*g-s
!            a(i,l)=f-g
!            do 19 k=l,n
!              rv1(k)=a(i,k)/h
!19          continue
!            do 23 j=l,m
!              s=0.0
!              do 21 k=l,n
!                s=s+a(j,k)*a(i,k)
!21            continue
!              do 22 k=l,n
!                a(j,k)=a(j,k)+s*rv1(k)
!22            continue
!23          continue
!            do 24 k=l,n
!              a(i,k)=scale*a(i,k)
!24          continue
!          endif
!        endif
!        anorm=max(anorm,(Abs(w(i))+Abs(rv1(i))))
!25    continue
!      do 32 i=n,1,-1
!        if(i.lt.n)then
!          if(g.ne.0.0)then
!            do 26 j=l,n
!              v(j,i)=(a(i,j)/a(i,l))/g
!26          continue
!            do 29 j=l,n
!              s=0.0
!              do 27 k=l,n
!                s=s+a(i,k)*v(k,j)
!27            continue
!              do 28 k=l,n
!                v(k,j)=v(k,j)+s*v(k,i)
!28            continue
!29          continue
!          endif
!          do 31 j=l,n
!            v(i,j)=0.0
!            v(j,i)=0.0
!31        continue
!        endif
!        v(i,i)=1.0
!        g=rv1(i)
!        l=i
!32    continue
!      do 39 i=min(m,n),1,-1
!        l=i+1
!        g=w(i)
!        do 33 j=l,n
!          a(i,j)=0.0
!33      continue
!        if(g.ne.0.0)then
!          g=1.0/g
!          do 36 j=l,n
!            s=0.0
!            do 34 k=l,m
!              s=s+a(k,i)*a(k,j)
!34          continue
!            f=(s/a(i,i))*g
!            do 35 k=i,m
!              a(k,j)=a(k,j)+f*a(k,i)
!35          continue
!36        continue
!          do 37 j=i,m
!            a(j,i)=a(j,i)*g
!37        continue
!        else
!          do 38 j= i,m
!            a(j,i)=0.0
!38        continue
!        endif
!        a(i,i)=a(i,i)+1.0
!39    continue
!      do 49 k=n,1,-1
!        do 48 its=1,30
!          do 41 l=k,1,-1
!            nm=l-1
!            if((Abs(rv1(l))+anorm).eq.anorm)  goto 2
!            if((Abs(w(nm))+anorm).eq.anorm)  goto 1
!41        continue
!1         c=0.0
!          s=1.0
!          do 43 i=l,k
!            f=s*rv1(i)
!            rv1(i)=c*rv1(i)
!            if((Abs(f)+anorm).eq.anorm) goto 2
!            g=w(i)
!            h=Sqrt(f*f+g*g)
!            w(i)=h
!            h=1.0/h
!            c= (g*h)
!            s=-(f*h)
!            do 42 j=1,m
!              y=a(j,nm)
!              z=a(j,i)
!              a(j,nm)=(y*c)+(z*s)
!              a(j,i)=-(y*s)+(z*c)
!42          continue
!43        continue
!2         z=w(k)
!          if(l.eq.k)then
!            if(z.lt.0.0)then
!              w(k)=-z
!              do 44 j=1,n
!                v(j,k)=-v(j,k)
!44            continue
!            endif
!            goto 3
!          endif
!          if(its.eq.30) pause 'no convergence in svdcmp'
!          x=w(l)
!          nm=k-1
!          y=w(nm)
!          g=rv1(nm)
!          h=rv1(k)
!          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
!          g=Sqrt(f*f+1.0)
!          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
!          c=1.0
!          s=1.0
!          do 47 j=l,nm
!            i=j+1
!            g=rv1(i)
!            y=w(i)
!            h=s*g
!            g=c*g
!            z=Sqrt(f*f+h*h)
!            rv1(j)=z
!            c=f/z
!            s=h/z
!            f= (x*c)+(g*s)
!            g=-(x*s)+(g*c)
!            h=y*s
!            y=y*c
!            do 45 jj=1,n
!              x=v(jj,j)
!              z=v(jj,i)
!              v(jj,j)= (x*c)+(z*s)
!              v(jj,i)=-(x*s)+(z*c)
!45          continue
!            z=Sqrt(f*f+h*h)
!            w(j)=z
!            if(z.ne.0.0)then
!              z=1.0/z
!              c=f*z
!              s=h*z
!            endif
!            f= (c*g)+(s*y)
!            x=-(s*g)+(c*y)
!            do 46 jj=1,m
!              y=a(jj,j)
!              z=a(jj,i)
!              a(jj,j)= (y*c)+(z*s)
!              a(jj,i)=-(y*s)+(z*c)
!46          continue
!47        continue
!          rv1(l)=0.0
!          rv1(k)=f
!          w(k)=x
!48      continue
!3       continue
!49    continue
!      return
!      END	SUBROUTINE svdcmp
!C
!C -----------------------------------------------------------------
!C  Solve A.X=B for a vector X, where A is specified by the arrays 
!C  U, W, V as returned by SVDCMP. M and N are the logical dimensions 
!C  of A, and will be equal for square matrices. MP and NP are the 
!C  physical dimensions of A. B is the input right-hand side. X is 
!C  the output solution vector. No input quantities are destroyed, 
!C  so the routine may be called sequentially with different B's. '
!C  M must be greater or equal to N; see SVDCMP
!C
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      IMPLICIT  REAL*8(A-H, O-Z)

      REAL*8    b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=5000)
      DIMENSION tmp(NMAX) 
      
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END SUBROUTINE svbksb





end module

