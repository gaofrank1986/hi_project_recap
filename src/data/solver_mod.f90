module solver_mod
    use kinds
    type solver
        ! @var : [amata] lhs matrix in Ax = B
        ! @var : [bmata] rhs matrix in Ax = B
        ! @var : [cmata] used for generating bmata
        real(rk),allocatable                        :: A(:,:,:),B(:,:)
        real(rk),allocatable                        :: C(:,:,:)
        ! @var : [indx] store pivot info in LU decomposision
        integer, allocatable,private                :: piv(:,:)
        integer, private                            :: nsys,n,m
    contains
        procedure :: init=>init_solver
        procedure :: LUDecomp
        procedure :: Solve=>getUnknown
        !generic :: Sover => g1,g2
    end type
    !
    private :: RLUDCMP,rlubksb
contains
   subroutine init_solver(this,nsys,nnode,nnrml)
       implicit none
       class(solver) :: this
       integer,intent(in) :: nsys,nnode,nnrml
       real(rk) :: dsign

       this%nsys=nsys
       this%n=nnode
       this%m=nnrml

       if (allocated(this%A)) deallocate(this%A)
       if (allocated(this%B)) deallocate(this%B)
       if (allocated(this%C)) deallocate(this%C)
       if (allocated(this%piv)) deallocate(this%piv)

       allocate(this%A(nnode,nnode,nsys),&
           &    this%B(nnode,nsys), this%piv(nnode,nsys))
       allocate(this%C(nnode,nnrml,nsys))
       this%A=0.0d0
       this%B=0.0d0
       this%C = 0.0d0

   end subroutine

   subroutine LUDecomp(this)
       implicit none
       class(solver) :: this
       integer ip
       real(8) :: dsign

        do ip=1, this%nsys
            call rludcmp(ip,this%A,this%n,this%n,this%nsys,this%piv,dsign)  
        enddo

   end subroutine

   function getUnknown(this,bc) result(ans)
       use proj_cnst,only:rsn
       implicit none
       class(solver) :: this
       real(rk),dimension(this%m,this%nsys) :: ans
       real(8) bmat(this%nsys),bc(this%n,this%nsys)
       integer :: is,ind,ip

!        do is=1, this%nsys
            !do ind=1, this%n
                !!---------------loop-body--------------------------
                !this%B(ind,is)=dot_product(this%C(ind,:,is),bc(:,is))
                !! .. cmat d\phi/dp   ... phi
                !! .. cmata
                !!if (ind<=nnf) then!potential only
                !!bmata(ind,is) = bmata(ind,is)-fterm(ind,is,1)*cmat(ind,is)&
                !!&-fterm(ind,is,2)*dpoxyz_save(1,is,ind) &
                !!&-fterm(ind,is,3)*dpoxyz_save(2,is,ind)
                !!end if
                !!---------------------------------------------------
            !enddo
        !enddo 

        !Solving Ax=B      // B is time varying
        ! result saved in bmata
        do  is=1, this%nsys   
            call rlubksb(is,this%A,this%n,this%n, 1,this%nsys, 1,this%piv,this%B)
        enddo

        do ip=1, this%nsys
            do ind=1, this%n

                bmat(ip)=(0.0d0, 0.0d0)
                do  is=1, this%nsys
                    bmat(ip)=bmat(ip)+this%B(ind,is)*rsn(ip,is)
                enddo
                bmat(ip)=bmat(ip)/this%nsys
                ans(ind,ip)=bmat(ip)        
            enddo
        enddo

    end function
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
        PARAMETER (NMAX=20000, TINY=1.0E-20) 
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

end module
