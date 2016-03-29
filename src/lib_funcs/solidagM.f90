
    ! @params   : inode
    ! @params   : nnode
    ! @params   : nelem
    ! @params   : ncn
    ! @params   : ncon
    ! @params   : nodqua,h,xyz,dxyze
    ! @return   : sangle
    subroutine solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
            h,xyz,dxyze,sangle)

        use shape_funcs

        implicit   none 

        integer,intent(in):: inode,nnode,nelem
        integer,intent(in):: ncn(nelem),ncon(nelem,8),nodqua(nnode)

        real*8,intent(in)::  h,xyz(3,nnode),dxyze(3,8,nelem) 
        real*8,intent(out):: sangle

        integer  ielem,i,j,lk,lj,li,iel,minside
        integer  mcele,mcnt(100),mend(100),mnod(100)
        integer  ncont(100),nside(100,2)
        integer  nods_even(4),nods_odd(8),nodt_even(3),nodt_odd(6)

        real*8   pi
        real*8   dr,si,eta,dum
        real*8   sf(8),dsf(2,8),siq(12),etaq(12),dsiq(12,2),xxx(3,8)

        real*8   dxyzn(100,3),txyzn(100,2,3),txyz(2,3),xj(2,3)
        !
        data pi/3.14159265359/

        data  siq/-1.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0,  &
            1.0d0, 0.0d0, 0.0d0,-1.0d0,-1.0d0,-1.0d0/ 

        data etaq/-1.0d0,-1.0d0,-1.0d0,-1.0d0, 0.0d0, 0.0d0,  &
            1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0/ 
        !
        data  dsiq/ 1.0d0,-1.0d0, 1.0d0,-1.0d0,-1.0d0,-1.0d0, &
            -1.0d0, 1.0d0,-1.0d0, 1.0d0, 1.0d0, 1.0d0, &    
            1.0d0, 1.0d0, 1.0d0, 1.0d0,-1.0d0, 1.0d0, &
            -1.0d0,-1.0d0,-1.0d0,-1.0d0, 1.0d0,-1.0d0/ 
        !

        data nods_even/1,4,7,10/
        data nods_odd/2,3,5,6,8,9,11,12/
        !
        data nodt_even/1,2,3/
        data nodt_odd/4,5,6,7,8,9/



        dxyzn=0.0d0
        txyzn=0.0d0

        if(inode ==1) then
            print *
            print *,'  computing solid angle by direct method'
            print *,'the program is developed by prof. bin teng at dut, china'
            print *     
        endif


        mcele=0
        do ielem=1, nelem;do i=1,ncn(ielem)

            if(inode .eq. ncon(ielem,i)) then
                if(mcele .gt. 25) then
                    print *,' there are ', mcele,' elements which have the node'
                    print *,' to inlarge the size of some arrays'
                    pause
                endif

                if(ncn(ielem) .eq. 8)  then

                    select case(i)
                    case(1,3,5,7)
                        mcele=mcele+1       ! num of links
                        mcnt(mcele)=ielem   !record connected elem id
                        mend(mcele)=i       ! record node pos
                        mnod(mcele)=nods_even((i+1)/2) ! new node pos id given
                    case(2,4,6,8)
                        mcele=mcele+1
                        mcnt(mcele)=ielem
                        mend(mcele)=i
                        mnod(mcele)=nods_odd(i-1)

                        mcele=mcele+1
                        mcnt(mcele)=ielem
                        mend(mcele)=i
                        mnod(mcele)=nods_odd(i)
                    end select
                elseif (ncn(ielem)==6) then
                    select case(i)
                    case(1,3,5)
                        mcele=mcele+1
                        mcnt(mcele)=ielem
                        mend(mcele)=i
                        mnod(mcele)=nodt_even(i)
                    case(2,4,6)
                        mcele=mcele+1
                        mcnt(mcele)=ielem
                        mend(mcele)=i
                        mnod(mcele)=nodt_odd(2*(i-3)-1)

                        mcele=mcele+1
                        mcnt(mcele)=ielem
                        mend(mcele)=i
                        mnod(mcele)=nodt_odd(2*(i-3))
                    end select
                endif
            endif
        enddo;enddo

        
        do iel=1, mcele

            ielem=mcnt(iel)
            ! mend(iel) is pos of inode at ielem
            dxyzn(iel,1:3) = dxyze(1:3,mend(iel),ielem)!

            IF(NCN(IELEM).EQ.8) THEN

                do li=1,8
                    xxx(1:3,li) = xyz(1:3,ncon(ielem,li))
                end do

                si =siq(mnod(iel))
                eta=etaq(mnod(iel))

                call spfunc8(si,eta,sf,dsf)

                xj(1:2,1:3)=matmul(dsf(1:2,1:8),transpose(xxx))

                do li=1,2
                    dr = norm2(xj(li,1:3))
                    txyzn(iel,li,1:3)  = dsiq(mnod(iel),li)*xj(li,1:3)/dr
                end do

            else if(ncn(ielem).eq.6) then
                call tangel6(ielem,mnod(iel),nnode,nelem,ncon,xyz,txyz)
                txyzn(iel,:,:)=txyz(:,:)
            endif
        enddo


        ! todo?? mirrored?
        if(abs(xyz(3,inode)+h) .lt. 1.0d-6)  then
            call xysym(dxyzn,txyzn,mcele)
        endif

        if(nodqua(inode) .eq. 2) then
            call zxsym(dxyzn,txyzn,mcele)
        else if(nodqua(inode) .eq. 4) then
            call yzsym(dxyzn,txyzn,mcele)
        else if(nodqua(inode) .eq. 5) then
            call yzzxsym(dxyzn,txyzn,mcele)
        endif

        !< order linked elements
        call domainfas1(dxyzn,txyzn,mcele,ncont,nside,minside)

        if(minside .gt.0) then
            call angle0(dxyzn,txyzn,mcele,ncont,nside,sangle)
            sangle=1.0d0-sangle/(4.0d0*pi)

        else
            sangle=0.5
        endif

    end subroutine


    !< Order linked element
    !! 
    !< @params dxyzn :  eleent normal vector
    !< @params txyzn :  element face vectors(2 for each)
    !< @return nside : order of face vecotr
    !< @return ncont : ordered list of element
    !< @return min?  :
    subroutine  domainfas1(dxyzn,txyzn,mcele,ncont,nside,minside)
        use linalg,only:cross_product
        implicit   none 

        integer  i,j,k,n,na,naxis,next_edge(1)
        integer  mcele,ncont(100),nside(100,2)

        integer  lct(100),lsd(100),minside

        real*8   dxyzn(100,3),txyzn(100,2,3)
        real*8   vec(3),dn(100),dmin,dot

        i=1  
        ncont(i)=1
        vec=cross_product(txyzn(i,1,:),txyzn(i,2,:))
        dot=dot_product(vec,dxyzn(1,1:3))
        ! decide if normal vector decided by txyzn and dxyzn are same
        if (dot .le. 0.0d0) then 
            !fixme why .le.
            nside(i,:) = [1,2]
        else
            nside(i,:) = [2,1]
        end if


        !print *,"mcele=",mcele
        do i=2,mcele
            na = 0

            do j=2,mcele

                if (any(ncont(1:i-1)==j) ) cycle

                do k=1,2

                    na=na+1
                    dn(na)=norm2(txyzn(j,k,:)-txyzn(ncont(i-1),nside(i-1,2),:))
                    ! distant to the last face vector
                    lct(na)=j! ielem
                    lsd(na)=k! kth face vector

                    !print *,"na=",na,dn(na),lct(na),lsd(na)
                end do
            end do
            next_edge = minloc(dn(1:na))
            !print *,next_edge(1)
            minside= next_edge(1) 
            ! check if min value is unrealistic
            ! todo if min too big (0.5) stop program
            ncont(i) = lct(next_edge(1)) 
            nside(i,1) = lsd(next_edge(1)) 
            nside(i,2) = 3 - nside(i,1)
            !print *,"selected",ncont(i),nside(i,1)
        end do
        !print *,"nside",nside(1:mcele,1)
        !print *,"nside",nside(1:mcele,2)
        !print *,"ncont",ncont(1:mcele)
        !print *,"min",minside

    end subroutine


    !< Evaluate angle given info
    !!
    ! @params dxyzn : normal vecotr
    ! @params txyzn : face vecotr
    ! @params mcele : num of connected edge
    ! @params ncont : ordered connected elem
    ! @params nside : edge order
    ! @return angle :

    subroutine  angle0(dxyzn,txyzn,mcele,ncont,nside,angle)
        use linalg,only:cross_product
        implicit none

        real(8),intent(out) :: angle

        integer  i,n,m
        integer  mcele,mcnt(100),mnod(100)
        integer  ncont(100),nside(100,2)
        real*8   dxyzn(100,3),txyzn(100,2,3)
        real*8   ang(100),vec(3),dot1,dot2,sign,pi


        angle=0.0
        pi=4.0d0*datan(1.0d0)

        !DO I=1, MCELE-1
        do i =1,mcele
            N=NCONT(I)
            if (i==mcele) then
                m=ncont(1)
            else
                m=ncont(i+1)
            endif
            !DOT1=DXYZN(N,1)*DXYZN(M,1)+DXYZN(N,2)*DXYZN(M,2)+ 

            !1	    DXYZN(N,3)*DXYZN(M,3)

            !< angle between two edge
            dot1=dot_product(dxyzn(n,1:3),dxyzn(m,1:3))
            vec=cross_product(dxyzn(n,1:3),dxyzn(m,1:3))
            dot2=dot_product(vec(1:3),txyzn(n,nside(i,2),1:3))
            if(abs(dot2) .lt. 1.0d-10) then
                sign=1.0d0
            else
                sign=dot2/abs(dot2)
            endif

            if(dot1  .gt. 1.0d0)   dot1= 1.0d0
            if(dot1  .lt.-1.0d0)   dot1=-1.0d0

            ang(i)=sign*acos(dot1)

        enddo

        angle=2.0d0*pi-sum(ang(1:mcele))

    end subroutine

    !
    subroutine symcore(dxyzn,txyzn,mcele,flag)
        implicit   none 

        integer  mcele,flag
        real*8   dxyzn(100,3),txyzn(100,2,3)

        real(8) :: pre(3)
        integer :: i,j
        select case(flag)
        case(1)!xy
            pre=[1.0d0,1.0d0,-1.0d0]
        case(2)!xz
            pre=[1.0d0,-1.0d0,1.0d0]
        case(3)!yz
            pre=[-1.0d0,1.0d0,1.0d0]
        end select

        do i=1, mcele
            j=i+mcele
            dxyzn(j,1:3) = pre*dxyzn(i,1:3)
            txyzn(j,1,1:3)= pre*txyzn(i,1,1:3)
            txyzn(j,2,1:3)= pre*txyzn(i,2,1:3)
        end do
        mcele=mcele+mcele

    end subroutine

    subroutine  xysym(dxyzn,txyzn,mcele)
        implicit   none 

        integer  mcele
        real*8   dxyzn(100,3),txyzn(100,2,3)

        call symcore(dxyzn,txyzn,mcele,1)
    end subroutine

    subroutine  yzsym(dxyzn,txyzn,mcele)
        implicit   none 

        integer  mcele
        real*8   dxyzn(100,3),txyzn(100,2,3)

        call symcore(dxyzn,txyzn,mcele,2)
    end subroutine

    subroutine  zxsym(dxyzn,txyzn,mcele)
        implicit   none 

        integer  mcele
        real*8   dxyzn(100,3),txyzn(100,2,3)

        call symcore(dxyzn,txyzn,mcele,3)
    end subroutine
    subroutine  yzzxsym(dxyzn,txyzn,mcele)
        implicit   none 

        !integer  i,j
        integer  mcele

        real*8   dxyzn(100,3),txyzn(100,2,3)
        !real(8) :: p(3,3),p2(3),p3(3)


        call zxsym(dxyzn,txyzn,mcele)
        call yzsym(dxyzn,txyzn,mcele)
    end subroutine
! *****************************************************************
!
!   对于三角形单元，计算各节点处边的切向量
!
! *****************************************************************
!
        SUBROUTINE TANGEL6(IELEM,JNEW,NNODE,NELEM,NCON,XYZ,TXYZ)
!	  USE MVAR_MOD
!	   
        IMPLICIT   NONE 
!	  
       INTEGER,INTENT(IN):: IELEM,JNEW
       INTEGER,INTENT(IN):: NNODE,NELEM
       INTEGER,INTENT(IN):: NCON(NELEM,8)
      
       REAL*8,INTENT(IN)::  XYZ(3,NNODE) 
	 REAL*8,INTENT(OUT):: TXYZ(2,3)

	  REAL*8   DX,DY,DZ,DR
!            
!
! ** Trianglar element **
!
        IF (JNEW .EQ. 1)  THEN
          DX=XYZ(1,NCON(IELEM,4))-XYZ(1,NCON(IELEM,1))
          DY=XYZ(2,NCON(IELEM,4))-XYZ(2,NCON(IELEM,1))
          DZ=XYZ(3,NCON(IELEM,4))-XYZ(3,NCON(IELEM,1))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,6))-XYZ(1,NCON(IELEM,1))
          DY=XYZ(2,NCON(IELEM,6))-XYZ(2,NCON(IELEM,1))
          DZ=XYZ(3,NCON(IELEM,6))-XYZ(3,NCON(IELEM,1))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        

        
        ELSE IF (JNEW .EQ. 2)  THEN
        
          DX=XYZ(1,NCON(IELEM,5))-XYZ(1,NCON(IELEM,2))
          DY=XYZ(2,NCON(IELEM,5))-XYZ(2,NCON(IELEM,2))
          DZ=XYZ(3,NCON(IELEM,5))-XYZ(3,NCON(IELEM,2))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,4))-XYZ(1,NCON(IELEM,2))
          DY=XYZ(2,NCON(IELEM,4))-XYZ(2,NCON(IELEM,2))
          DZ=XYZ(3,NCON(IELEM,4))-XYZ(3,NCON(IELEM,2))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        

        ELSE IF (JNEW .EQ. 3)  THEN
          DX=XYZ(1,NCON(IELEM,6))-XYZ(1,NCON(IELEM,3))
          DY=XYZ(2,NCON(IELEM,6))-XYZ(2,NCON(IELEM,3))
          DZ=XYZ(3,NCON(IELEM,6))-XYZ(3,NCON(IELEM,3))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,5))-XYZ(1,NCON(IELEM,3))
          DY=XYZ(2,NCON(IELEM,5))-XYZ(2,NCON(IELEM,3))
          DZ=XYZ(3,NCON(IELEM,5))-XYZ(3,NCON(IELEM,3))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        

        ELSE IF (JNEW .EQ. 4)  THEN
          DX=XYZ(1,NCON(IELEM,3))-XYZ(1,NCON(IELEM,4))
          DY=XYZ(2,NCON(IELEM,3))-XYZ(2,NCON(IELEM,4))
          DZ=XYZ(3,NCON(IELEM,3))-XYZ(3,NCON(IELEM,4))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,1))-XYZ(1,NCON(IELEM,4))
          DY=XYZ(2,NCON(IELEM,1))-XYZ(2,NCON(IELEM,4))
          DZ=XYZ(3,NCON(IELEM,1))-XYZ(3,NCON(IELEM,4))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 5)  THEN
          DX=XYZ(1,NCON(IELEM,2))-XYZ(1,NCON(IELEM,4))
          DY=XYZ(2,NCON(IELEM,2))-XYZ(2,NCON(IELEM,4))
          DZ=XYZ(3,NCON(IELEM,2))-XYZ(3,NCON(IELEM,4))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,3))-XYZ(1,NCON(IELEM,4))
          DY=XYZ(2,NCON(IELEM,3))-XYZ(2,NCON(IELEM,4))
          DZ=XYZ(3,NCON(IELEM,3))-XYZ(3,NCON(IELEM,4))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 6)  THEN
          DX=XYZ(1,NCON(IELEM,1))-XYZ(1,NCON(IELEM,5))
          DY=XYZ(2,NCON(IELEM,1))-XYZ(2,NCON(IELEM,5))
          DZ=XYZ(3,NCON(IELEM,1))-XYZ(3,NCON(IELEM,5))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,2))-XYZ(1,NCON(IELEM,5))
          DY=XYZ(2,NCON(IELEM,2))-XYZ(2,NCON(IELEM,5))
          DZ=XYZ(3,NCON(IELEM,2))-XYZ(3,NCON(IELEM,5))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 7)  THEN
          DX=XYZ(1,NCON(IELEM,3))-XYZ(1,NCON(IELEM,5))
          DY=XYZ(2,NCON(IELEM,3))-XYZ(2,NCON(IELEM,5))
          DZ=XYZ(3,NCON(IELEM,3))-XYZ(3,NCON(IELEM,5))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,1))-XYZ(1,NCON(IELEM,5))
          DY=XYZ(2,NCON(IELEM,1))-XYZ(2,NCON(IELEM,5))
          DZ=XYZ(3,NCON(IELEM,1))-XYZ(3,NCON(IELEM,5))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 8)  THEN
          DX=XYZ(1,NCON(IELEM,2))-XYZ(1,NCON(IELEM,6))
          DY=XYZ(2,NCON(IELEM,2))-XYZ(2,NCON(IELEM,6))
          DZ=XYZ(3,NCON(IELEM,2))-XYZ(3,NCON(IELEM,6))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,3))-XYZ(1,NCON(IELEM,6))
          DY=XYZ(2,NCON(IELEM,3))-XYZ(2,NCON(IELEM,6))
          DZ=XYZ(3,NCON(IELEM,3))-XYZ(3,NCON(IELEM,6))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR        
        ELSE IF (JNEW .EQ. 9)  THEN
          DX=XYZ(1,NCON(IELEM,1))-XYZ(1,NCON(IELEM,6))
          DY=XYZ(2,NCON(IELEM,1))-XYZ(2,NCON(IELEM,6))
          DZ=XYZ(3,NCON(IELEM,1))-XYZ(3,NCON(IELEM,6))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(1,1)=DX/DR
	    TXYZ(1,2)=DY/DR
	    TXYZ(1,3)=DZ/DR
!  
          DX=XYZ(1,NCON(IELEM,2))-XYZ(1,NCON(IELEM,6))
          DY=XYZ(2,NCON(IELEM,2))-XYZ(2,NCON(IELEM,6))
          DZ=XYZ(3,NCON(IELEM,2))-XYZ(3,NCON(IELEM,6))
          DR=SQRT(DX*DX+DY*DY+DZ*DZ)
      	  TXYZ(2,1)=DX/DR
	    TXYZ(2,2)=DY/DR
	    TXYZ(2,3)=DZ/DR              
!       
        ENDIF
        
!
!	  print *, TXYZN(IEL,LI,1),TXYZN(IEL,LI,2),TXYZN(IEL,LI,3)

       RETURN
       END           
