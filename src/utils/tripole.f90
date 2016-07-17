!-------------------------------------------------------------------------------
! DUTWAV
!-------------------------------------------------------------------------------
!  MODULE: tripole_transform
!
!> @brief
!! <Triple transform to prepare for weak singular integral>
!!
!! @author
!! DUT 
!!
!! @date
!! 07 Mar 2015
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Mar 2015 -- Rewrite tripole for clarity.
!
!-------------------------------------------------------------------------------
module tripole_transform
    use kinds
    implicit none
    
    integer :: nosamp
    real(rk) :: xynod(3,50),dxynod(6,50),samnod(50,0:8)

contains
    !! ===================================================
    !    !C  Tri-pole transformation for removing singularity
    !
    subroutine tripol(nodnum,ncne,xyzt,dxyzt,xc,yc,zc)
        use shape_funcs
        use linalg,only:cross_product

        implicit none

        integer,intent(in) :: nodnum,ncne
        real(rk),intent(in) :: xyzt(3,8),dxyzt(3,8),xc,zc,yc

        integer :: ixrho,iyrho,notri,i,ub
        real(rk) :: xiq(4),wiq(4)
        real(rk) :: sf(8),dsf(2,8),xj(2,3)
        real(rk) :: sitri(3,3),etatri(3,3)
        real(rk) :: si,eta,xrho,yrho,ag,det,dum,dum1,dum2,dum3 
        real(rk) :: p(3),np(3),tp(3) 

        !C ** matrix XIQ store the Gauss-Legendre sampling points (4) for the 
        !C    quadrilateral element 

        DATA XIQ/-0.861136311594053D0,-0.339981043584856D0,&
                  0.339981043584856D0, 0.861136311594053D0/
        !C ** matrix WIQ store the Gauss-Legendre weighting factors for the 
        !C    quadrilateral element 
        DATA WIQ/0.347854845137454D0,0.652145154862546D0,&
                 0.652145154862546D0,0.347854845137454D0/
             
        !C ** identity the local nodal position of INODE IN IELEM

        if(ncne.eq.8) then
            if(nodnum.eq.1.or.nodnum.eq.3.or.&
                nodnum.eq.5.or.nodnum.eq.7) then
                notri=2
            else
                notri=3
            endif
        else
            if(nodnum.le.3)then
                notri=1
            else
                notri=2
            endif
        endif

        !C ** identity the local coordinate(SI,ETA) of the trianges

        if(ncne.eq.8) then

            !c ** quadrilateral element
            if(nodnum.eq.1) then

                sitri(1:2,1:3)  =  1.0d0 
                sitri(3,1:3)    =  0.0d0
                sitri(1,1)      = -1.0d0
                sitri(1,3)      = -1.0d0
                sitri(2,1)      = -1.0d0
                
                etatri(1:2,1:3)  =  1.0d0 
                etatri(3,1:3)    =  0.0d0
                etatri(1,1)      = -1.0d0
                etatri(2,1)      = -1.0d0
                etatri(2,2)      = -1.0d0

            else if(nodnum.eq.3) then

                sitri(1:2,1:3)  =  1.0d0 
                sitri(3,1:3)    =  0.0d0
                sitri(1,2)      = -1.0d0
                sitri(1,3)      = -1.0d0
                sitri(2,3)      = -1.0d0

                etatri(1:2,1:3)  =  1.0d0 
                etatri(3,1:3)    =  0.0d0
                etatri(1,1)      = -1.0d0
                etatri(1,3)      = -1.0d0
                etatri(2,1)      = -1.0d0


            else if(nodnum.eq.5) then

                sitri(1:2,1:3)  =  1.0d0 
                sitri(3,1:3)    =  0.0d0
                sitri(1,2)      = -1.0d0
                sitri(1,3)      = -1.0d0
                sitri(2,2)      = -1.0d0

                etatri(1:2,1:3)  =  1.0d0 
                etatri(3,1:3)    =  0.0d0
                etatri(1,3)      = -1.0d0
                etatri(2,2)      = -1.0d0
                etatri(2,3)      = -1.0d0

            else if(nodnum.eq.7) then

                sitri(1:2,1:3)  =  1.0d0 
                sitri(3,1:3)    =  0.0d0
                sitri(1,1)      = -1.0d0
                sitri(1,2)      = -1.0d0
                sitri(2,1)      = -1.0d0

                etatri(1:2,1:3)  =  1.0d0 
                etatri(3,1:3)    =  0.0d0
                etatri(1,2)      = -1.0d0
                etatri(1,3)      = -1.0d0
                etatri(2,2)      = -1.0d0


            else if(nodnum.eq.2) then

                sitri(1:3,2:3)  =  1.0d0 
                sitri(1:3,1)    =  0.0d0
                sitri(1,2)      = -1.0d0
                sitri(1,3)      = -1.0d0
                sitri(2,3)      = -1.0d0

                etatri(1:3,1:3)  =  1.0d0 
                etatri(1,1)      = -1.0d0
                etatri(1,3)      = -1.0d0
                etatri(2,1)      = -1.0d0
                etatri(3,1)      = -1.0d0

            ELSE IF(NODNUM.EQ.4) THEN
                
                sitri(1:3,1:3)  =  1.0d0 
                sitri(1,3)      = -1.0d0
                sitri(2,2)      = -1.0d0
                sitri(2,3)      = -1.0d0
                sitri(3,2)      = -1.0d0

                etatri(1:3,2:3)  =  1.0d0 
                etatri(1:3,1)    =  0.0d0 
                etatri(2,3)      = -1.0d0
                etatri(3,2)      = -1.0d0
                etatri(3,3)      = -1.0d0

            else if(nodnum.eq.6) then

                sitri(1:3,2:3)  =  1.0d0 
                sitri(1:3,1)    =  0.0d0 
                sitri(1,2)      = -1.0d0
                sitri(1,3)      = -1.0d0
                sitri(2,2)      = -1.0d0

                etatri(1:3,1:3)  =  1.0d0 
                etatri(1,3)      = -1.0d0
                etatri(2,2)      = -1.0d0
                etatri(2,3)      = -1.0d0
                etatri(3,2)      = -1.0d0


            ELSE IF(NODNUM.EQ.8) THEN
                ! why 5  -1
                sitri(1:3,1:3)  =  1.0d0 
                sitri(1,1)      = -1.0d0
                sitri(1,3)      = -1.0d0
                sitri(2,1)      = -1.0d0
                sitri(3,1)      = -1.0d0
                sitri(3,2)      = -1.0d0

                etatri(1:3,1)    =  0.0d0 
                etatri(1:3,2:3)  =  1.0d0 
                etatri(2,2)      = -1.0d0
                etatri(3,3)      = -1.0d0
                etatri(3,2)      = -1.0d0

            endif

        else
            !C
            !C ** triangular element ***************
            !C          
            IF(NODNUM.EQ.1) THEN

                sitri(1:3,1:3)  =  0.0d0
                sitri(1,2)      =  1.0d0

                etatri(1:3,1:3) =  0.0d0 
                etatri(1,3)     =  1.0d0 

            ELSE IF(NODNUM.EQ.2) THEN

                sitri(1:3,1:3)  =  0.0d0
                sitri(1,1)      =  1.0d0

                etatri(1:3,1:3) =  0.0d0 
                etatri(1,2)     =  1.0d0 

            ELSE IF(NODNUM.EQ.3) THEN

                sitri(1:3,1:3)  =  0.0d0
                sitri(1,3)      =  1.0d0

                etatri(1:3,1:3) =  0.0d0 
                etatri(1,1)     =  1.0d0 


            ELSE IF(NODNUM.EQ.4) THEN
                
                sitri(1:3,1:3)  =  0.0d0
                sitri(1:2,1)    =  0.5d0
                sitri(2,2)      =  1.0d0

                etatri(1:3,1:3) =  0.0d0 
                etatri(1,2)     =  1.0d0 
                etatri(2,3)     =  1.0d0 

            ELSE IF(NODNUM.EQ.5) THEN
                
                sitri(1:3,1:3)  =  0.0d0
                sitri(1:2,1)    =  0.5d0
                sitri(1,3)      =  1.0d0

                etatri(1:3,1:3) =  0.0d0 
                etatri(1:2,1)   =  0.5d0 
                etatri(2,2)     =  1.0d0 

            ELSE IF(NODNUM.EQ.6) THEN

                sitri(1:3,1:3)  =  0.0d0
                sitri(2,2)      =  1.0d0
                sitri(1,3)      =  1.0d0

                etatri(1:3,1:3) =  0.0d0 
                etatri(1:2,1)   =  0.5d0 
                etatri(2,3)     =  1.0d0 

            ENDIF

        ENDIF

        nosamp=0
        do i=1,notri

            if(ncne.eq.8) then
                ag=2.0d0
                if(notri.eq.3) then
                    if(i.eq.1.or.i.eq.3) ag=1.0d0
                endif
            else
                ag=0.5d0
                if(notri.eq.2) ag=0.25d0
            endif

            do ixrho=1,4;do iyrho=1,4
                nosamp=nosamp+1
                !c
                !c **  calculate the shape function at the sampling points(xrho,yrho)
                !c
                xrho = xiq(ixrho)
                yrho = xiq(iyrho)

                dum1=0.5d0*(1.0d0-xrho)
                dum2=0.25d0*(1.0d0+xrho)*(1.0d0-yrho)
                dum3=0.25d0*(1.0d0+xrho)*(1.0d0+yrho)

                si =dum1*sitri (i,1)+dum2*sitri (i,2)+dum3*sitri (i,3)
                eta=dum1*etatri(i,1)+dum2*etatri(i,2)+dum3*etatri(i,3)

                if(ncne.eq.8) then
                    call spfunc8(si,eta,sf,dsf)
                    ! @var ub :  upper bound
                    ub=8
                else
                    call spfunc6(si,eta,sf,dsf)
                    ub=6
                endif
                !c
                !c ** evaluate the jacobian matrix at (si,eta)

                xj=matmul(dsf(1:2,1:ub),transpose(xyzt(1:3,1:ub)))
                det = norm2(cross_product(xj(1,:),xj(2,:)))
                !< global pos
                p = matmul(sf(1:ub),transpose(xyzt(1:3,1:ub)))
                !< norm vector
                np = matmul(sf(1:ub),transpose(dxyzt(1:3,1:ub)))
                np = np/norm2(np)
                !< tagent vector
                tp = cross_product(p-(/xc,yc,zc/),np)
                tp = tp/norm2(tp)

                xynod(1:3,nosamp)=p
                dxynod(1:3,nosamp)=np
                dxynod(4:6,nosamp) =tp

                samnod(nosamp,0)=0.25d0*ag*(1.0+xrho)*wiq(ixrho)*wiq(iyrho)*det
                samnod(nosamp,1:ub)=sf(1:ub)*samnod(nosamp,0)

            end do;end do
        end do
    end subroutine tripol
end module
