module tripole_transform
    use kinds
    implicit none
    
    integer :: nosamp
    real(rk) :: xynod(3,50),dxynod(6,5),samnod(50,0:8)

contains
    !! ===================================================
    !    !C  Tri-pole transformation for removing singularity
    !
    subroutine tripol(nodnum,ncne,xyzt,dxyzt,xc,yc,zc)
        !use mesh,only: 
        use shape_funcs
        use proj_cnst,only:cross_product
        !use body_property,only:xc,yc,zc

        implicit none

        integer,intent(in) :: nodnum,ncne
        real(rk),intent(in) :: xyzt(3,8),dxyzt(3,8),xc,zc,yc

        integer :: ixrho,iyrho,notri,i
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
                !sitri (1,1)=-1.0d0
                !sitri (1,2)= 1.0d0
                !sitri (1,3)=-1.0d0

                !sitri (2,1)=-1.0d0
                !sitri (2,2)= 1.0d0
                !sitri (2,3)= 1.0d0

                !sitri (3,1)= 0.0d0
                !sitri (3,2)= 0.0d0
                !sitri (3,3)= 0.0d0

                !ETATRI(1,1)=-1.0D0
                !ETATRI(1,2)= 1.0D0
                !ETATRI(1,3)= 1.0D0
                !C
                !ETATRI(2,1)=-1.0D0
                !ETATRI(2,2)=-1.0D0
                !ETATRI(2,3)= 1.0D0
                !C
                !ETATRI(3,1)= 0.0D0
                !ETATRI(3,2)= 0.0D0
                !ETATRI(3,3)= 0.0D0
                !C
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

                !C
                !SITRI (1,1)= 1.0D0
                !SITRI (1,2)=-1.0D0
                !SITRI (1,3)=-1.0D0
                !C
                !SITRI (2,1)= 1.0D0
                !SITRI (2,2)= 1.0D0
                !SITRI (2,3)=-1.0D0
                !C
                !SITRI (3,1)= 0.0D0
                !SITRI (3,2)= 0.0D0
                !SITRI (3,3)= 0.0D0
                !C
                !ETATRI(1,1)=-1.0D0
                !ETATRI(1,2)= 1.0D0
                !ETATRI(1,3)=-1.0D0
                !C
                !ETATRI(2,1)=-1.0D0
                !ETATRI(2,2)= 1.0D0
                !ETATRI(2,3)= 1.0D0
                !C
                !ETATRI(3,1)= 0.0D0
                !ETATRI(3,2)= 0.0D0
                !ETATRI(3,3)= 0.0D0
                !C
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
                !C
                !SITRI (1,1)= 1.0D0
                !SITRI (1,2)=-1.0D0
                !SITRI (1,3)=-1.0D0
                !C
                !SITRI (2,1)= 1.0D0
                !SITRI (2,2)=-1.0D0
                !SITRI (2,3)= 1.0D0
                !C
                !SITRI (3,1)= 0.0D0
                !SITRI (3,2)= 0.0D0
                !SITRI (3,3)= 0.0D0
                !C
                !ETATRI(1,1)= 1.0D0
                !ETATRI(1,2)= 1.0D0
                !ETATRI(1,3)=-1.0D0
                !C
                !ETATRI(2,1)= 1.0D0
                !ETATRI(2,2)=-1.0D0
                !ETATRI(2,3)=-1.0D0
                !C
                !ETATRI(3,1)= 0.0D0
                !ETATRI(3,2)= 0.0D0
                !ETATRI(3,3)= 0.0D0
                !C
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

                !C
                !SITRI (1,1)=-1.0D0
                !SITRI (1,2)=-1.0D0
                !SITRI (1,3)= 1.0D0
                !C
                !SITRI (2,1)=-1.0D0
                !SITRI (2,2)= 1.0D0
                !SITRI (2,3)= 1.0D0
                !C
                !SITRI (3,1)= 0.0D0
                !SITRI (3,2)= 0.0D0
                !SITRI (3,3)= 0.0D0
                !C
                !ETATRI(1,1)= 1.0D0
                !ETATRI(1,2)=-1.0D0
                !ETATRI(1,3)=-1.0D0
                !C
                !ETATRI(2,1)= 1.0D0
                !ETATRI(2,2)=-1.0D0
                !ETATRI(2,3)= 1.0D0
                !C
                !ETATRI(3,1)= 0.0D0
                !ETATRI(3,2)= 0.0D0
                !ETATRI(3,3)= 0.0D0
                !C
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
                !C
                !SITRI (1,1)= 0.0D0
                !SITRI (1,2)=-1.0D0
                !SITRI (1,3)=-1.0D0
                !C
                !SITRI (2,1)= 0.0D0
                !SITRI (2,2)= 1.0D0
                !SITRI (2,3)=-1.0D0
                !C
                !SITRI (3,1)= 0.0D0
                !SITRI (3,2)= 1.0D0
                !SITRI (3,3)= 1.0D0
                !C
                !ETATRI(1,1)=-1.0D0
                !ETATRI(1,2)= 1.0D0
                !ETATRI(1,3)=-1.0D0
                !C
                !ETATRI(2,1)=-1.0D0
                !ETATRI(2,2)= 1.0D0
                !ETATRI(2,3)= 1.0D0
                !C
                !ETATRI(3,1)=-1.0D0
                !ETATRI(3,2)=-1.0D0
                !ETATRI(3,3)= 1.0D0
                !C
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
                !C
                !SITRI (1,1)= 1.0D0
                !SITRI (1,2)= 1.0D0
                !SITRI (1,3)=-1.0D0
                !C
                !SITRI (2,1)= 1.0D0
                !SITRI (2,2)=-1.0D0
                !SITRI (2,3)=-1.0D0
                !C
                !SITRI (3,1)= 1.0D0
                !SITRI (3,2)=-1.0D0
                !SITRI (3,3)= 1.0D0
                !C
                !ETATRI(1,1)= 0.0D0
                !ETATRI(1,2)= 1.0D0
                !ETATRI(1,3)= 1.0D0
                !C
                !ETATRI(2,1)= 0.0D0
                !ETATRI(2,2)= 1.0D0
                !ETATRI(2,3)=-1.0D0
                !C
                !ETATRI(3,1)= 0.0D0
                !ETATRI(3,2)=-1.0D0
                !ETATRI(3,3)=-1.0D0
                !C
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

                !C
                !SITRI (1,1)= 0.0D0
                !SITRI (1,2)=-1.0D0
                !SITRI (1,3)=-1.0D0
                !C
                !SITRI (2,1)= 0.0D0
                !SITRI (2,2)=-1.0D0
                !SITRI (2,3)= 1.0D0
                !C
                !SITRI (3,1)= 0.0D0
                !SITRI (3,2)= 1.0D0
                !SITRI (3,3)= 1.0D0
                !C
               !ETATRI(1,1)= 1.0D0
                !ETATRI(1,2)= 1.0D0
                !ETATRI(1,3)=-1.0D0
                !C
                !ETATRI(2,1)= 1.0D0
                !ETATRI(2,2)=-1.0D0
                !ETATRI(2,3)=-1.0D0
                !C
                !ETATRI(3,1)= 1.0D0
                !ETATRI(3,2)=-1.0D0
                !ETATRI(3,3)= 1.0D0
                !C
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

                !C
                !SITRI (1,1)=-1.0D0
                !SITRI (1,2)= 1.0D0
                !SITRI (1,3)=-1.0D0
                !C
                !SITRI (2,1)=-1.0D0
                !SITRI (2,2)= 1.0D0
                !SITRI (2,3)= 1.0D0
                !C
                !SITRI (3,1)=-1.0D0
                !SITRI (3,2)=-1.0D0
                !SITRI (3,3)= 1.0D0
                !C
                !ETATRI(1,1)= 0.0D0
                !ETATRI(1,2)= 1.0D0
                !ETATRI(1,3)= 1.0D0
                !C
                !ETATRI(2,1)= 0.0D0
                !ETATRI(2,2)=-1.0D0
                !ETATRI(2,3)= 1.0D0
                !C
                !ETATRI(3,1)= 0.0D0
                !ETATRI(3,2)=-1.0D0
                !ETATRI(3,3)=-1.0D0
                !C
            endif

        ELSE
            !C
            !C ** triangular element ***************
            !C          
            IF(NODNUM.EQ.1) THEN
                
                SITRI (1,1)= 0.0D0
                SITRI (1,2)= 1.0D0
                SITRI (1,3)= 0.0D0
                
                SITRI (2,1)= 0.0D0
                SITRI (2,2)= 0.0D0
                SITRI (2,3)= 0.0D0
                
                SITRI (3,1)= 0.0D0
                SITRI (3,2)= 0.0D0
                SITRI (3,3)= 0.0D0
                
                ETATRI(1,1)= 0.0D0
                ETATRI(1,2)= 0.0D0
                ETATRI(1,3)= 1.0D0

                ETATRI(2,1)= 0.0D0
                ETATRI(2,2)= 0.0D0
                ETATRI(2,3)= 0.0D0

                ETATRI(3,1)= 0.0D0
                ETATRI(3,2)= 0.0D0
                ETATRI(3,3)= 0.0D0

            ELSE IF(NODNUM.EQ.2) THEN

                SITRI (1,1)= 1.0D0
                SITRI (1,2)= 0.0D0
                SITRI (1,3)= 0.0D0

                SITRI (2,1)= 0.0D0
                SITRI (2,2)= 0.0D0
                SITRI (2,3)= 0.0D0

                SITRI (3,1)= 0.0D0
                SITRI (3,2)= 0.0D0
                SITRI (3,3)= 0.0D0

                ETATRI(1,1)= 0.0D0
                ETATRI(1,2)= 1.0D0
                ETATRI(1,3)= 0.0D0

                ETATRI(2,1)= 0.0D0
                ETATRI(2,2)= 0.0D0
                ETATRI(2,3)= 0.0D0

                ETATRI(3,1)= 0.0D0
                ETATRI(3,2)= 0.0D0
                ETATRI(3,3)= 0.0D0

            ELSE IF(NODNUM.EQ.3) THEN

                SITRI (1,1)= 0.0D0
                SITRI (1,2)= 0.0D0
                SITRI (1,3)= 1.0D0

                SITRI (2,1)= 0.0D0
                SITRI (2,2)= 0.0D0
                SITRI (2,3)= 0.0D0

                SITRI (3,1)= 0.0D0
                SITRI (3,2)= 0.0D0
                SITRI (3,3)= 0.0D0

                ETATRI(1,1)= 1.0D0
                ETATRI(1,2)= 0.0D0
                ETATRI(1,3)= 0.0D0

                ETATRI(2,1)= 0.0D0
                ETATRI(2,2)= 0.0D0
                ETATRI(2,3)= 0.0D0

                ETATRI(3,1)= 0.0D0
                ETATRI(3,2)= 0.0D0
                ETATRI(3,3)= 0.0D0

            ELSE IF(NODNUM.EQ.4) THEN

                SITRI (1,1)= 0.5D0
                SITRI (1,2)= 0.0D0
                SITRI (1,3)= 0.0D0

                SITRI (2,1)= 0.5D0
                SITRI (2,2)= 1.0D0
                SITRI (2,3)= 0.0D0

                SITRI (3,1)= 0.0D0
                SITRI (3,2)= 0.0D0
                SITRI (3,3)= 0.0D0

                ETATRI(1,1)= 0.0D0
                ETATRI(1,2)= 1.0D0
                ETATRI(1,3)= 0.0D0

                ETATRI(2,1)= 0.0D0
                ETATRI(2,2)= 0.0D0
                ETATRI(2,3)= 1.0D0

                ETATRI(3,1)= 0.0D0
                ETATRI(3,2)= 0.0D0
                ETATRI(3,3)= 0.0D0

            ELSE IF(NODNUM.EQ.5) THEN

                SITRI (1,1)= 0.5D0
                SITRI (1,2)= 0.0D0
                SITRI (1,3)= 1.0D0

                SITRI (2,1)= 0.5D0
                SITRI (2,2)= 0.0D0
                SITRI (2,3)= 0.0D0

                SITRI (3,1)= 0.0D0
                SITRI (3,2)= 0.0D0
                SITRI (3,3)= 0.0D0

                ETATRI(1,1)= 0.5D0
                ETATRI(1,2)= 0.0D0
                ETATRI(1,3)= 0.0D0

                ETATRI(2,1)= 0.5D0
                ETATRI(2,2)= 1.0D0
                ETATRI(2,3)= 0.0D0

                ETATRI(3,1)= 0.0D0
                ETATRI(3,2)= 0.0D0
                ETATRI(3,3)= 0.0D0

            ELSE IF(NODNUM.EQ.6) THEN

                SITRI (1,1)= 0.0D0
                SITRI (1,2)= 0.0D0
                SITRI (1,3)= 1.0D0

                SITRI (2,1)= 0.0D0
                SITRI (2,2)= 1.0D0
                SITRI (2,3)= 0.0D0

                SITRI (3,1)= 0.0D0
                SITRI (3,2)= 0.0D0
                SITRI (3,3)= 0.0D0

                ETATRI(1,1)= 0.5D0
                ETATRI(1,2)= 0.0D0
                ETATRI(1,3)= 0.0D0

                ETATRI(2,1)= 0.5D0
                ETATRI(2,2)= 0.0D0
                ETATRI(2,3)= 1.0D0

                ETATRI(3,1)= 0.0D0
                ETATRI(3,2)= 0.0D0
                ETATRI(3,3)= 0.0D0

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
                else
                    call spfunc6(si,eta,sf,dsf)
                endif
                !c
                !c ** evaluate the jacobian matrix at (si,eta)
                !c
                !do 130 li=1,2
                !do 130 lj=1,3
                !dum=0.0d0
                !do 140 lk=1,ncne
                !dum=dum+dsf(li,lk)*xyzt(lj,lk)
                !140   	continue
                !130   	xj(li,lj)=dum
                xj=matmul(dsf(1:2,1:8),transpose(xyzt(1:3,1:8)))
                det = norm2(cross_product(xj(1,:),xj(2,:)))

                p = matmul(sf(1:8),transpose(xyzt(1:3,1:8)))
                !< norm vector
                np = matmul(sf(1:8),transpose(dxyzt(1:3,1:8)))
                np = np/norm2(np)
                !< tagent vector
                tp = cross_product(p-(/xc,yc,zc/),np)
                tp = tp/norm2(tp)


                xynod(1:3,nosamp)=p
                dxynod(1:3,nosamp)=np
                dxynod(4:6,nosamp) =tp

                samnod(nosamp,0)=0.25d0*ag*(1.0+xrho)*wiq(ixrho)*wiq(iyrho)*det
                samnod(nosamp,1:8)=sf(1:8)*samnod(nosamp,0)


                !c
                !c ** compute the determinant of the jacobian maxtix at (xrho,yrho)
                !c
                !det1=xj(1,2)*xj(2,3)-xj(1,3)*xj(2,2) 
                !det2=xj(1,1)*xj(2,3)-xj(1,3)*xj(2,1) 
                !det3=xj(1,1)*xj(2,2)-xj(1,2)*xj(2,1) 
                !det=dsqrt(det1*det1+det2*det2+det3*det3)

                !X=0.0D0
                !Y=0.0D0
                !Z=0.0D0

                !DX=0.0D0
                !DY=0.0D0
                !DZ=0.0D0

                !DO 150 LK=1, NCNE
                !X=X+SF(LK)*XYZT(1,LK)
                !Y=Y+SF(LK)*XYZT(2,LK)
                !Z=Z+SF(LK)*XYZT(3,LK)
                !DX=DX+SF(LK)*DXYZT(1,LK)
                !DY=DY+SF(LK)*DXYZT(2,LK)
                !DZ=DZ+SF(LK)*DXYZT(3,LK)
                !150   	CONTINUE
                !C
                !XYNOD(1,NOSAMP)=X
                !XYNOD(2,NOSAMP)=Y
                !XYNOD(3,NOSAMP)=Z
                !DXYNOD(1,NOSAMP)=DX
                !DXYNOD(2,NOSAMP)=DY
                !DXYNOD(3,NOSAMP)=DZ 
                !tmp = norm2(DXYNOD(1:3,NOSAMP))
                !DXYNOD(1:3,NOSAMP)=DXYNOD(1:3,NOSAMP)/tmp

                !DXYNOD(4,NOSAMP)=(Y-YC)*DZ-(Z-ZC)*DY
                !DXYNOD(5,NOSAMP)=(Z-ZC)*DX-(X-XC)*DZ
                !DXYNOD(6,NOSAMP)=(X-XC)*DY-(Y-YC)*DX
                !C
                !SAMNOD(NOSAMP,0)=0.25D0*AG*(1.0+XRHO)*WIQ(IXRHO)*WIQ(IYRHO)*DET
                !C
                !DO 160 J=1,NCNE
                !SAMNOD(NOSAMP,J)=SF(J)*SAMNOD(NOSAMP,0)
                !160   	CONTINUE
                !C
            end do;end do
        end do
    end subroutine tripol
end module
