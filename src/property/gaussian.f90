module gaussian_info

    use mesh

    implicit none
    real(8),allocatable :: samb(:,:,:),sambxy(:,:,:),dsamb(:,:,:)

contains
    !C
    !C *********************************************************************
    !C *                                                                   *
    !C *  Identify the Gaussian sampling points and evaluate the           *
    !C *  corresponding values for the shape functions and Jacobian        *
    !C *  matrices                                                         *
    !C *                                                                   *
    !C *********************************************************************
    !C
    subroutine get_gaussian_data(xc,yc,zc)

        use mfunc_mod,only:spfunc8,spfunc6
        use proj_cnst,only:cross_product


        implicit   none  
        real(8),intent(in) :: xc,yc,zc
        integer ie,k,j,lk,li,lj,isi,ieta,nsamb,ind 
        real*8 xitsi(4),xiteta(4),wit(4),xiq(4),wiq(4)
        real*8 sf(8),dsf(2,8),xj(3,3) 

        real*8 det,det1,det2,det3,dum,si,eta 
        real(8) :: p(3)

        !C                      
        !C
        !C ** matrix XITSI store the Gauss-Legendre sampling points(3) for the 
        !C    triangular element in SI coordinate
        !C
        DATA XITSI/0.333333333333333D0, 0.6D0, 0.2D0, 0.2D0/
        !C
        !C ** matrix XITETA store the Gauss-Legendre sampling points(3) for the 
        !C    triangular element in ETA coordinate                         
        !C
        DATA XITETA/0.333333333333333D0, 0.2D0, 0.6D0, 0.2D0/

        !C
        !C ** matrix WIT store the Gauss-Legendre weighting factors for the 
        !C    triangular element in SI coordinate
        !C
        DATA WIT/-0.28125D0,.260416666666667D0,.260416666666667D0,&
            &          .260416666666667D0/
        !C
        !C
        !C ** matrix XIQ store the Gauss-Legendre sampling points(3) for the 
        !C    quadrilateral element 
        !c   
        DATA XIQ/-0.861136311594053D0,-0.339981043584856D0,&
            &          0.339981043584856D0, 0.861136311594053D0/ 
        !C
        !C ** matrix WIQ store the Gauss-Legendre weighting factors for the 
        !C    quadrilateral element 
        !c
        DATA WIQ/0.347854845137454D0,0.652145154862546D0,&
            &         0.652145154862546D0,0.347854845137454D0/

        allocate(samb(nelem,16,0:8),sambxy(nelem,16,3))
        allocate(dsamb(nelem,16,6))
        p=(/xc,yc,zc/)

        do  ie=1, nelem

            !c       NSAMB:  codes of sampling points inside an element
            nsamb=0
            !C ** Quadrilateral element **
            IF(NCN(IE).EQ.8) THEN
                do isi=1,4;do ieta=1,4
                    nsamb=nsamb+1

                    si =xiq(isi)
                    eta=xiq(ieta)
                    call spfunc8(si,eta,sf,dsf)
                    !C ** evaluate the Jacobian matrix at (SI,ETA),  XJ(2,3)
                    !c       LI: 1--SI,  2--ETA
                    !c       LJ: 1--X,   2--Y,   3--Z
                    xj(1:2,1:3) = matmul(dsf(1:2,1:ncn(ie)), &
                        & transpose(xyze(1:3,1:ncn(ie),ie)))
                    !** compute the determinant of the Jacobian maxtix at (SI,ETA), DET
                    det=norm2(cross_product(xj(1,:),xj(2,:)))
                    ! ** transform the local coordinates of the sampling points to 
                    !    global coordinates
                    !                                   
                    !       SAMBXY=??
                    !       I: code of the element,  LI:1--X, 2--Y, 3--Z  
                    !       NCON:code of the node in the global mesh
                    !       NSAMB:   
                    !       SF: shape function              

                    sambxy(ie,nsamb,1:3) = matmul(sf(1:ncn(ie)),transpose(xyze(1:3,1:ncn(ie),ie))) 
                    dsamb(ie,nsamb,1:3) = matmul(sf(1:ncn(ie)),transpose(dxyze(1:3,1:ncn(ie),ie)))
                    ! normlised normal vector
                    dsamb(ie,nsamb,1:3)=dsamb(ie,nsamb,1:3)/norm2(dsamb(ie,nsamb,1:3))
                    ! compute tangent vector
                    dsamb(ie,nsamb,4:6)=cross_product(sambxy(ie,nsamb,1:3)-p,&
                        & dsamb(ie,nsamb,1:3))
                    !C ** calculate the free surface boundary condition*WIT*DET
                    samb(ie,nsamb,0) = wiq(isi)*wiq(ieta)*det
                    samb(ie,nsamb,1:ncn(ie)) =sf(1:ncn(ie))*samb(ie,nsamb,0)
                end do; end do


            else
                !c ** for triangular element **
                print *,"triangle is depreciated"
                do  j=1, 4
                    nsamb=nsamb+1
                    !c
                    !c **  calculate the shape function at the sampling points
                    !c
                    si =xitsi(j)
                    eta=xiteta(j)
                    call spfunc6(si,eta,sf,dsf)

                end do
            endif
        end do



        do ind=1, nnoded
            dxyz(4:6,ind)=cross_product(xyz(1:3,ind)-p,dxyz(1:3,ind))
        enddo


      end subroutine
end module

