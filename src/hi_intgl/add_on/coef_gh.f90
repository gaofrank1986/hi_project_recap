
     
    subroutine compute_coeff_gh(ndim,nbdm,npw,node,npowg,xp,xip,xiq,  &
     &                    coefg,coefh)
        implicit none

        integer,intent(in) ::  ndim,nbdm,node,npw,npowg
        real(8),intent(in) :: xp(ndim),xip(nbdm),xiq(nbdm)
        real(8),intent(out) :: coefg(0:npowg),coefh(0:npw)
        real(8) :: drdn

        real(8) ::cosn(ndim),gcd(ndim,nbdm),xi(ndim),ri(ndim),shap(node),     &
        &         coefc(0:npw)

        integer :: i,n,ip,jp

        call compute_coeff_g(ndim,nbdm,node,npowg,xp,xip,xiq,    &
        &           coefg)
        ! determine coefficients cn using eq. (3-6-37)

        coefc(0)=dsqrt(coefg(0))
        do n=1,npw
            coefc(n)=0.d0
            do i=1,n-1
                coefc(n)=coefc(n)-coefc(i)*coefc(n-i)     
            enddo
            if(n.le.npowg)coefc(n)=coefg(n)+coefc(n)   ! eq. (3-6-37)
            coefc(n)=coefc(n)/(2.d0*coefc(0))
        enddo
        ! determine coefficients hi using eq.(3-6-47a)
        coefc(1:npw)=coefc(1:npw)/coefc(0)    ! eq.(3-6-47b)
        coefh(0)=1.d0/coefc(0)
        coefh(1)=coefc(1)
        coefh(2)=2.*coefc(2)-coefc(1)*coefc(1)
        coefh(3)=3.*coefc(3)-3.*coefc(1)*coefc(2)+coefc(1)**3
        coefh(4)=4.*coefc(4)+4.*coefc(1)*coefc(1)*coefc(2)                &
        &        -4.*coefc(1)*coefc(3)-2.*coefc(2)*coefc(2)                &
        &        -coefc(1)**4
   
    end subroutine
    subroutine compute_coeff_g(ndim,nbdm,node,npowg,xp,xip,xiq,         &
     &                  coefg)
      ! THIS ROUTINE DETERMINES COEFFICIENTS Gm USING Eqs.(3-6-58)~(3-6-60)
        implicit none 

        integer,intent(in) ::  ndim,nbdm,node,npowg
        real(8),intent(in) :: xp(ndim),xip(nbdm),xiq(nbdm)
        real(8) :: drdn,xbar2,xxbar,vk,rhoq,rho,r2,fjcb

        real(8) ::     cosn(ndim),        &
            &          gcd(ndim,nbdm),xi(ndim),slop(nbdm),ri(ndim),     &
            &          shap(node),rmat(npowg,npowg),coefg(0:npowg)

        integer :: ip,jp


      IF(NODE.EQ.2) THEN         ! Eq.(24)
!        COEFG(0)=0.25*((cnr_glb_mtx(1,1)-cnr_glb_mtx(1,2))**2+(cnr_glb_mtx(2,1)-cnr_glb_mtx(2,2))**2)  
!       ELSEIF(NODE.EQ.3) THEN             ! FORM Eq.(3-6-33)
!        RI=cnr_glb_mtx(:,1)+cnr_glb_mtx(:,2)-2.*cnr_glb_mtx(:,3)     ! Xbar
!        XI=cnr_glb_mtx(:,1)-cnr_glb_mtx(:,2)                ! X1i-X2i 
!        XBAR2=RI(1)*RI(1)+RI(2)*RI(2)     ! Xbar^2
!        XXBAR=XI(1)*RI(1)+XI(2)*RI(2)
!        COEFG(0)=0.25*(XI(1)*XI(1)+XI(2)*XI(2))-XXBAR*XIP(1)             &
!      &         +XBAR2*XIP(1)*XIP(1)
!        COEFG(1)=(-0.5*XXBAR+XBAR2*XIP(1))*XIQ(1)
!        COEFG(2)=0.25*XBAR2
        ELSE      ! QUADRILATERAL ELEMENTS
            slop=xiq-xip
            !rhoq=dsqrt(dot_product(slop,slop))
            rhoq = norm2(slop)
            slop=slop/rhoq
            vk=rhoq/dble(npowg)
            do 20 ip=0,npowg
            rho=vk*dble(ip)
            xi(1:nbdm)=xip+rho*slop
            call shapef(ndim,node,cnr_lcl_mtx,cnr_glb_mtx,xi,xp,ri,shap)
            if(ip.ne.0) goto 10
            call dshape(ndim,node,cnr_lcl_mtx,cnr_glb_mtx,xi,cosn,fjcb,gcd)   !dx/dxi
            ri=matmul(gcd,slop); coefg(0)=dot_product(ri,ri)  ! eq.(3-6-58)
            goto 20
!  10    coefg(ip)=((dsqrt(r2)/rho)**2-coefg(0))/rho     ! eq.(3-6-60)

 10    coefg(ip)=((norm2(ri)/rho)**2-coefg(0))/rho     ! eq.(3-6-60)

       rmat(ip,1)=1.d0        ! eq.(3-6-42)
       do jp=2,npowg; rmat(ip,jp)=rmat(ip,jp-1)*rho; enddo
 20    continue
       
       call invsolvr(npowg,npowg,rmat,npowg,-1)       ! inverse matrix [r] in eq.(3-6-59)
       coefg(1:npowg)=matmul(rmat,coefg(1:npowg))     ! solving eq.(3-6-59) for {g}
        endif
    end subroutine compute_coeff_g
