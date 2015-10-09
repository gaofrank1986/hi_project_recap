        subroutine compute_coeff_b(ndim,nf,lamda,node,npowg,npowf,xp,xip,&
                                & xiq,coefg,coefb)
      
        implicit none

        integer,intent(in) ::  node,npowg,npowf,ndim,nf
        real(8),intent(in) :: lamda,xp(ndim),xip(ndim - 1)&
                    & ,xiq(ndim - 1), coefg(0:npowg)

        real(8),intent(out) :: coefb(0:11,nf)

        real(8):: ri(ndim), cosn(ndim), gcd(3,ndim-1)

        real(8) :: rho_q,slop(ndim - 1)

        real(8) :: drdx(ndim),xi(ndim - 1),x(ndim)&
                    &,rmat(npowf,npowf), &
                    &   sf_iter(node),&
                    &   fq(nf),a(ndim)      
        integer :: i,j,step_n,m,jp
        integer :: nbdm
        real(8) :: rho_step,fjcb,robar
        real(8) :: rho,r2,r,gm,drdn

        NBDM = num_dim - 1

        slop = xiq - xip
        rho_q = norm2(slop)
        slop = slop/rho_q ! normalized vector, cos(theta),sin(theta)

        rho_step=rho_q/dble(npowf) ! divide rho_q to npowf parts

        do 20 step_n=0,npowf ! this the k value in Gao,XW note 

            rho=rho_step*dble(step_n) 
            xi=xip+rho*slop !!!!!!!!!!!!!!!!! xi updated here!!!!!!!!!!!!!!!!!!
            robar=coefg(0)

            do m=1,npowg
                 robar=robar+coefg(m)*rho**m
            enddo

            ROBAR=DSQRT(ROBAR)       ! Eq.(3-6-28)

            call shapef(ndim,node,cnr_lcl_mtx,cnr_glb_mtx,xi,xp,ri,sf_iter)
            !shape func based on xi
            !ri give glb vector from xp(it is src_glb) to xi (update along rho)
            R = norm2(RI)

            X=XP+RI

            call dshape(ndim,node,cnr_lcl_mtx,cnr_glb_mtx,xi,cosn,fjcb,gcd)
            ! GCD gives the two tangent vector
            ! cosn gives a normalized normal vector
            ! dshape based on xi

            if(rho.gt.1.0d-10)then    
                drdx=ri/r
            else
                a = 0.d0
     !             do 10 i=1,ndim; a(i)=0.d0
     !             do 10 j=1,nbdm        
                forall (i = 1:ndim)    
                    a(i)=a(i)+dot_product(gcd(i,1:nbdm),slop(1:nbdm))
                end forall                  
                gm=dsqrt(dot_product(a,a))
                drdx=a/gm                 ! eq.(3-6-74)
            endif 

            ! dr/dx is defined above
            drdn = dot_product(cosn,drdx)  !!!!  notice dr/dn is defined  here
            !CALL F_BAR(NDIM,NBDM,DRDX,COSN,R,DRDN,XI,SF_iter,XP,X,NF,FQ)
            call f_integrand(ndim,nf,cosn,drdx,drdn,sf_iter,fq)
            COEFB(step_n,:) = FQ*FJCB/ROBAR**lamda

     
            IF(step_n.EQ.0) GOTO 20
            COEFB(step_n,:)=(COEFB(step_n,:)-COEFB(0,:))/RHO

            RMAT(step_n,1)=1.D0
            DO JP=2,NPOWF; RMAT(step_n,JP)=RMAT(step_n,JP-1)*RHO; ENDDO

20      CONTINUE
        


        CALL INVSOLVR(NPOWF,NPOWF,RMAT,NPOWF,-1)
   

        forall (i = 1:npowf,j=1:NF)

        coefb(i,j) = DOT_PRODUCT(rmat(i,:),coefb(1:npowf,j))

        end forall
        !COEFB(1:NPOWF,:)=MATMUL(RMAT,COEFB(1:NPOWF,:))    ! Bk IN Eq.(3-6-62)
        ! coefb has size 12xnum_intgd,but only from index 1 to 7 has useful info
      
        

     end subroutine compute_coeff_B
