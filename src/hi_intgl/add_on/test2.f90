
        SUBROUTINE integrate_RHO(ndim,nf,lamda,npw,n_pwr_g,src_lcl,pt_intg,coef_g,&
                &coef_h,hiresult)      
      ! changed cnr_glb_mtx to private variable shared in module
        implicit none 

        real(8),intent(in)  :: src_lcl(ndim-1),pt_intg(ndim-1)
        integer,intent(in)  :: n_pwr_g,ndim,nf,npw,lamda
        real(8),intent(out) :: hiresult(nf)

        real(8)  :: cosn(ndim),ri(ndim),gcd(ndim,ndim-1)
        real(8)  :: coef_g(0:n_pwr_g),coef_h(0:npw),coef_b(0:11,nf)

        integer :: k

        real(8) :: slop(ndim-1),rho_q,e_k,pw,nbeta

        integer :: n_pwr_k

        !=============================================

        hiresult=0.D0 ! output initialization 
        SLOP = pt_intg - src_lcl
        rho_q = norm2(SLOP)
        SLOP = SLOP/rho_q ! normalized vector, cos(theta),sin(theta)

        !!! -------------Compute n_pwr_k
        n_pwr_k = int(3+2.1214*RHO_Q)    ! NPOWF IS FROM 3 TO 9
        !order of power expansion, for parameter K in equation (3-6-62)
        if (n_pwr_k.LT.(lamda-ndim+1)) then
            n_pwr_k = int(lamda)-ndim+1
        end if       

        !!! - ----------End computing pwr_k

        call compute_coeff_B(ndim,nf,elem_type,n_pwr_g,n_pwr_k, &
                & src_glb,src_lcl,pt_intg,SLOP,RHO_Q,COEF_G,COEF_B)     

        ! Case 1 for Ek, 0 <= k <= lamda - 3        

        do k =0 ,int(lamda)-ndim

            PW= lamda - k - (ndim - 1) ! the power coefficient of rho_q
            E_k = (1.D0/RHO_Q**PW-COEF_H(INT(PW)))/(-PW)
            hiresult = hiresult + E_k*COEF_B(K,:)
!              print *,"case 1"
             !print *,"E_k = ",E_k

        end do

        ! Case 2 for Ek, k = lamda - 2
        k = int(lamda) - (ndim - 1)
        IF (k.GE.0) then
            E_k = (DLOG(rho_q) - DLOG(COEF_H(0)))
            hiresult = hiresult + COEF_B(k,:)*E_k
        end if 

        ! Case 3 for Ek,  lamda - 2 <= k <= n_pwr_k
        do k = int(int(lamda)-(ndim-1)+1),n_pwr_k

            PW= lamda - k - (ndim - 1) ! the power coefficient of rho_q

            !PW = k - lamda +(ndim - 1) ! Equ (3-6-64) Ek, power k+2 - lamda
            E_k = rho_q**(-PW)/(-PW)
            hiresult = hiresult + E_k*COEF_B(k,:)
        end do


        END SUBROUTINE 
