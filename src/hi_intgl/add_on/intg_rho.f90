
    subroutine integrate_rho(flag,ndim,nf,lamda,npw,n_pwr_g,src_lcl,pt_intg,coef_g,&
                &coef_h,hiresult)      
            ! changed cnr_glb_mtx to private variable shared in module
            implicit none 

        real(8),intent(in)  :: src_lcl(ndim-1),pt_intg(ndim-1),lamda
        integer,intent(in)  :: n_pwr_g,ndim,nf,npw,flag!flag control which kernel to use
        real(8),intent(out) :: hiresult(nf)

        real(8)  :: cosn(ndim),ri(ndim),gcd(ndim,ndim-1)
        real(8)  :: coef_g(0:n_pwr_g),coef_h(0:npw),coef_b(0:11,nf)

        integer :: k

        real(8) :: slop(ndim-1),rho_q,e_k,pw,nbeta

        integer :: n_pwr_k

        !=============================================

        hiresult=0.d0 ! output initialization 
        slop = pt_intg - src_lcl
        rho_q = norm2(slop)
        slop = slop/rho_q ! normalized vector, cos(theta),sin(theta)

        !!! -------------Compute n_pwr_k
        n_pwr_k = int(3+2.1214*RHO_Q)    ! NPOWF IS FROM 3 TO 9
        !order of power expansion, for parameter K in equation (3-6-62)
        if (n_pwr_k.LT.(lamda-ndim+1)) then
            n_pwr_k = int(lamda)-ndim+1
        end if       

        !!! - ----------End computing pwr_k

        call compute_coeff_B(flag,ndim,nf,lamda,elem_type,n_pwr_g,n_pwr_k, &
                & src_glb,src_lcl,pt_intg,coef_g,coef_b)     
        if (flag .eq. 2) then
            do k=1,8
                write(502,230) coef_b(0:11,k)
                end do
        end if
        230 format(12f10.6)
        ! Case 1 for Ek, 0 <= k <= lamda - 3        
        if (flag .eq. 1) then
            do k=1,8
                write(503,230) coef_b(0:11,k)
                end do
        end if

        do k =0 ,int(lamda)-ndim

            pw= lamda - k - (ndim - 1) ! the power coefficient of rho_q
            e_k = (1.d0/rho_q**pw-coef_h(int(pw)))/(-pw)
            hiresult = hiresult + e_k*coef_b(k,:)

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


    end subroutine 
