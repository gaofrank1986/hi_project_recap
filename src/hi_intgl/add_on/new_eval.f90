subroutine new_eval_singular(cnr_glb_mtx,ctr_glb,src_glb,src_lcl,hi_result)
      implicit none
      
      integer,parameter :: ndim = 3
      integer,parameter :: nf = 8
        real(8),intent(in) :: ctr_glb(ndim),src_glb(ndim),src_lcl(ndim-1)
        real(8),intent(inout) :: cnr_glb_mtx(ndim,nf)
        real(8),intent(out) :: hi_result(nf)
        ! nf : num of kernel funcs

        real(8) :: pt_intg(ndim-1),pt_intg_tmp(ndim - 1)  
        ! src point, integration point , temporary integration point     

        real(8) :: end_nodes(2,2),ri(3),RINT(nf)
       ! end nodes recorder, shape function
        ! rho integration result

        real(8),allocatable :: coef_g(:),coef_h(:),sf_src(:),gpl(:),gwl(:)

        !====================================================================
        integer :: id,ie,tmp,i_edge,ks,src_identifier,igl,num_converge,tmp_m(3)

        integer :: unfixed_cmp,fixed_cmp,num_edge

        real(8) :: vlc2,rho_q,comt
        real(8) :: fk,wfa,sign_value,ri_tmp,full_step_size,step_size,drdnp
        real(8) :: addup,diff_1,diff_2
        integer :: debug_flag,debug_file_id

        debug_file_id = 109
        debug_flag = 0
        !==================================

        allocate(coef_g(0:n_pwr_g),coef_h(0:npw),sf_src(elem_nd_count))
        allocate(gpl(iabs(ngl)),gwl(iabs(ngl)))


        num_edge = 2 * (num_dim - 1 ) ! 4 -----how many edges
        !print *,"============in eval_SINGULAR_ELEM=============="

        hi_result = 0.

        ri = src_glb - src_ctr_glb 

        if ( nf .eq. 8) then
                tmp_m = cnr_glb_mtx(1:3,2)
                cnr_glb_mtx(1:3,2) = cnr_glb_mtx(1:3,3)
                cnr_glb_mtx(1:3,3) = cnr_glb_mtx(1:3,5)
                cnr_glb_mtx(1:3,5) = tmp_m
                tmp_m = cnr_glb_mtx(1:3,4)
                cnr_glb_mtx(1:3,4) = cnr_glb_mtx(1:3,7)
                cnr_glb_mtx(1:3,7) = cnr_glb_mtx(1:3,6)
                cnr_glb_mtx(1:3,6) = tmp_m
        end if
        write(500,*) cnr_glb_mtx
        if (ndim == 2) then
            print *,"2d case not implemented"
        else 
            !-----------------------------------------------------------------------
            call gaussv(iabs(ngl),gpl,gwl)

            wfa=dsqrt(hi_beta*2.d0/3.d0+0.4d0)*dlog(dabs(tolgp)/2.d0)   
            fk=3.d0/8.d0*(-10.d0*dble(iabs(ngl))/wfa-1.d0)**(4.d0/3.d0)

            do i_edge = 1,num_edge ! ITERATE through each edge

                KS=KSB(i_edge)
                IF(DABS(src_lcl(IABS(KS))-DBLE(KS)/DABS(DBLE(KS))).LT.TOL) then
                    !print *,"Current edge iteration skipped! elem_id = ",this_elem_id," edge =",i_edge
                    goto 100
                end if

                do id = 1,2
                    tmp = node_grp_by_edge(3*(i_edge - 1) + ID)! determine which group of node to used
                    end_nodes(1:2,ID)=cnr_lcl_mtx(2*tmp-1 : 2*tmp) !get local node from corner table
                end do

                !==================================================
                ! for first and third edge, vertical component unchanged,1st comp change
                ! for second and forth edge, horizontal component unchanged,2rd comp change
                !       4----7----3
                !       |         |
                !       8    9    6
                !       |         |
                !       1----5----2

                unfixed_cmp=1
                if (i_edge/2*2 .eq. i_edge) unfixed_cmp=2 
                ! i_edge can is mutiple of 2
                fixed_cmp=3-unfixed_cmp
                !==================================================
                sign_value = dsign(1.d0,end_nodes(unfixed_cmp,2)-end_nodes(unfixed_cmp,1))
                !dsign(a,b) a time sign of b,end_node(:,id)
                ! sign_value is sign of second node minus first node,actually shows the direction

                vlc2=(end_nodes(fixed_cmp,2)-src_lcl(fixed_cmp))**2

                pt_intg_tmp=end_nodes(:,1) 
                !tmp integration point
                pt_intg(fixed_cmp)=pt_intg_tmp(fixed_cmp)
                ! the other xiq component is not initialised

                do num_converge=1,500 
                !control the maxium step to go thru one edge
                !also controled by step-size, example finished in less than 10 step
                !pt_intg(unfixed_cmp) is updated each time

                    diff_1 = src_lcl(unfixed_cmp)-pt_intg_tmp(unfixed_cmp)!between tmp src and end node
                    diff_2 = end_nodes(unfixed_cmp,2)-pt_intg_tmp(unfixed_cmp)!between two nodes

                    if (sign_value*diff_2 < 1.D-8) then
                        !print *, "if pt_intg_tmp coordinate out of range,exit num_converge loop"
                        !if pt_intg_tmp goes out of the edge,then stop
                        goto 100
                    end if
                    
                    ri_tmp = norm2(pt_intg_tmp-src_lcl) 
                    ! this is r from temporary integration point to local src point

                    if (0..GE.sign_value*diff_1) then
                        full_step_size = FK*ri_tmp ! fk is a factor?
                    else 
                        addup=VLc2*(1.-FK**2)+(diff_1)**2
                        IF(addup < 0.D0) THEN
                            full_step_size=sign_value*diff_1
                        ELSE
                            full_step_size=FK*(FK*sign_value*(diff_1)-DSQRT(addup))/(FK*FK-1.)
                        ENDIF
                    endif

                
                    IF (sign_value*sign_value*full_step_size+1.D-8 > sign_value*diff_2) then
                            full_step_size=sign_value*(diff_2)
                        endif
                    step_size = 0.5D0*sign_value*full_step_size
                    !step size control,step size is half the full step size,also consider sign
                   
                    !print *,"step_size = ",step_size
                    !==================================================================
                    !--- compute integral below

                    do IGL = 1,iabs(NGL) ! gaussian sampling points
                        ! this method change the double integral to line integral
                        
                        pt_intg(unfixed_cmp)=pt_intg_tmp(unfixed_cmp)+step_size*(1.D0+GPL(IGL))
                        ! update integration point position
                        
                        RHO_Q=norm2(pt_intg-src_lcl)
                        ! recaluate rho_q
                        
                        DRDNP=DABS(pt_intg(fixed_cmp)-src_lcl(fixed_cmp))/RHO_Q !sin(theta)

                               
                        call compute_coeff_GH(num_dim,num_dim - 1,npw,elem_nd_count,n_pwr_g,src_glb &
                                                & ,src_lcl,pt_intg,COEF_G,COEF_H)
                        


                        call integrate_rho(ndim,nf,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,RINT)

                        hi_result = hi_result + (dabs(step_size)*gwl(igl)*drdnp/rho_q)*rint
                        ! Equation (3-6-50)

                    end do ! igl =1,iabs(ngl)
                    
                    pt_intg_tmp(unfixed_cmp)=pt_intg_tmp(unfixed_cmp)+sign_value*full_step_size
                end do ! num_converge = 1,500


        100          end do ! i_edge = 1,num_edge
        end if


        !             write (110,*) "==========final result========="
            call swap_result(hi_result)
        !             write (110,*) result
        !         end if


    end subroutine
