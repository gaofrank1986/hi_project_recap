
    subroutine eval_singular_elem(passed_mtx,hiresult)
        implicit none
        integer,parameter :: nf = 8 
        integer,parameter :: ndim = 3
        real(8),intent(in) :: passed_mtx(3,8)
        real(8),intent(out) :: hiresult(nf)
        ! nf : num of kernel funcs

        real(8) :: src_lcl(ndim-1),pt_intg(ndim-1),buf_pnt(ndim - 1)  
        ! src point, integration point , temporary integration point     

        real(8) :: end_nodes(2,2),ri(3),RINT(nf)
       ! end nodes recorder, shape function
        ! rho integration result

        real(8),allocatable :: coef_g(:),coef_h(:),gpl(:),gwl(:)

        !====================================================================
        integer :: id,ie,tmp,i_edge,ks,src_identifier,igl,num_converge

        integer :: unfixed,fixed,num_edge

        real(8) :: vlc2,rho_q,comt
        real(8) :: fk,wfa,edge_direct,ri_tmp,buf_step,intg_pnt_step,drdn_p
        real(8) :: addup,diff_1,diff_2,bool_expr 
            
        integer :: debug_flag,debug_file_id

        cnr_glb_mtx = passed_mtx
        debug_file_id = 109
        debug_flag = 0
        !==================================

        allocate(coef_g(0:n_pwr_g),coef_h(0:npw))
        allocate(gpl(iabs(ngl)),gwl(iabs(ngl)))

        num_edge = 2 * (ndim - 1 ) ! 4 -----how many edges
        hiresult = 0.

            src_lcl = src_lcl_preset
            src_glb = src_glb_preset
            ri = src_glb - src_ctr_glb 


        if (ndim == 2) then
            print *,"2d case not implemented"
        else 
            !-----------------------------------------------------------------------
            call gaussv(iabs(ngl),gpl,gwl)!guassion_point_list, gaussian_weight_list

            wfa=dsqrt(hi_beta*2.d0/3.d0+0.4d0)*dlog(dabs(tolgp)/2.d0)  
            fk=3.d0/8.d0*(-10.d0*dble(iabs(ngl))/wfa-1.d0)**(4.d0/3.d0)

            do i_edge = 1,num_edge ! ITERATE through each edge

                ks=ksb(i_edge)
                if(dabs(src_lcl(iabs(ks))-dble(ks)/dabs(dble(ks))).lt.tol) then
                    !print *,"Current edge iteration skipped! elem_id = ",this_elem_id," edge =",i_edge
                    goto 100
                end if
                !find end node local position
                do id = 1,2
                    tmp = node_grp_by_edge(3*(i_edge - 1) + ID)! determine which group of node to used
                    end_nodes(1:2,ID)=cnr_lcl_mtx(2*tmp-1 : 2*tmp) !get local node from corner table
                end do

                call get_fixed_id(i_edge,fixed,unfixed)!get fixed and unfixed id

                edge_direct = dsign(1.d0,end_nodes(unfixed,2)-end_nodes(unfixed,1))
                !dsign(a,b) a time sign of b,end_node(:,id)

                VLc2=(end_nodes(fixed,2)-src_lcl(fixed))**2

                !move point goes along the four edges
                !it start from  end node 1 for each edge
                buf_pnt=end_nodes(:,1) 
                pt_intg(fixed)=buf_pnt(fixed)

                do num_converge=1,500 
                !control the maxium step to go thru one edge
                !also controled by step-size, example finished in less than 10 step
                !pt_intg(unfixed) is updated each time

                    diff_1 = edge_direct*(src_lcl(unfixed)-buf_pnt(unfixed))
                    !src point to moving pt distance
                    diff_2 = edge_direct*(end_nodes(unfixed,2)-buf_pnt(unfixed))
                    !edge end to moving point distance

                    if (diff_2 < 1.D-8) then
                        !print *, "if buf_pnt located out of edge,exit num_converge loop"
                        !if buf_pnt goes out of the edge,then stop
                        goto 100
                    end if

                    ri_tmp = norm2(buf_pnt-src_lcl) 
                    ! this is r from moving  point to local src point
                        
                    ! compute step size====================
                    if (diff_1 < 1.d-8) then
                        buf_step = FK*ri_tmp ! fk is a factor?
                    else 
                        addup=VLc2*(1.-FK**2)+(diff_1)**2
                        IF(addup < 0.D0) THEN
                            buf_step=diff_1
                        ELSE
                            buf_step=FK*(FK*diff_1-DSQRT(addup))/(FK*FK-1.)
                        ENDIF
                    endif

                    bool_expr = edge_direct*(edge_direct*buf_step+buf_pnt(unfixed))+1.d-8 

                    if (bool_expr > diff_2) then
                            buf_step = diff_2
                    endif
                    intg_pnt_step = 0.5D0*edge_direct*buf_step
                   ! step size finished ======================

                    do igl = 1,iabs(ngl) ! gaussian sampling points

                        pt_intg(unfixed)=buf_pnt(unfixed)+intg_pnt_step*(1.d0+gpl(igl))
                        rho_q=norm2(pt_intg-src_lcl)
                        drdn_p=dabs(pt_intg(fixed)-src_lcl(fixed))/rho_q !sin(theta)

                        call compute_coeff_gh(num_dim,num_dim - 1,npw,elem_type,n_pwr_g,src_glb &
                                                & ,src_lcl,pt_intg,coef_g,coef_h)

                        call integrate_rho(ndim,nf,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,rint)

                        hiresult = hiresult + (dabs(intg_pnt_step)*gwl(igl)*drdn_p/rho_q)*rint
                        ! equation (3-6-50)
                    end do ! igl =1,iabs(ngl)

                    buf_pnt(unfixed)=buf_pnt(unfixed)+edge_direct*buf_step
                end do ! num_converge = 1,500


        100          end do ! i_edge = 1,num_edge
        end if

        call swap_result(hiresult)


    end subroutine


    subroutine get_fixed_id(i_edge,fixed_id,unfixed_id)
        implicit none

        integer,intent(in) :: i_edge
        integer,intent(out) :: fixed_id,unfixed_id
        ! id can be either 1 or 2,indicating which coordinate  is no changed
        ! along the edge
        !==================================================
        ! for first and third edge, vertical component unchanged,1st comp change
        ! for second and forth edge, horizontal component unchanged,2rd comp change
        !       4----7----3
        !       |         |
        !       8    9    6
        !       |         |
        !       1----5----2

        !==================================================
        unfixed_id=1!x direction changed, horizontal case
        if (i_edge/2*2 .eq. i_edge) then
        !check if i_edge is even
            unfixed_id=2 ! y direction changed,vertical case
        end if
        ! i_edge can is mutiple of 2
        ! edge 1 and 3 is horizontal, while 2 and 4 are vertical
        fixed_id=3-unfixed_id
    end subroutine
