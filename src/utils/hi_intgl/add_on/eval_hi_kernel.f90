
    subroutine eval_singular_elem(passed_mtx,passed_nrml,hiresult,str_result,wak_result)
        implicit none
        integer,parameter :: nf = 8 
        integer,parameter :: ndim = 3
        real(8),intent(in) :: passed_mtx(ndim,nf), passed_nrml(ndim,nf)
        real(8),intent(out) :: hiresult(nf),str_result(nf),wak_result(nf)
        ! nf : num of kernel funcs

        real(8) :: src_lcl(ndim-1),pt_intg(ndim-1),seg_start(ndim - 1)  
        ! src point, integration point , temporary integration point     

        real(8) :: end_nodes(2,2),ri(3),RINT(nf),rint2(nf),rint3(nf)
       ! end nodes recorder, shape function
        ! rho integration result

        real(8),allocatable :: coef_g(:),coef_h(:),gpl(:),gwl(:)

        !====================================================================
        integer :: id,ie,tmp,i_edge,ks,src_identifier,igl,num_converge

        integer :: unfixed,fixed,num_edge

        real(8) :: dist_fixed2,rho_q,comt
        real(8) :: fk,wfa,edge_direct,seg_step,half_step,drdn_p
        real(8) :: sub,dist_unfixed,dist_to_end,bool_expr 
            
        integer :: debug_flag,debug_file_id

        cnr_glb_mtx = passed_mtx
        cnr_nrml=passed_nrml

        allocate(coef_g(0:n_pwr_g),coef_h(0:npw))
        allocate(gpl(iabs(ngl)),gwl(iabs(ngl)))

        num_edge = 4!2 * (ndim - 1 ) ! 4 -----how many edges
        hiresult = 0.
        str_result = 0
        wak_result = 0
            src_lcl = src_lcl_pre
            src_glb = src_glb_pre
            ri = src_glb - src_ctr_glb 

            wfa=dsqrt(hi_beta*2.d0/3.d0+0.4d0)*dlog(dabs(tolgp)/2.d0)  
            fk=3.d0/8.d0*(-10.d0*dble(iabs(ngl))/wfa-1.d0)**(4.d0/3.d0)

        if (ndim == 2) then
            print *,"2d case not implemented"
        else 
            !-----------------------------------------------------------------------
            call gaussv(iabs(ngl),gpl,gwl)!guassion_point_list, gaussian_weight_list

            do i_edge = 1,num_edge ! ITERATE through each edge

                ks=ksb(i_edge)
                if(dabs(src_lcl(iabs(ks))-dble(ks)/dabs(dble(ks))).lt.tol) then
                    !print *,"src on current edge,Current edge skipped!"
                    goto 100
                end if
                !get the two end nodes coordinates for current edge
                do id = 1,2
                    tmp = node_grp_by_edge(3*(i_edge - 1) + ID)! determine which group of node to used
                    end_nodes(1:2,ID)=cnr_lcl_mtx(2*tmp-1 : 2*tmp) !get local node from corner table
                end do
                !see on which direction ( 1 for x,2 for y) along the edge remain unchanged 
                call get_fixed_id(i_edge,fixed,unfixed)!get fixed and unfixed id

                edge_direct = dsign(1.d0,end_nodes(unfixed,2)-end_nodes(unfixed,1))
                !dsign(a,b) a time sign of b,end_node(:,id)

                !line segment start from  end node 1 for each edge
                seg_start=end_nodes(:,1) 

                pt_intg(fixed)=seg_start(fixed)
                dist_fixed2=(seg_start(fixed)-src_lcl(fixed))**2

                do num_converge=1,500 

                    dist_to_end = edge_direct*(end_nodes(unfixed,2)-seg_start(unfixed))
                    !edge end to moving point distance
                    if (dist_to_end < 1.D-8) then
                        !print *, "if seg_start located out of edge,exit num_converge loop"
                        goto 100
                    end if

                    dist_unfixed = edge_direct*(src_lcl(unfixed)-seg_start(unfixed))
                    ! compute step size====================
                    if (dist_unfixed < 1.d-8) then!if passed src unfixed 
                        seg_step = FK*norm2(seg_start-src_lcl)  ! fk is a factor?
                    else 
                        sub=(dist_fixed2+(dist_unfixed)**2)-dist_fixed2*FK**2
                        if(sub < 0.d0) then
                            seg_step=dist_unfixed
                        else
                            seg_step=fk*(fk*dist_unfixed-dsqrt(sub))/(fk**2-1.)
                        endif
                    endif

                    bool_expr = edge_direct*(edge_direct*seg_step+seg_start(unfixed))+1.d-8 

                    if (bool_expr > dist_to_end) then
                    ! if next step goes beyond the edge, use the end node
                            seg_step = dist_to_end
                    endif
                    half_step = 0.5D0*edge_direct*seg_step
                   ! =====================

                    do igl = 1,iabs(ngl) !cutting each intg segment to ngl parts 

                        pt_intg(unfixed)=seg_start(unfixed)+half_step*(1.d0+gpl(igl))
                        ! 1+gpl , is shifted gpl to start at least from 0
                        ! gpl span from -1 to 1, so intg-seg-size if half buffer point step
                        rho_q=norm2(pt_intg-src_lcl)
                        drdn_p=dabs(pt_intg(fixed)-src_lcl(fixed))/rho_q !sin(theta)

                        call compute_coeff_gh(num_dim,num_dim - 1,npw,elem_type,n_pwr_g,src_glb &
                                                & ,src_lcl,pt_intg,coef_g,coef_h)

                        call integrate_rho(1,ndim,nf,hi_beta,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,rint)
                        call integrate_rho(2,ndim,nf,2.d0,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,rint2)
                        call integrate_rho(3,ndim,nf,1.d0,npw,n_pwr_g,src_lcl,pt_intg,coef_g,coef_h,rint3)

                        hiresult = hiresult + (dabs(half_step)*gwl(igl)*drdn_p/rho_q)*rint
                        str_result = str_result + (dabs(half_step)*gwl(igl)*drdn_p/rho_q)*rint2
                        wak_result = wak_result + (dabs(half_step)*gwl(igl)*drdn_p/rho_q)*rint3
                        ! equation (3-6-50)
                    end do ! igl =1,iabs(ngl)

                    seg_start(unfixed)=seg_start(unfixed)+edge_direct*seg_step
                end do ! num_converge = 1,500


        100          end do ! i_edge = 1,num_edge
        end if

        call swap_result(hiresult)
        call swap_result(str_result)
        call swap_result(wak_result)


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
