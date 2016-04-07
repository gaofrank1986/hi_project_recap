
   
      subroutine RIM_ELEMS()

        implicit none 

        integer  :: ie,id

        IF(num_dim == 2) n_pwr_g = (elem_type/3)*2               ! 0,  2
        IF(num_dim == 3) n_pwr_g = elem_type/2+(elem_type/9)*2        ! 2,  4,  6 
        ! refer to  Equ. 3-6-56 for parameter m
        num_intgd = 8   

        ! NDSID=2+elem_type/8 !3 !removed not used=== July 25th====

        src_ctr_glb = 0. !src center global, define the center of src for calculating r

        do ie = 1,num_elem ! can introduce parallel here!!!!!!!!!!!!!!!!!
             If (src_flag(ie) == 0) THEN    
                 !CALL ADAPTINT_ELEM(ie,src_ctr_glb,cnr_lcl_mtx,value_list(ie),GPR,GWR,GPL,GWL,INT_ELEM)
                 print *," need evaluate integral over element,src not on element"            
             else     
                 !CALL eval_SINGULAR_ELEM(ie,8,num_dim,value_list(ie,:),0)
             end if 
         end do
        !print *,value_list(1,:)
        print *,"==========final result========="
        print *,sum(value_list)
        print *,"==============================="
        print *,"end of RIM_ELEMS"
      end subroutine


