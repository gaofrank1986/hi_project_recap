module data_all
    use mesh
    use wave
    use time_mod
    use potential_mod
    use motion
    use matrix_mod
    use misc_var
    implicit none
contains
    subroutine init_data(nsys,nnode,nnrml,nnf)
        implicit none
        integer,intent(in) :: nsys,nnode,nnrml,nnf
        call init_matrix_var(nsys,nnode,nnrml)
        call init_pot_var(nsys,nnf,nnode,nnrml)
        call init_misc_var(nsys,nnode,nnrml)
    end subroutine
end module

