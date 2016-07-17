module potential_mod
    use kinds
    implicit none
    real(rk),allocatable :: bkn(:,:),bkn_o(:,:)
    real(rk),allocatable :: unkn(:,:),unkn_o(:,:)
    real(rk),allocatable :: et(:,:),et_o(:,:)
    real(rk),allocatable :: dpdt(:,:)
    real(rk),allocatable :: dh(:,:,:),dp(:,:,:),dposi(:,:)

contains
    subroutine init_pot_var(nsys,nnf,nnode,nnrml)
        implicit none
        integer,intent(in) :: nnf,nnode,nnrml,nsys

    allocate(unkn(nnode,nsys),unkn_o(nnode,nsys))
    allocate(bkn(nnrml,nsys),bkn_o(nnrml,nsys),&
        &    et(nnf,nsys),et_o(nnf,nsys), dpdt(nnode,nsys))
    allocate(dh(4,nnf,nsys),dp(4,nnf,nsys),dposi(4,6))
    unkn=0.
    unkn_o=0.
    bkn_o=0.
    et=0.
    et_o=0.
    dh=0.
    dp=0.
    dpdt=0.
    end subroutine
    
end module