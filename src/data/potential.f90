module potential_mod
    use kinds

    implicit none

    ! @var [bkn] : surface dp/dz
    ! @var [unkn]: solved solution
    ! @var [et]  : wave height
    ! @var [dpdt]: d \phi /dt 
    ! @var [dh,dp]: dh/dt ,dp/dt
    ! @var [dposi]: difference of position at six direction

    real(rk),allocatable :: bkn(:,:),bkn_o(:,:)
    real(rk),allocatable :: unkn(:,:),unkn_o(:,:)
    real(rk),allocatable :: et(:,:),et_o(:,:)
    real(rk),allocatable :: dpdt(:,:)
    real(rk),allocatable :: dh(:,:,:),dp(:,:,:),dposi(:,:)

contains

    subroutine init_pot_var(nsys,nnf,nnode,nnrml)
        implicit none
        integer,intent(in) :: nnf,nnode,nnrml,nsys

        if (allocated(unkn)) deallocate(unkn,unkn_o)
        if (allocated(bkn)) deallocate(bkn,bkn_o,et,et_o,dpdt,dh,dp,dposi)

        allocate(unkn(nnode,nsys),unkn_o(nnode,nsys))
        allocate(bkn(nnrml,nsys),bkn_o(nnrml,nsys),&
            &    et(nnf,nsys),et_o(nnf,nsys), dpdt(nnode,nsys))

        allocate(dh(4,nnf,nsys),dp(4,nnf,nsys),dposi(0:4,6))
        !allocate(dh(4,nnf,nsys),dp(4,nnf,nsys),dposi(4,6))
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
