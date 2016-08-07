!-------------------------------------------------------------------------------
! DUTWAV
!-------------------------------------------------------------------------------
!  MODULE: free term
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Song Gao,DUT 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Mar 2015 - 
!
!-------------------------------------------------------------------------------
module gradient
    use kinds
    use linalg,only: inverse
    use shape_funcs,only:spfunc8
    implicit none
    ! @var : [param_pos] local node position of a rectangle
    real(rk),parameter :: param_pos(8,2) = &
        reshape((/-1. , 0. , 1. ,1. ,1. ,0. ,-1. ,-1. ,&
                  -1. ,-1. ,-1. ,0. ,1. ,1. , 1. , 0. /),(/8,2/))
    ! @var : [sf_infp
    real(rk),allocatable,save :: sf_info(:,:,:),jacob_info(:,:,:)

contains

    ! @param nodele:connected elem
    ! @param nodelj:node is jth node in connected elem
    ! @param : [xyze(2,8,nnf)] x,y info for each node of each elem
    ! @param : [nodele] first connected elem
    ! @param : [nodej] first connected node position in elem

    subroutine init_gradient(nnf,nelemf,xyze,nodele,nodelj,debug)
        integer,intent(in) :: nnf,nelemf,nodele(nnf),nodelj(nnf)
        logical,intent(in),optional ::debug
        real(rk),intent(in) :: xyze(2,8,nnf)
        real(rk) :: p_xy(2,8),jacobian_mt(2,2)
        real(rk):: sf(8),dsf(2,8),tmp(2,2),xi(2)

        integer :: i
        allocate(sf_info(nnf,3,8),jacob_info(nnf,2,2))
        do i = 1,nnf
            ! get local pos
            xi= param_pos(nodelj(i),:)
            ! get global pos
            p_xy = xyze(1:2,1:8,nodele(i))
            ! get shape func and first deriv of shape func
            call spfunc8(xi(1),xi(2),sf,dsf)
            !nodelj is the pos of node i in the element

            !record the info in sf_info
            sf_info(i,1,:) = sf
            sf_info(i,2:3,:) = dsf

            !cal the jacobian mtx at node i
            jacobian_mt = matmul(dsf,transpose(p_xy))!(global over local jacobain) 

            !inverse jacob
            call inverse(jacobian_mt,tmp,2,2)!FIXME inverse will modify jacobian_mt,and return a tranposed inverse
            jacob_info(i,:,:) = transpose(tmp)

        end do


    end subroutine

    subroutine eval_gradient(inode,surface_value,global_grad)
        implicit none
        integer,intent(in) :: inode
        real(rk),intent(out) :: global_grad(2,1)
        real(rk),dimension(8,1),intent(in) :: surface_value
        real(rk) :: local_grad(2,1),tmp(2,8)

        tmp = sf_info(inode,2:3,:)
        !dsf info of node i
        ! local gradient = dsf.*z_value 
        local_grad = matmul(tmp,surface_value)
        global_grad  = matmul(jacob_info(inode,:,:),local_grad)
        ! jacobian inverse(identically local over global jacobian) * local grad

    end subroutine 

    subroutine test_gradient()
        !call  init_gradient(1,1,2,transpose(param_pos),(/1/),(/2/),.true.)
        print *,"jacob after inversion"
        !print *,jacob_info
    end subroutine




end module
