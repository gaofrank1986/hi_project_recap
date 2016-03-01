module gradient
    use mfunc_mod,only: inverse
    implicit none
    real(8),parameter :: param_pos(8,2) = &
        reshape((/-1. ,0.  ,1.  ,1.  ,1. ,0.  ,-1. ,-1. ,&
                  -1. ,-1. ,-1. ,0.  ,1. ,1.  ,1.  ,0./),(/8,2/))
    real(8),allocatable,save :: sf_info(:,:,:),jacob_info(:,:,:)
    
contains
        
subroutine init_gradient(nnf,nelemf,xyze,nodele,nodelj,debug)
    integer,intent(in) :: nnf,nelemf,nodele(nnf),nodelj(nnf)
    logical,intent(in),optional ::debug
    real(8),intent(in) :: xyze(2,8,nnf)
    real(8) :: global_pos_value(2,8),jacobian_mt(2,2)
    real(8):: sf(8),dsf(2,8),tmp(2,2)

    integer :: i
    allocate(sf_info(nnf,3,8),jacob_info(nnf,2,2))
    do i = 1,nnf
        global_pos_value = xyze(1:2,1:8,nodele(i))
        if (present(debug) .and. debug) then
        print *,"global pos"
        !print *,global_pos_value
        endif
        call spfunc8(param_pos(nodelj(i),:),sf,dsf)
        sf_info(i,1,:) = sf
        sf_info(i,2:3,:) = dsf
        if (present(debug) .and. debug) then
        print *,"dsf"
        !print *,sf_info(1,2,:)
        !print *,sf_info(1,3,:)
        end if
        jacobian_mt = matmul(dsf,transpose(global_pos_value))!(global over local jacobain) 
        if (present(debug) .and. debug) then
        print *,"jacobian before inversion"
        print *,jacobian_mt
        end if
        !inverse jacob
        call inverse(jacobian_mt,tmp,2,2)!FIXME inverse will modify jacobian_mt,and return a tranposed inverse
        jacob_info(i,:,:) = transpose(tmp)
!        print *,i
        !write (*,100) jacob_info(i,:,:)
        !100 format(4f10.6)
        end do


end subroutine
subroutine eval_gradient(inode,surface_value,global_grad)
    implicit none
    integer,intent(in) :: inode
    real(8),intent(out) :: global_grad(2,1)
    real(8),dimension(8,1),intent(in) :: surface_value
    real(8) :: local_grad(2,1),tmp(2,8)
    tmp = sf_info(inode,2:3,:)

    
    local_grad = matmul(tmp,surface_value)
    global_grad  = matmul(jacob_info(inode,:,:),local_grad)
    ! jacobian inverse(identically local over global jacobian) * local grad
    
end subroutine 

subroutine test_gradient()
!call  init_gradient(1,1,2,transpose(param_pos),(/1/),(/2/),.true.)
  print *,"jacob after inversion"
  !print *,jacob_info
end subroutine



        subroutine spfunc8(pos,sf,dsf)
        implicit  none
        real(8),intent(in) :: pos(2)

	  !REAL*8,INTENT(IN) :: SI,ETA
	  real*8,intent(out):: sf(8),dsf(2,8)
        real(8) :: si,eta
        si = pos(1)
        eta = pos(2)
        SF(1)=-0.25D0*(1.0D0-SI)*(1.0D0-ETA)*(1.0D0+SI+ETA)
        SF(2)= 0.5D0 *(1.0D0-SI*SI)*(1.0-ETA)
        SF(3)= 0.25D0*(1.0D0+SI)*(1.0D0-ETA)*(SI-ETA-1.0D0)
        SF(4)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0D0+SI)
        SF(5)= 0.25D0*(1.0D0+SI)*(1.0D0+ETA)*(SI+ETA-1.0D0)
        SF(6)= 0.5D0 *(1.0-SI*SI)*(1.0D0+ETA)
        SF(7)=-0.25D0*(1.0D0-SI)*(1.0D0+ETA)*(SI-ETA+1.0D0)
        SF(8)= 0.5D0 *(1.0D0-ETA*ETA)*(1.0-SI)            
      
      
        DSF(1,1)= 0.25D0*(2.0D0*SI+ETA)*(1.0D0-ETA)
        DSF(1,2)=-SI*(1.0D0-ETA)     	
        DSF(1,3)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0-ETA)
        DSF(1,4)= 0.5D0*(1.0D0-ETA*ETA)
        DSF(1,5)= 0.25D0*(2.0*SI+ETA)*(1.0D0+ETA)
        DSF(1,6)=-SI*(1.0D0+ETA)     
        DSF(1,7)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0+ETA)
        DSF(1,8)=-0.5D0*(1.0-ETA*ETA)

        !DSF(1,1)= 0.25D0*(2.0D0*SI+ETA)*(1.0D0-ETA)
        !DSF(1,2)=-SI*(1.0D0-ETA)     	
        !DSF(1,3)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0-ETA)
        !DSF(1,4)= 0.5D0*(1.0D0-ETA*ETA)
        !DSF(1,5)= 0.25D0*(2.0*SI+ETA)*(1.0D0+ETA)
        !DSF(1,6)=-SI*(1.0D0+ETA)     
        !DSF(1,7)= 0.25D0*(2.0D0*SI-ETA)*(1.0D0+ETA)
        !DSF(1,8)=-0.5D0*(1.0-ETA*ETA)

        DSF(2,1)= 0.25D0*(SI+2.0D0*ETA)*(1.0D0-SI)
        DSF(2,2)=-0.5D0*(1.0-SI*SI)
        DSF(2,3)= 0.25D0*(1.0D0+SI)*(2.0D0*ETA-SI)
        DSF(2,4)=-(1.0D0+SI)*ETA
        DSF(2,5)= 0.25*(1.0D0+SI)*(SI+2.0D0*ETA)
        DSF(2,6)= 0.5D0*(1.0D0-SI*SI)
        DSF(2,7)=-0.25D0*(1.0D0-SI)*(SI-2.0D0*ETA)
        DSF(2,8)=-(1.0D0-SI)*ETA   

	RETURN
	END SUBROUTINE SPFUNC8
        end module
