!!-------------------------------------------------------------------------------
! DUTWAV
!-------------------------------------------------------------------------------
!  MODULE: motion 
!
!> @brief
!! 
!!
!! @author
!! DUT 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!

!
!------------------------------------------------------------------------------- 
module motion
    use kinds
    implicit none

    integer  nnt
    parameter (nnt=2000)    
    
    ! @var : [xc,yc,zc] body center
    real(rk) :: xc,yc,zc 
    ! @var : [xtc,ytc,ztc] rotation center
    real(rk) :: xtc,ytc,ztc 

    ! @var : [are] area
    ! @var : [xf,yf]
    ! @var : [xk2,yk2]
    ! @var : [xcf]
    real(rk) :: are,xf,yf,xk2,yk2,xcf !

    ! @var : [volm] displacement volume
    ! @var : [xb,yb,zb] buocy center
    real(rk) :: volm,xb,yb,zb

    real(rk) :: amas(6,6),bdmp(6,6)  
    real(rk) :: rmas(6,6),cdmp(6,6),crs(6,6),stkm(6,6),xia(3,3)

    real(rk) :: forcew(6),force0(6),forscd(6),ampj(6) !wave force,...
    real(rk) :: force(6),force_o(6)

    real(rk) :: trmas(6,6),visc(6,6)

    ! @var : [disp,disp_o] displacement,disp at initial time
    ! @var : [dst,dsdt_0] velocity,disp at initial time
    real(rk) ::  disp(6),disp_o(6)
    real(rk) ::  dsdt(6),dsdt_o(6)

    real(rk) ::  dsddtl(6),trxyz(3,3)   

    real(rk)  respr(6),respi(6)  

contains
    subroutine init_motion_var(nnf,nnode,nnrml)
        integer,intent(in) :: nnf,nnode,nnrml

        ! pass

    end subroutine
end module

