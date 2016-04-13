module motion
    use kinds
    implicit none
    real(8) :: xc,yc,zc ! body center

    integer  nnt
    parameter (nnt=2000)    

    real*8 xtc,ytc,ztc ! rotation center
    real*8 are,xf,yf,xk2,yk2,xcf !
    real*8 volm,xb,yb,zb!displacement volume, buocy center

    real*8 amas(6,6),bdmp(6,6)  
    real*8 rmas(6,6),cdmp(6,6),crs(6,6),stkm(6,6),xia(3,3)

    real*8 forcew(6),force0(6),forscd(6),ampj(6) !wave force,...
    real*8 force(6),force_o(6)

    real*8 trmas(6,6),visc(6,6)

    real*8  disp(6),dsdt(6),disp_o(6)!displacement,velocity,disp at initial time
    real(8) ::  dsdt_o(6),dsddtl(6),trxyz(3,3)   

    real*8  respr(6),respi(6)  
    !      xc = 1
    !      yc = 2
    !      zc = 3
contains
    !subroutine init_motion_var(nnf,nnode,nnrml)
        !integer,intent(in) :: nnf,nnode,nnrml
    !end subroutine



end module

