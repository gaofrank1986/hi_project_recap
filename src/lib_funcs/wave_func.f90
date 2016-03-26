module wave_funcs_simple
    use kinds
    use wave,only:amp,wk,beta,w1,g,h
    use wave,only:timerk,rampf

contains
    
    !< Wave profile function
    function eti(x,y)
        implicit  none
        real(rk),intent(in):: x,y
        real(rk) ::   eti
        eti=amp*dcos(wk*(x*dcos(beta)+y*dsin(beta)) -w1*timerk)
        eti=rampf*eti

        return
    end function eti

    !< temporal Derivate of wave elevation
    function deti(x,y)
        implicit  none

        real(rk),intent(in):: x,y
        real(rk)  :: deti

        deti=w1*amp*dsin(wk*(x*dcos(beta)+y*dsin(beta)) -w1*timerk)
        deti=rampf*deti

    end function deti

    ! < incident wave potential at a spatial point,at timerk
    function poxy(x,y,z)
        implicit  none

        real(rk),intent(in):: x,y,z
        real(rk)  poxy
        real(rk)  wkx,dum

        !< shallow water case
        if(h.gt.0.0d0) then
            dum=amp*g/w1*dcosh(wk*(z+h))/dcosh(wk*h)
        else
        !<  deep water case
            dum=amp*g/w1*dexp(wk*z)
        endif

        wkx=wk*(x*dcos(beta)+y*dsin(beta))
        poxy=dum*dsin(wkx-w1*timerk)

        poxy = rampf*poxy

    end function poxy
    
    ! < time derivative of inicidnet wave potential at given spartial point at timerk
    function dpot(x,y,z)

        implicit  none

        real(rk),intent(in):: x,y,z
        real(rk)  dpot
        real(rk)  wkx,dum

        if(h.gt.0.0d0) then
            dum=-amp*g*dcosh(wk*(z+h))/dcosh(wk*h)
        else
            dum=-amp*g*dexp(wk*z)
        endif

        wkx=wk*( x*dcos(beta)+y*dsin(beta))
        dpot=dum*dcos(wkx-w1*timerk)

        dpot = rampf*dpot

    end function dpot


    !< spacial derivatives of incident wave potential
    
    subroutine  dinp(x,y,z,dpox,dpoy,dpoz)

        implicit    none

        real(rk),intent(in)::   x,y,z
        real(rk),intent(out)::  dpox,dpoy,dpoz

        real*8 dum,wkx

        dum=amp*g/w1
        wkx=wk*(x*dcos(beta)+y*dsin(beta))

        if (h.gt.0.0d0) then           
            dpox= dum*wk*dcos(beta)* &
               dcosh(wk*(z+h))/dcosh(wk*h)*dcos(wkx-w1*timerk)
            dpoy= dum*wk*dsin(beta)* &
               dcosh(wk*(z+h))/dcosh(wk*h)*dcos(wkx-w1*timerk)
            dpoz= dum*wk*&
               dsinh(wk*(z+h))/dcosh(wk*h)*dsin(wkx-w1*timerk)

        else

            dpox= dum*wk*dcos(beta)*dexp(wk*z)*dcos(wkx-w1*timerk)
            dpoy= dum*wk*dsin(beta)*dexp(wk*z)*dcos(wkx-w1*timerk)
            dpoz= dum*wk*dexp(wk*z)*dsin(wkx-w1*timerk)

        endif

        dpox=rampf*dpox
        dpoy=rampf*dpoy
        dpoz=rampf*dpoz


    end subroutine  dinp 

    !< Constructed potential function for solving fterm coefficients
    subroutine  dinp1(i,p,phi,dpdx)
        implicit    none
        integer,intent(in) :: i
        real(8),intent(in)::  p(3)
        real*8,intent(out)::  phi,dpdx(3)


        if (i==1) then
            phi=1.0d0           
            dpdx=0.0d0
        else if (i==2) then
            phi=p(1)       
            dpdx=(/1.0d0,0.0d0,0.0d0/)
        else if (i==3) then
            phi=p(2)           
            dpdx=(/0.0d0,1.0d0,0.0d0/)
        else if (i==4) then
            phi=p(3)           
            dpdx=(/0.0d0,0.0d0,1.0d0/)
        endif
    end subroutine

end module

