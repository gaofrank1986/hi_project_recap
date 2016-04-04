!-------------------------------------------------------------------------------
! DUTWAV
!-------------------------------------------------------------------------------
!  MODULE: shape_funcs
!
!> @brief
!! <shape function>
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
! 07 Mar 2015 - Added gcombo to wrap mirrored source or sink. 
!
!-------------------------------------------------------------------------------
module shape_funcs

contains
    !c
    !c *********************************************************
    !c *                                                       *
    !c * calculate the shape functions and their derivatives   *
    !c * for 6 node triangular elements                        *
    !c *                                                       *
    !c *********************************************************     
    !c
    subroutine spfunc6(si,eta,sf,dsf)
        implicit  none
        
        real*8,intent(in) :: si,eta
        real*8,intent(out):: sf(8),dsf(2,8)
        
        sf(1)=(1.0d0-si-eta)*(1.0d0-2.0d0*si-2.0d0*eta)
        sf(2)=si *(2.0d0*si -1.0d0)
        sf(3)=eta*(2.0d0*eta-1.0d0)
        sf(4)=4.0d0*si*(1.0d0-si-eta)
        sf(5)=4.0d0*si*eta
        sf(6)=4.0d0*eta*(1.0d0-si-eta)

        dsf(1,1)=4.0d0*si+4.0*eta-3.0d0
        dsf(1,2)=4.0d0*si-1.0d0
        dsf(1,3)=0.0d0
        dsf(1,4)=4.0d0-8.0d0*si-4.0d0*eta
        dsf(1,5)=4.0d0*eta
        dsf(1,6)=-4.0*eta

        dsf(2,1)=4.0d0*si+4.0d0*eta-3.0d0
        dsf(2,2)=0.0d0
        dsf(2,3)=4.0d0*eta-1.0d0
        dsf(2,4)=-4.0*si
        dsf(2,5)=4.0*si
        dsf(2,6)=4.0-8.0d0*eta-4.0*si
        
    end subroutine spfunc6
    
    !c *********************************************************
    !c *                                                       *
    !c * calculate the shape functions and their derivatives   *
    !c * for 8 node quadrilaterial elements                    *
    !c *                                                       *
    !c *********************************************************     
    
    subroutine spfunc8(si,eta,sf,dsf)
        implicit  none
        
        real*8,intent(in) :: si,eta
        real*8,intent(out):: sf(8),dsf(2,8)
        
        sf(1)=-0.25d0*(1.0d0-si)*(1.0d0-eta)*(1.0d0+si+eta)
        sf(2)= 0.5d0 *(1.0d0-si*si)*(1.0-eta)
        sf(3)= 0.25d0*(1.0d0+si)*(1.0d0-eta)*(si-eta-1.0d0)
        sf(4)= 0.5d0 *(1.0d0-eta*eta)*(1.0d0+si)
        sf(5)= 0.25d0*(1.0d0+si)*(1.0d0+eta)*(si+eta-1.0d0)
        sf(6)= 0.5d0 *(1.0-si*si)*(1.0d0+eta)
        sf(7)=-0.25d0*(1.0d0-si)*(1.0d0+eta)*(si-eta+1.0d0)
        sf(8)= 0.5d0 *(1.0d0-eta*eta)*(1.0-si)            
        
        dsf(1,1)= 0.25d0*(2.0d0*si+eta)*(1.0d0-eta)
        dsf(1,2)=-si*(1.0d0-eta)     	
        dsf(1,3)= 0.25d0*(2.0d0*si-eta)*(1.0d0-eta)
        dsf(1,4)= 0.5d0*(1.0d0-eta*eta)
        dsf(1,5)= 0.25d0*(2.0*si+eta)*(1.0d0+eta)
        dsf(1,6)=-si*(1.0d0+eta)     
        dsf(1,7)= 0.25d0*(2.0d0*si-eta)*(1.0d0+eta)
        dsf(1,8)=-0.5d0*(1.0-eta*eta)
        
        dsf(2,1)= 0.25d0*(si+2.0d0*eta)*(1.0d0-si)
        dsf(2,2)=-0.5d0*(1.0-si*si)
        dsf(2,3)= 0.25d0*(1.0d0+si)*(2.0d0*eta-si)
        dsf(2,4)=-(1.0d0+si)*eta
        dsf(2,5)= 0.25*(1.0d0+si)*(si+2.0d0*eta)
        dsf(2,6)= 0.5d0*(1.0d0-si*si)
        dsf(2,7)=-0.25d0*(1.0d0-si)*(si-2.0d0*eta)
        dsf(2,8)=-(1.0d0-si)*eta   
        
    end subroutine spfunc8
    
    
    !c *********************************************************
    !c *                                                       *
    !c * calculate the shape functions and their derivatives   *
    !c * for 6 node triangular elements                        *
    !c *                                                       *
    !c *********************************************************     

    subroutine spfunc6_1(si,eta,sf,dsf,ddsf)
        implicit  none
        
        real*8,intent(in) :: si,eta
        real*8,intent(out):: sf(8),dsf(2,8),ddsf(3,8)

        sf(1)=(1.0d0-si-eta)*(1.0d0-2.0d0*si-2.0d0*eta)
        sf(2)=si *(2.0d0*si -1.0d0)
        sf(3)=eta*(2.0d0*eta-1.0d0)
        sf(4)=4.0d0*si*(1.0d0-si-eta)
        sf(5)=4.0d0*si*eta
        sf(6)=4.0d0*eta*(1.0d0-si-eta)
        
        dsf(1,1)=4.0d0*si+4.0*eta-3.0d0
        dsf(1,2)=4.0d0*si-1.0d0
        dsf(1,3)=0.0d0
        dsf(1,4)=4.0d0-8.0d0*si-4.0d0*eta
        dsf(1,5)=4.0d0*eta
        dsf(1,6)=-4.0*eta

        dsf(2,1)=4.0d0*si+4.0d0*eta-3.0d0
        dsf(2,2)=0.0d0
        dsf(2,3)=4.0d0*eta-1.0d0
        dsf(2,4)=-4.0*si
        dsf(2,5)=4.0*si
        dsf(2,6)=4.0-8.0d0*eta-4.0*si
        
        
        ddsf(1,1)= 4.d0
        ddsf(1,2)= 4.0d0     	
        ddsf(1,3)= 0.0d0
        ddsf(1,4)=-8.0d0
        ddsf(1,5)= 0.0d0
        ddsf(1,6)= 0.0d0     
        
        ddsf(2,1)= 4.0d0
        ddsf(2,2)= 0.0d0
        ddsf(2,3)= 4.0d0
        ddsf(2,4)= 0.0d0
        ddsf(2,5)= 0.0d0
        ddsf(2,6)=-8.0d0
        
        ddsf(3,1)= 4.0d0
        ddsf(3,2)= 0.0d0   	
        ddsf(3,3)= 0.0d0
        ddsf(3,4)=-4.0d0
        ddsf(3,5)= 4.0d0
        ddsf(3,6)=-4.0d0     
        

    end subroutine spfunc6_1
    
    
    
    subroutine spfunc8_1(si,eta,sf,dsf,ddsf)
        implicit  none
        
        real*8,intent(in) :: si,eta
        real*8,intent(out):: sf(8),dsf(2,8),ddsf(3,8)

        sf(1)=-0.25d0*(1.0d0-si)*(1.0d0-eta)*(1.0d0+si+eta)
        sf(2)= 0.5d0 *(1.0d0-si*si)*(1.0-eta)
        sf(3)= 0.25d0*(1.0d0+si)*(1.0d0-eta)*(si-eta-1.0d0)
        sf(4)= 0.5d0 *(1.0d0-eta*eta)*(1.0d0+si)
        sf(5)= 0.25d0*(1.0d0+si)*(1.0d0+eta)*(si+eta-1.0d0)
        sf(6)= 0.5d0 *(1.0-si*si)*(1.0d0+eta)
        sf(7)=-0.25d0*(1.0d0-si)*(1.0d0+eta)*(si-eta+1.0d0)
        sf(8)= 0.5d0 *(1.0d0-eta*eta)*(1.0-si)            

        
        dsf(1,1)= 0.25d0*(2.0d0*si+eta)*(1.0d0-eta)
        dsf(1,2)=-si*(1.0d0-eta)     	
        dsf(1,3)= 0.25d0*(2.0d0*si-eta)*(1.0d0-eta)
        dsf(1,4)= 0.5d0*(1.0d0-eta*eta)
        dsf(1,5)= 0.25d0*(2.0*si+eta)*(1.0d0+eta)
        dsf(1,6)=-si*(1.0d0+eta)     
        dsf(1,7)= 0.25d0*(2.0d0*si-eta)*(1.0d0+eta)
        dsf(1,8)=-0.5d0*(1.0-eta*eta)

        dsf(2,1)= 0.25d0*(si+2.0d0*eta)*(1.0d0-si)
        dsf(2,2)=-0.5d0*(1.0-si*si)
        dsf(2,3)= 0.25d0*(1.0d0+si)*(2.0d0*eta-si)
        dsf(2,4)=-(1.0d0+si)*eta
        dsf(2,5)= 0.25*(1.0d0+si)*(si+2.0d0*eta)
        dsf(2,6)= 0.5d0*(1.0d0-si*si)
        dsf(2,7)=-0.25d0*(1.0d0-si)*(si-2.0d0*eta)
        dsf(2,8)=-(1.0d0-si)*eta   
        
        
        ddsf(1,1)=0.5d0*(1.0d0-eta)
        ddsf(1,2)=-(1.0d0-eta)     	
        ddsf(1,3)= 0.5d0*(1.0d0-eta)
        ddsf(1,4)= 0.0d0
        ddsf(1,5)= 0.5d0*(1.0d0+eta)
        ddsf(1,6)=-(1.0d0+eta)     
        ddsf(1,7)= 0.5d0*(1.0d0+eta)
        ddsf(1,8)= 0.0d0
        
        ddsf(2,1)= 0.50d0*(1.0d0-si)
        ddsf(2,2)= 0.0d0
        ddsf(2,3)= 0.50d0*(1.0d0+si)
        ddsf(2,4)=-(1.0d0+si)
        ddsf(2,5)= 0.50d0*(1.0d0+si)
        ddsf(2,6)= 0.0d0
        ddsf(2,7)= 0.50d0*(1.0d0-si)
        ddsf(2,8)=-(1.0d0-si) 
        
        ddsf(3,1)= 0.25d0*(-2.0d0*si-eta*2.0d0+1.0d0)
        ddsf(3,2)= si     	
        ddsf(3,3)= 0.25d0*(-2.0d0*si+eta*2.0d0-1.0d0)
        ddsf(3,4)= -eta
        ddsf(3,5)= 0.25d0*(2.0d0*si+eta*2.0d0+1.0d0)
        ddsf(3,6)=-si     
        ddsf(3,7)= 0.25d0*(2.0d0*si-eta*2.0d0-1.0d0)
        ddsf(3,8)= eta
        

    end subroutine spfunc8_1
end module
