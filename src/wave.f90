module wave
  implicit none
      
    real(8),protected :: wk,h,beta,amp
    real(8),protected :: w1,tper,v,wl 
    real(8),protected :: timerk,rampf
    real(8),parameter :: g = 9.87
    real(8),parameter :: rho = 1.023e3
    real(8),parameter :: pi = 3.14159265359 
    !w1 => angular freq
    !v => wav num deep water
    !wk =>wave num
contains

    subroutine read_wav_data()

        implicit none

        integer :: ifwko
!        real(8) :: wl

        OPEN(1, FILE='INPUT/DATIN.txt',      STATUS='OLD') 

        READ(1,*)      IFWKO
        READ(1,*)      H, AMP, WK, BETA

        if (ifwko.eq.0) then
            if (h .le. 0.0d0) then
                w1 = dsqrt(g*wk)
            else
                w1 = dsqrt(g*wk*dtanh(wk*h))
            end if
         else
             w1 = wk
             if (h .le. 0.0d0) then
                 wk = w1**2/g
             else 
                 call waveck(w1,h,wk)
             end if
          end if

          tper = 2.*pi/w1
          beta = beta*pi/180.0d0
          v = w1**2/g
          wl = 2.0d0*pi/wk
          close(1)
    end subroutine

    subroutine waveck(sigma,h,wk)
        implicit  none

        real(8) sigma,h,wk,b,g,a,y,c
        !c  h: water depth    sigma: wave frequency
        !c  c: wave celerity  wk: wave number
        if( sigma .le. 0.0d0 )  then
            print *,' in waveck:  w=',sigma
		stop
	
        else 

            if( h .gt. 0.0d0) then
                b = g * h
                y = sigma * sigma *h/g
                a=1.0d0/(1.0d0+y*(0.66667d0+y*(0.35550d0+&
                    &       y*(0.16084d0+y*(0.063201d0+y*(0.02174d0+&
                    &       y*(0.00654d0+y*(0.00171d0+&
                    &       y*(0.00039d0+y*0.00011d0)))))))))
                c=sqrt(b/(y+a))
                wk=sigma/c
            else if( h .le. 0.0d0) then
                wk=sigma**2/g
            end if
         end if


    end subroutine
        
end module 
