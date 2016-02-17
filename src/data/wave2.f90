module wave
    implicit none
      real(8),protected :: wf_temp(200)
      real(8),protected :: wk,h,beta,amp
      real(8),protected :: w1,tper,v,wl
      real(8),public :: timerk,rampf
      real(8)  AMPN(200),Freq(200),WKN(200),Phi_w(200),TP_wv(200) 
      real(8) :: g,rho,pi4
      integer,parameter :: nfreq = 200
      real(8),parameter :: pi = 3.14159265359 
      integer :: iorder,nwave
      real(8) :: time,tstep


contains
      subroutine read_wav_data()
              implicit none
              integer :: i,m,kp,ifwko,wvsimu,ntnum 
              real(8) :: ntime
        OPEN(1, FILE='INPUT/DATIN.txt',      STATUS='OLD') 
       READ(1,*)      Rho, G  
       READ(1,*)      H,  BETA
       READ(1,*)    IORDER
       READ(1,*)    WVSimu, NTNUM 
       READ(1,*)    IFWKO, Nwave
!
       DO I=1, Nwave
       READ(1,*) M,  WF_Temp(I), AMPN(I), PHi_w(I)     
       ENDDO     
!
! WVSimu: number of waves to simulate, 
! NtNUM: number of time steps in a wave
! Nwave: number of compositing waves for irregular incident waves
!

        IF (IFWKO .EQ. 0)  THEN
                WKN(:)=WF_Temp(:)

                DO I=1, Nwave
                IF(H .LE. 0.0d0) THEN
                        FREQ(I)=SQRT(G*WKN(I))
                ELSE
                        FREQ(I)=SQRT(G*WKN(I)*TANH(WKN(I)*H))
                END IF
                TP_wv(I)=2.0*PI/FREQ(I)

                END DO

! --------------------------------
 
        ELSE IF (IFWKO .EQ. 1)  THEN

                FREQ(:)=WF_Temp(:)

                DO I=1, NWave
                IF(H .LE. 0.0D0) THEN
                        WKN(I)=FREQ(I)*FREQ(I)/G
                   ELSE     
                           CALL WAVECK(FREQ(I),H,G,WKN(I))               ! COMPUTE WAVE NUMBER
                   END IF

                   TP_wv(I)=2.0*PI/FREQ(I)

                   END DO

           ENDIF
              
! =============================================================
       
!
!  Waves: number of waves to simulate (may be 10 or 0.5)
!  NTnum: number of time steps in one wave
!  NPLOUT=1-5, Output the wave profile with step intervals 
!
          PI4=4.0*PI
        
        WRITE(6,1113)
!
! 
        BETA=BETA*PI/180.0D0 
! 
       IF(Nwave==1)  KP=1
       IF(Nwave .GT. 1)  KP=Nwave/2

! 
       Amp=Ampn(KP)
       W1=FREQ(KP)
       WK=WKN(KP)
       WL=2.0d0*PI/WK !wave number to wave len
       Tper=2.0*PI/FREQ(KP)! wave frewq to wave period
     
       TSTEP=TPER/NTNUM   ! a perid divded to ntnum 
         NTIME=INT(WVSimu*TPER/TSTEP)!how many steps total to simulate
         !wvsimu is how many periods of wave to simulate
      print *,"time info",tstep,ntime 

        WRITE(6,*) 
        WRITE(9,*)
          WRITE(6,*) '                   ================='
        WRITE(9,*) '                   ================='
!
        WRITE(6,1111)  H,Amp,WK,WL,W1,TPER,BETA*180./PI 
        WRITE(9,1111)  H,Amp,WK,WL,W1,TPER,BETA*180./PI 
        1111    format(//,'  water depth=',f9.3,'    wave amplitude=', f6.2,/,&
            &    '  wave number=',f9.5,'  k0=',f9.5,'  wave length=',f9.4,/, &
            &    '  angular frequ.=',f9.5,'   wave period=',f7.3,/,      &
            &    '  wave direction:',f7.3,'  degree',/)
!1111    FORMAT(//,'  WATER DEPTH=',F9.3,'    WAVE AMPLITUDE=', F6.2,/, &
     !&    '  WAVE NUMBER=',F9.5,'  WAVE LENGTH=',F9.4,/, 
     !&    '  ANGULAR FREQU.=',F9.5,'   WAVE PERIOD=',F7.3,/,      
     !&    '  WAVE DIRECTION:',F7.3,'  Degree',/)
! 
1113    FORMAT(/,15x,' FIRST  ORDER PROBLEM')  

        end subroutine

    subroutine waveck(sigma,h,g,wk)
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
