module wave
    implicit none

    ! @var [wk] : wave number
    ! @var [h]  : water depth
    ! @var [beta] : incident wave angle
    ! @var [amp] : wave amplitude
    ! @var [w1]  : wave frequency
    ! @var [tper] : wave period
    ! @var [v] : wave number deep water ??
    ! @var [wl] : wave length ?

    real(8),protected :: wk,h,beta,amp
    real(8),protected :: w1,tper,v,wl 
    real(8),parameter :: g = 9.807
    real(8),parameter :: rho = 1.023e3
    real(8),parameter :: pi = 3.14159265359 
    real(8) :: timerk,rampf
    !v => wav num deep water
    !type :: single_wav
        !real(rk) :: amp,beta,freq,per,wlen

contains

    subroutine read_wav_data()

        implicit none

        integer :: ifwko

        open(1, file='INPUT/DATIN.txt', status='old') 

        read(1,*)      ifwko
        read(1,*)      h, amp, wk, beta

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

    subroutine output_wav_data()

        implicit none

        open(1, file='output/wave_info.txt',      status='unknown') 
        write(1,*) 
        write(1,*) '                   ================='
        write(1,*) " gravity =",g
        write(1,*) " rho     =",rho
        write(1,*) " pi      =",pi  
        write(1,1111)  h,amp,wk,v,wl,w1,tper,beta*180./pi 

        1111    format(//,'  water depth=',f9.3,'    wave amplitude=', f6.2,/,&
            &    '  wave number=',f9.5,'  wave length=',f9.4,/, &
            &    '  angular frequ.=',f9.5,'   wave period=',f7.3,/,      &
            &    '  wave direction:',f7.3,'  degree',/)
        close(1)
    end subroutine 

    ! @func : compute wk given (h,sigma)
    ! @param : [h] water depth, for infinite depth use negative value
    ! @param :[sigma] wave frequency
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
!C  ******************************************************* 
!C  *                                                     * 
!C  *    The subroutine computes the real and imaginary   * 
!C  *  roots of dispersion equation by wave frequency     * 
!C  *  W and water depth H.                               * 
!C  *                                                     * 
!C  ******************************************************* 

      SUBROUTINE DISPERS(WAVENO,MNO,W,H) 
          !C 
          !CCC   EVALUATION OF THE ROOTS OF THE FOLLOWING EQUATIONS 
          !CCC   BY NEWTON-RAPHSON METHOD,RESULTS ARE GIVEN IN ARRAY WAVENO 
          !CCC   MNO ARE THE NUMBER OF ROOTS REQUIRED, FIRST ROOT IS FROM EQN. (I) 
          !CCC   THE REST ARE THE FIRST +VE (MNO-1) ROOTS OF (II) 
          !CCC   I) W*W/G = K TANH( KH ) 
          !CCC   II) -W*W/G = M TAN( MH ) 
          !C 
 
          IMPLICIT NONE 
          INTEGER,INTENT(IN)::MNO 
          REAL*8, INTENT(IN)::W,H 
          REAL*8, INTENT(OUT)::WAVENO(800) 

          INTEGER I,M,MM,IFLAG 
          REAL*8 UPLIM,LOLIM,G,PI 
          REAL*8 WWH,FUN,DFUN,TRIAL,EXX,EXXOR,CHECK 


          DATA G,PI/9.807d0,3.141592653589793d0/ 

          IF(MNO .GT. 799) THEN 
              Print *,' MNO=',MNO,' To enlarge WAVENO(*)' 
              STOP 
          END IF 


          DO 1 I=1,9 
              1 WAVENO(I)=0.D0 
              WWH=W*W*H 

              !CCC   CALCULATION OF WAVE NUMBER (ROOT OF EQN. (I)) 

      TRIAL=WWH/G 
      M=0 
      IF (TRIAL.GT.1.D1) GO TO 20 
      EXX=0.D0 
      IFLAG=0 
   10 FUN=G*TRIAL - WWH/DTANH(TRIAL) 
      DFUN=G + WWH/DSINH(TRIAL)/DSINH(TRIAL) 
      TRIAL=TRIAL - FUN/DFUN 
      EXXOR=DABS(TRIAL - EXX) 
      IF (EXXOR.LE.1.0D-10) GO TO 20 
      EXX=TRIAL 
      GO TO 10 
   20 MM=M + 1 
      WAVENO(MM)=TRIAL/H 
      CHECK=DABS(W*W/G - WAVENO(MM)*DTANH(TRIAL)) 
      IF (CHECK.GT.1.0D-5) GO TO 999 
      IF (MNO.LE.1) RETURN 

!CCC   CALCULATION OF FIRST +VE (MNO-1) ROOTS OF EQN. (II) 

      M=1 
      IFLAG=1 
      EXX=0.D0 
      IF (WWH.LE.2.25D2) GO TO 120 
      GO TO 110 
  100 M=MM 
      EXX=0.D0 
      IF (MM.EQ.MNO) GO TO 9999 
      IF (IFLAG.EQ.2) GO TO 120 
  110 TRIAL=(DBLE(FLOAT(M)) - 1.D-1)*PI 
      GO TO 140 
  120 IFLAG=2 
      TRIAL=(DBLE(FLOAT(M)) - 5.D-1)*PI + 1.D-1 
  140 IF(IFLAG .EQ. 1)GO TO 160 
      IF(IFLAG .EQ. 2)GO TO 170 
  150 TRIAL=TRIAL - FUN/DFUN 
      EXXOR=DABS(TRIAL - EXX) 
      IF (EXXOR.LE.1.D-10) GO TO 180 
      EXX=TRIAL 
      IF(IFLAG .EQ. 2)GO TO 170 
  160 FUN=G*TRIAL + WWH/DTAN(TRIAL) 
      DFUN=G - WWH/DSIN(TRIAL)/DSIN(TRIAL) 
      GO TO 150 
  170 FUN=WWH/(G*TRIAL) + DTAN(TRIAL) 
      DFUN=-WWH/(G*TRIAL*TRIAL) + 1.D0/DCOS(TRIAL)/DCOS(TRIAL) 
      GO TO 150 
  180 UPLIM=(DBLE(FLOAT(M)) + 5.D-1)*PI 
      LOLIM=UPLIM - PI 
      IF ((TRIAL.GT.UPLIM).OR.(TRIAL.LT.LOLIM)) GO TO 190 
      MM=M + 1 
      WAVENO(MM)=TRIAL/H 
      CHECK=DABS(W*W/G + WAVENO(MM)*DTAN(TRIAL)) 
      IF (CHECK.GT.1.0D-5) GO TO 999 
      GO TO 100 
  190 IF (IFLAG.EQ.1) GO TO 120 
      WRITE(6,200)M 
  200 FORMAT('  **ERROR**  OCCURS AT M',I3) 
      STOP 
  999 WRITE(6,1000)CHECK 
      GO TO 190 
 1000 FORMAT('  CHECK=',D11.4) 
 9999 CONTINUE 
 
      RETURN 
      END 

end module 
