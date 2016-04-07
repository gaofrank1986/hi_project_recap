module motion
      implicit none
      real(8) :: xc,yc,zc

          INTEGER  NNT
          PARAMETER (NNT=2000)    

        REAL*8 XTC,YTC,ZTC
        REAL*8 ARE,XF,YF,XK2,YK2,XCF
        REAL*8 VOLM,XB,YB,ZB

        REAL*8 AMAS(6,6),BDMP(6,6)  
        REAL*8 RMAS(6,6),CDMP(6,6),CRS(6,6),STKM(6,6),XIA(3,3)

        REAL*8 FORCEW(6),FORCE0(6),FORSCD(6),AMPJ(6)
          REAL*8 FORCE(6),FORCE_O(6)
          
        REAL*8 TRMAS(6,6),VISC(6,6)

        REAL*8  DISP(6),DSDT(6),DISP_O(6)
        real(8) ::  DSDT_O(6),DSDDTL(6),TRXYZ(3,3)   

        REAL*8  RESPR(6),RESPI(6)  
!      xc = 1
!      yc = 2
!      zc = 3
end module
