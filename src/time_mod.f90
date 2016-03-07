module time_mod
    implicit none
    real(8) :: ampn(1),phi_w(1),wkn(1),freq(1)
    real(8) :: time,tstep
    integer :: nfreq,nwave,iorder
contains
    subroutine init_time_var()
        use wave,only:wk,w1,amp
        implicit none
        ampn(1) = amp
        phi_w(1) = 0.0d0
        wkn(1) = wk
        freq(1) = w1
        nwave = 1
        iorder = 1
        nfreq = 1

    end subroutine
end module
