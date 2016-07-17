program hi_project

    use data_all
    use gaussian_info
    use free_term,only:get_free_term,init_ft
    use hi_intg
    use gradient,only:init_gradient
    use wave_funcs_simple,only:poxy,eti

    implicit  none  

    integer :: inode
    real(8) :: xp,yp,zp

    open(9,  file='output/output1.txt',    status='unknown')
    open(10, file='output/output.txt' ,    status='unknown')
    open(101,file='output/outamt.txt' ,    status='unknown')          
    open(102,file='output/outbmt.txt' ,    status='unknown')

         
    call read_wav_data()
    call init_time_var
    call output_wav_data()
    xc=0
    yc=0
    zc=0

    print *,xc,yc,zc
    call read_mesh()
!  --------------------------------------------
    call init_ft(nsys,nnf) 
    call init_data(nsys,nnode,nnoded,nnf) 

    call get_gaussian_data(xc,yc,zc)                  
    call init_hi_var() 
    !call get_free_term()
    !call tassb0_freq
    call tassb0
    call init_gradient(nnf,nelemf,xyze(1:2,:,1:nelemf),nodele(1:nnf,1),nodelj(1:nnf,1))
    !call tassbt




    !!! >>===========================================<<

    time=0.0d0
    tstep=0.05
    
    do itime = 0,20
        print *,itime,'/200'
        time = itime*tstep
    call time_intg_rk4
  !  do inode =1,nnf
        !xp = xyz(1,inode)
        !yp = xyz(2,inode)
        !zp = xyz(3,inode)
     !write(4003,1202) bkn(inode,1),poxy(xp,yp,zp)
     !write(4004,1202) et(inode,1),eti(xp,yp)
     !end do
     !if (mod(itime,10).eq.0) then
        do inode =1,nnf
            xp = xyz(1,inode)
            yp = xyz(2,inode)
            zp = xyz(3,inode)
            
            write(7000+itime,1202) inode,bkn(inode,1),poxy(xp,yp,zp)
            
            write(8000+itime,1202) inode,et(inode,1),eti(xp,yp)
        end do
      !endif
        !do inode =nnf+1,nnoded
            !xp = xyz(1,inode)
            !yp = xyz(2,inode)
            !zp = xyz(3,inode)
            !write(7000+itime,1202) bkn(inode,1),poxy(xp,yp,zp)
        !end do
        bkn_o=bkn
        et_o=et
        !do inode =1,nnf
            !xp = xyz(1,inode)
            !yp = xyz(2,inode)
            !zp = xyz(3,inode)
            !write(7000+itime,1202) bkn(inode,1),poxy(xp,yp,zp)
            
            !write(8000+itime,1202) et(inode,1),eti(xp,yp)
        !end do

        !do inode =nnf+1,nnode
            !xp = xyz(1,inode)
            !yp = xyz(2,inode)
            !zp = xyz(3,inode)
            !write(8500+itime,1202) unkn(inode,1),poxy(xp,yp,zp)
        !end do
    end do

     1202 format(i6,5x,2f14.8)
    print *,"============== main program ends ==============="
end  program      

