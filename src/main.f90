program hi_project

    use mvar_mod
    use free_term,only:get_free_term,init_ft
    use body_property
    use hi_intg
    use gradient,only:init_gradient
    use time_mod
    use mfunc_mod,only:poxy,eti

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
    print *,xc,yc,zc
    xc=0
    yc=0
    zc=0

    call read_mesh()
!  --------------------------------------------
    call init_ft(nsys,nnf) 
 

    allocate(angle(nnode))
    allocate(fra3(nnode,nsys),&
        &          frc31(nnode,nsys),frc32(nnode,nsys),frc33(nnode,nsys))

    allocate(amata(nnode,nnode,nsys),&
        &                 bmata(nnode,nsys), indx(nnode,nsys))
    allocate(cmata(nnode,nnoded,nsys))
    allocate(nrml_2_node(nnoded))
    amata=0
    bmata=0
    cmata = 0
    allocate(unkn(nnode,nsys))


    allocate(bkn(nnoded,nsys),&
        &                 unkn_o(nnode,nsys),bkn_o(nnoded,nsys),&
        &            et(nnf,nsys),et_o(nnf,nsys), dpdt(nnode,nsys))
    allocate(dh(4,nnf,nsys),dp(4,nnf,nsys),dposi(4,6))
    unkn=0.
    unkn_o=0.
    bkn_o=0.
    et=0.
    et_o=0.
    dh=0.
    dp=0.
    dpdt=0.
    

    call get_gaussian_data(xc,yc,zc)                  
    call init_hi_var() 
    call get_free_term()
    !call tassb0_freq
    call tassb0
    call init_gradient(nnf,nelemf,xyze(1:2,:,1:nelemf),nodele(1:nnf,1),nodelj(1:nnf,1))
    !call tassbt
    time=0.0d0
    tstep=0.05
    call time_intg_rk4
  !  do inode =1,nnf
        !xp = xyz(1,inode)
        !yp = xyz(2,inode)
        !zp = xyz(3,inode)
     !write(4003,1202) bkn(inode,1),poxy(xp,yp,zp)
     !write(4004,1202) et(inode,1),eti(xp,yp)
     !end do
        do inode =1,nnf
            xp = xyz(1,inode)
            yp = xyz(2,inode)
            zp = xyz(3,inode)
            write(7000+itime,1202) bkn(inode,1),poxy(xp,yp,zp)
            write(8000+itime,1202) et(inode,1),eti(xp,yp)
        end do

     1202 format(2f14.8)
    print *,"=================== main program ends ==============="
end  program      

