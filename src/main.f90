program hi_project

    use mvar_mod
    use pvar_mod
    use free_term,only:get_free_term,init_ft
    !use body_property
    use hi_intg
    use gradient,only:init_gradient

    implicit  none  
    integer inode

    open(9,  file='output/output1.txt',    status='unknown')
    open(10, file='output/output.txt' ,    status='unknown')
    open(101,file='output/outamt.txt' ,    status='unknown')          
    open(102,file='output/outbmt.txt' ,    status='unknown')
    open(103,file='INPUT/init_cond.txt' ,    status='old')

         
    call read_wav_data()
    !call output_wav_data()
    print *,xc,yc,zc
    xc=0
    yc=0
    zc=0

    call read_mesh()
!  --------------------------------------------
    call init_ft(nsys,nnf) 
 
    call bodmass()

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

    call get_gaussian_data(xc,yc,zc)                  
    call init_hi_var() 
    call get_free_term()
    !call tassb0_freq
    call tassb0
    call init_gradient(nnf,nelemf,xyze(1:2,:,1:nelemf),nodele(1:nnf,1),nodelj(1:nnf,1))
    !call tassbt
    print *,"tassb0 ended"
    itime=0
    time=0.0d0

    et(:,:)  =0.0
    bkn(:,:) =0.0
    unkn(:,:)=0.0
    unkn_o(:,:)=0.0
    force(:) =0.0
    disp(:)  =0.0
    dsdt(:)  =0.0

    do inode=1,nnf
        read(103,*) bkn(inode,1),et(inode,1)
        bkn(inode,1)=bkn(inode,1)*(1-dampf(inode)/w1)
        et(inode,1)=et(inode,1)*(1-dampf(inode)/w1)
    end do
    tstep=0.05 
    print *,"before time step"
    do 500 itime=0,0! endtime!0!ntime
    print *,"Calculating",itime,"/",endtime
    write(10,*)
    write(10,*) '  itime=',itime
    write(11,*)
    write(11,*) '  itime=',itime

    time=itime*tstep
    !if(mod(itime,10) .eq. 0)  then
            !write(6,1115) itime,time     
    !endif

    bkn_o(:,:) = bkn(:,:)
    et_o(:,:)  = et(:,:)  

    disp_o(:) =disp(:)
    dsdt_o(:)= dsdt(:)

    !c                  force_o=force
    call time_intg_rk4
    500      continue      
    print *,"=================== main program ends ==============="
end  program      

