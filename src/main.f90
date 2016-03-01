program hi_project

    use mvar_mod
    !use motion
    use free_term,only:fterm,get_free_term,init_ft
    use pvar_mod
    use hi_intg
    use gradient,only:init_gradient

    implicit  none  

    real(8) Stime,Etime
    integer :: i,j,k,inode

    open(9,  file='output/output1.txt',    status='unknown')
    open(10, file='output/output.txt' ,    status='unknown')
    open(101,file='output/outamt.txt' ,    status='unknown')          
    open(102,file='output/outbmt.txt' ,    status='unknown')

         
    call read_wav_data()
    !call output_wav_data()
    !print *,xc,yc,zc
    !xc=0
    !yc=0
    !zc=0

    call read_mesh()
!  --------------------------------------------

    call init_ft(nsys,nnf) 

    allocate(angle(nnode))
    !time domian: var def changed
    !allocate(fra3(nnode),&
        !&          frc31(nnode),frc32(nnode),frc33(nnode))

    allocate(fra3(nnode,nsys),&
        &          frc31(nnode,nsys),frc32(nnode,nsys),frc33(nnode,nsys))

    allocate(amata(nnode,nnode,nsys),&
        &                 bmata(nnode,nsys), indx(nnode,nsys))
    allocate(cmata(nnode,nnoded,nsys))
    allocate(unkn(nnode,nsys))


    allocate(bkn(nnoded,nsys),&
        &                 unkn_o(nnode,nsys),bkn_o(nnoded,nsys),&
        &            et(nnf,nsys),et_o(nnf,nsys), dpdt(nnode,nsys))

    allocate(dh(4,nnf,nsys),dp(4,nnf,nsys),dposi(4,6))

    call bodmass       
    call get_gaussian_data(xc,yc,zc)                  
    call init_hi_var() 
    call get_free_term()
    call tassb0   
    !do i=1,nnf 
        !print *,i
        !print *,nodelj(:,inode)
    call init_gradient(nnf,nelemf,xyze(1:2,:,1:nelemf),nodele(1:nnf,1),nodelj(1:nnf,1))
    !do inode =1,nnf
!  print *,"inode =",inode
        !print *,nodele(inode,1:5)
        !print *,nodelj(inode,1:5)
        !pause
        !end do

    print *,"tassb0 ended"
    !do ip = 1,nsys
            !fra3(inode,ip)=fterm_coef(0,ip,inode)
            !frc31(inode,ip)=fterm_coef(1,ip,inode)-fterm_coef(0,ip)*xp!
            !frc32(inode,ip)=fterm_coef(2,ip)-fterm_coef(0,ip)*yp!
            !frc33(inode,ip)=fterm_coef(3,ip)-fterm_coef(0,ip)*zp!
    !end do
    itime=0
    time=0.0d0

    et(:,:)  =0.0
    bkn(:,:) =0.0
    unkn(:,:)=0.0
    unkn_o(:,:)=0.0
    force(:) =0.0
    disp(:)  =0.0
    dsdt(:)  =0.0
    !
    !         CALL PLOTOUT8
    print *,"before time step"
    do 500 itime=0, 0!ntime
    write(10,*)
    write(10,*) '  itime=',itime
    write(11,*)
    write(11,*) '  itime=',itime

    time=itime*tstep
    if(mod(itime,10) .eq. 0)  then
            write(6,1115) itime,time     
    endif

    bkn_o(:,:) = bkn(:,:)
    et_o(:,:)  = et(:,:)  

    disp_o(:) =disp(:)
    dsdt_o(:)= dsdt(:)

    !c                  force_o=force
    call time_intg_rk4
    500      continue      
    deallocate(amata,cmata,bmata,indx)


    call cpu_time(etime)

    write(9,*)
    write(9,*) '    Ending Time=',etime
    write(9,*) '    Time used=',etime-stime


1115      FORMAT(' I_time=',I5,'    at the time:',F10.5,' s')
    print *,"=================== main program ends ==============="
end  program      

