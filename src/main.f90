program hi_project

    use data_all
    use gaussian_info
    use free_term,only:get_free_term,init_ft
    use hi_intg
    use gradient,only:init_gradient
    use wave_funcs_simple,only:poxy,eti
    use io

    implicit  none  
    type(Ostream) :: fstream
    character(len=1024) :: filename
    

    integer :: inode
    real(8) :: xp,yp,zp
    fstream = Ostream("main",[6])
    open(9,  file='output/output1.txt',    status='unknown')
    open(10, file='output/output.txt' ,    status='unknown')
    open(101,file='output/outamt.txt' ,    status='unknown')          
    open(102,file='output/outbmt.txt' ,    status='unknown')

         
    call read_wav_data()
    call fstream%fout('Finished reading wave info')
    call init_time_var
    call output_wav_data()
    xc=0
    yc=0
    zc=0

    !print *,xc,yc,zc
    call read_mesh()

    call fstream%fout('Finished reading mesh info')
    call fstream%fout('Total Surface Node Number:'//fstream%toString(nnf))
    call fstream%fout('Total Normal  Node Number:'//fstream%toString(nnoded))
    call fstream%fout('Total         Node Number:'//fstream%toString(nnode))

    call topology_analysis()
!  --------------------------------------------
    call init_ft(nsys,nnf) 

    call fstream%fout('Initilaized free term vars')
    call init_data(nsys,nnode,nnoded,nnf) 
    call fstream%fout('Initilaized data vars')
    call get_gaussian_data(xc,yc,zc)                  
    call fstream%fout('Generating guassian point info')
    call init_hi_var() 
    call fstream%fout('Initialising hi mod vars')
    !call get_free_term()
    !call tassb0_freq
    call tassb0
    call init_gradient(nnf,nelemf,xyze(1:2,:,1:nelemf),nodele(1:nnf,1),nodelj(1:nnf,1))
    call fstream%fout('Initialising info needed for surface gradient evaluation')
    !call tassbt




    !!! >>===========================================<<

    time=0.0d0
    tstep=0.05
    ntime=100
    
    do itime = 0,ntime
        call fstream%fout(fstream%toString(itime)//'/'//fstream%toString(ntime,'(i5)'))
        time = itime*tstep
        call time_intg_rk4

        !if (mod(itime,10).eq.0) then
        ! elem =491 node =633,1503

        do inode =1,nnf
            xp = xyz(1,inode)
            yp = xyz(2,inode)
            zp = xyz(3,inode)

            !//output all surface nodes for potential and wave elevation
            !write (filename, "(A5,I4)") "phi",itime
            write(7000+itime,1202) inode,bkn(inode,1),poxy(xp,yp,zp)
            !write(filename,1202) inode,bkn(inode,1),poxy(xp,yp,zp)
            !write (filename, "(A5,I4)") "phi",itime
            write(8000+itime,1202) inode,et(inode,1),eti(xp,yp)
            !write(9000+itime,1202) bkn(inode,1),poxy(xp,yp,zp)
        end do

        !//output error on given node
        !inode =288
        !xp = xyz(1,inode)
        !yp = xyz(2,inode)
        !zp = xyz(3,inode)
        !write(6001,'(1f14.6)') bkn(inode,1) - poxy(xp,yp,zp)
        !//------------------------------

        !end do
        !endif


        !/// OUTPUT BODY MESH
        !do inode =nnf+1,nnoded

        !fixme inode -> node
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

    1202 format(i7,5x,2f14.8)
    call fstream%fout("============== main program ends ===============")
end  program      

