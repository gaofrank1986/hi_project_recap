program hi_project

    use mvar_mod
    use body_property
    use hi_intg

    implicit  none  

    open(6,  file='output/outscreen.txt',    status='unknown')
    open(9,  file='output/output1.txt',    status='unknown')
    open(10, file='output/output.txt' ,    status='unknown')
    open(11, file='output/mesh_track.txt' ,    status='unknown')
    open(12, file='output/output_amatrix.txt' ,    status='unknown')


    open(101,file='output/outamt.txt' ,    status='unknown')          
    open(102,file='output/outbmt.txt' ,    status='unknown')
    open(103,file='output/outcmt.txt' ,    status='unknown')
    open(109,file='output/output9.txt' ,    status='unknown')
    open(108,file='output/output7.txt' ,    status='unknown')
    open(110,file='output/output10.txt' ,    status='unknown')    

         
    call read_wav_data()
    call output_wav_data()
    print *,xc,yc,zc
    xc=0
    yc=0
    zc=0

    call read_mesh()

    write(11,*) ' isys=',isys,' nsys=',nsys
    write(11,*) ' nelemb=',nelemb,' nelemf=',nelemf

    allocate(samb(nelem,16,0:8),sambxy(nelem,16,3))
    allocate(dsamb(nelem,16,6))

    call pre_mesh_2()!
!  --------------------------------------------
    allocate(angle(nnode),fra3(nnode),&
        &          frc31(nnode),frc32(nnode),frc33(nnode))

    allocate(amata(nnode,nnode,nsys),&
        &                 bmata(nnode,nsys), indx(nnode,nsys))
    !       
    allocate(unkn(nnode,nsys),  bkn(nnoded,nsys),&
        &                 unkn_o(nnode,nsys),bkn_o(nnoded,nsys),&
        &             dpdt(nnode,nsys))
    allocate(dh(4,nnf,nsys),dp(4,nnf,nsys),dposi(4,6))

    call bodyfd                  
    call init_hi_var() 
    call tassb0   

    print *,"=================== main program ends ==============="
end  program      

