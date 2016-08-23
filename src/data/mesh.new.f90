!!-------------------------------------------------------------------------------
! DUTWAV
!-------------------------------------------------------------------------------
!  MODULE: mesh
!
!> @brief
!! 
!!
!! @author
!! DUT 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 7.20.2016 : add read damp info  
!
!------------------------------------------------------------------------------- 
module mesh

    use wave 

    implicit none
    ! @var : [xyzb] store body node with dimension(3,node_id) 
    ! @var : [dxyzb] store body nrml with dimension (3,nrml_id)
    real(8),private,allocatable :: xyzb(:,:),dxyzb(:,:)

    integer,private,allocatable :: nconb(:,:),ncondb(:,:)   
    ! nconb => node list body mesh
    ! ncondb => normal list body
    real(8),protected,allocatable :: XYZE(:,:,:),DXYZE(:,:,:),TXYZE(:,:,:)
    ! @var : [xyze]  before combine, new full mesh for node (3,8,elem_id)
    ! @var : [dxyze] before combine, new full mesh for normal (3,8,elem,id)
    real(8),private,allocatable :: XYZTP(:,:),DXYZTP(:,:)
    ! xyztp => convsb, combined node mesh, xyztp(3,node_id)
    ! dxyztp = > combined normal mesh, dxyztp(3,nrml_id)
    integer,protected,allocatable :: NCN(:),NCON(:,:),NCOND(:,:),IETYPE(:)
    ! ncn => elem_type ,namely, the node number in a elem
    ! ncon => combined node list,ncon(8,elem_id)
    ! ncond => combined nrml list,ncond(8,elem_id)
    ! ietype => flag show if a elem is free surface mesh or a body mesh

    integer,private,allocatable :: NNORMN(:),NNORMC(:)
    real(8),protected,allocatable :: xyz(:,:)
    real(8),public,allocatable ::dxyz(:,:)
    integer,protected :: nsys,nelem,nnode,nnoded,isys
    ! nsys => about symmetr
    ! nelem => elem number in combined mesh
    ! nnode => node number in combined mesh
    ! nnoded => nrml nmber in combined mesh
    integer,protected :: nelemb,nnb,nnbd,nelemf,nnf
    ! nelemb => elem num in body
    ! nnb => node num  in body 
    ! nnbd => nrml num in body
    ! nelemf => elem num in fs
    ! nnf => node num in fs

    real(8),protected,allocatable :: dampe(:,:),dampf(:),damptp(:)
    integer,protected,allocatable :: nodele(:,:),nodnoe(:),nodelj(:,:),nodqua(:) 
contains
    subroutine read_mesh()
        implicit none

        INTEGER IFWKO,IPOL

    
        open(2, FILE='INPUT/DATBDMS.txt',    status='old') 
        open(3, FILE='INPUT/DATWPMS.txt',    status='old')

        !==================body mesh================================
        READ(2,*)   ISYS 
        READ(2,*)   NELEMB, NNB, NNBD, IPOL

        IF(ISYS.EQ.0) NSYS=1
        IF(ISYS.EQ.1) NSYS=2
        IF(ISYS.EQ.2) NSYS=4

        READ(3,*)   NELEMF

        NELEM=NELEMB+NELEMF

        allocate (nconb(nelemb,8),ncondb(nelemb,8))!body node/normal list
        allocate (xyzb(3,nnb),dxyzb(3,nnbd))    !body node/normal
        allocate(ncn(nelem),ncon(nelem,8),ncond(nelem,8),ietype(nelem))

        allocate(nnormn(8*nelem) )

        allocate( xyze(3,8,nelem),dxyze(3,8,nelem))
        allocate( txyze(3,8,nelem))
        allocate(xyztp(3,8*nelem),dxyztp(3,8*nelem))
        allocate(dampe(8,nelem),damptp(8*nelem))

        call read_surface(1) !read in data on free surface mesh
        call read_body(ipol) ! read in data on body mesh

        close(2)
        close(1)
        close(3)

        call merge_mesh()
        call assign_info()
    end subroutine


    subroutine assign_info()
        implicit none
        integer :: ndmax,ind
        ndmax=max(nnode,nnoded)
        allocate(nodele(nnode,64),nodnoe(nnode),nodelj(nnode,64),&
            &         nodqua(nnode),nnormc(nnoded),&
            &           xyz(3,ndmax),dxyz(6,ndmax),dampf(nnf))
        !

        do ind=1, nnf
            dampf(ind)=dampe(ind)
            xyz(:,ind)=xyztp(:,ind)
        end do

        do ind=nnf+1, nnode
            xyz(:,ind)=xyztp(:,ind)
        end do

        do ind=1, nnoded
            dxyz(:,ind)=dxyztp(:,ind)
            nnormc(ind)=nnormn(ind)
        end do

        deallocate(xyztp,dxyztp,nnormn)

        dampf(:)=w1*dampf(:)
    end subroutine

    ! @func : read in data file of free surface mesh
    ! @param: [flag] 0 if no damping info, 1 if damp info presents

    subroutine read_surface(flag)
        implicit none
        integer :: ie,j,m,flag

        xyze=0.0d0
        do  ie=1, nelemf
            ietype(ie)=2
            read(3, *)    m, ncn(ie)
            read(3, *) (xyze(1,j,ie), j=1, ncn(ie))
            read(3, *) (xyze(2,j,ie), j=1, ncn(ie))
            if (flag.eq.0) then
                dampe(:,:)=0.0d0
            else if(flag.eq.1) then
                read(3, *) (dampe(j,ie), j=1, ncn(ie))
            end if

        end do
        ! in dat file ,positive nrml is set as pointing into fulid field
        dxyze(1, 1:8, 1:nelemf)= 0.0d0
        dxyze(2, 1:8, 1:nelemf)= 0.0d0
        dxyze(3, 1:8, 1:nelemf)= 1.0d0 
        !     DXYZE(3, 1:8, 1:NELEMF)=-1.0d0 
    end subroutine

    !C ***************************************************************
    !C *                                                             *
    !C *  Generate nodal and element data on the body surface of     *
    !C *  an arbitrary body                                          *
    !C *                                                             *
    !C ***************************************************************
    !C 
    subroutine read_body(ncor)
        use kinds
        implicit none

        integer ncor,ipl,ipolar(50)
        integer i,ie,m,inode,ncnn,k,kk
        real*8  a1,a2,r1,r2,z1
        real*8  xofset(50),yofset(50),zofset(50)
        real(rk) :: pi=4*atan(1.0d0)

        do i=1, ncor
            read(2,*)  m, ipolar(i), xofset(i),yofset(i),zofset(i)
        end do

        do  inode=1, nnb
            read(2,*) m, ipl, r1, a1, z1
            if(ipolar(ipl).eq.0) then
                xyzb(1,inode)=r1+xofset(ipl)
                xyzb(2,inode)=a1+yofset(ipl)
            else if(ipolar(ipl).eq.1) then
                xyzb(1,inode)=r1*dcos(a1*pi/180.0d0)+xofset(ipl)
                xyzb(2,inode)=r1*dsin(a1*pi/180.0d0)+yofset(ipl)
            endif
            xyzb(3,inode)=z1+zofset(ipl)
        end do

        do  inode=1, nnbd
            read(2,*) m,ipl,r2,a2,dxyzb(3,inode)

            if(ipolar(ipl).eq.0) then
                dxyzb(1,inode)= r2
                dxyzb(2,inode)= a2
            else if(ipolar(ipl).ne.0) then
                dxyzb(1,inode)=r2*dcos( a2*pi/180.0d0 )
                dxyzb(2,inode)=r2*dsin( a2*pi/180.0d0 )
            endif
        end do

        do  i=1, nelemb
            ie=i+nelemf
            ietype(ie)=1
            read(2,*) m, ncn(ie)    
            read(2,*) (nconb(i,k), k=1, ncn(ie))
        end do

        do  i=1,  nelemb
            ie=i+nelemf
            read(2,*)  m, ncnn
            read(2,*) (ncondb(i,k), k=1, ncn(ie))
        end do!
        ! ===============================================================
        do i=1, nelemb
            ie=i+nelemf
            do  k=1, ncn(ie)
                xyze(:,k,ie)=xyzb(:, nconb(i,k))
                dxyze(:,k,ie)=dxyzb(:,ncondb(i,k))
            end do!
        end do!
    end subroutine
 
    ! @var nmat : node matrix shape (3,node num)
    ! @var [n] : current nodes stores in nmat
    function get_node_id(p,nmat,n) result(ans)
        use kinds
        implicit none
        real(rk) :: p(3),nmat(3,:),tol
        integer :: n,i,ans
        tol=1e-6 
        ans=-1

        do i=1,n
            if (norm2(p-nmat(:,i))<tol) then
                ans=i
                exit!stop closest loop
            end if
        end do
    end function

     function get_nrml_id(base,p,nmat,mn,n) result(ans)
        use kinds
        implicit none
        real(rk) :: p(3),nmat(3,:),tol
        integer :: n,i,ans,mn(:),base
        tol=1e-6 
        ans=-1

        do i=1,n
            if ((norm2(p-nmat(:,i))<tol).and.(base.eq.mn(i))) then
                ans=i
                exit!stop closest loop
            end if
        end do
    end function           
        



    subroutine merge_mesh(tol)

        use kinds
        implicit none

        integer ieb,iel,ind,i,l,m,nnd0,kk,n
        real(8) x,y,z,dx,dy,dz,dr2,tmp1,tmp2,tmp3,tol
        real(rk) :: p(3),dp(3)
        real(rk),optional :: tol
        integer :: id,cur_avail_nid

        if(present(tol)) then
        else
            tol=1.0e-6
        end if

        next_node_id=1
        next_nrml_id=1

        if (nelemf .eq. 0) then
            nnf=0
        else
            do iel=1,  nelemf
                do i=1, ncn(iel)
                    p=xyze(1:3,i,iel) ! eid element/ith node/xyz
                    !========check coincident node
                    id = get_node_id(p,xyztp,next_node_id-1)
                    ! cannot find existing node, create new node
                    if(id.ne.-1) then
                        ncon(iel, i)=next_node_id
                        xyztp(1:3,next_node_id) = p
                        if(iel.le.nelemf) then
                            damptp(next_node_id) = dampe(i,iel)
                        end if
                        next_node_id = next_node_id + 1
                    else
                        ncon(iel,i) = id
                        !if(dabs(dampe(i,iel)-damptp(id))>tol) then
                            !print *,"Error!damp info mismatch!"
                        !end if
                    end if

                    dp=dxyze(:,i,iel)
                    id2 = get_nrml_id(ncon(iel,i),dp,dxyztp,nnormn,next_nrml_id-1)

                    if(id2.ne.-1) then
                        ncond(iel, i)=next_nrml_id
                        dxyztp(1:3,next_nrml_id) = dp
                        nnormn(next_nrml_id)=ncon(iel,i)
                        !node nrml relationship nnormn(j) = i  nrml id j has base node i 
                        next_nrml_id = next_nrml_id + 1
                    else
                        ncond(iel,i) = id2
                        nnormn(next_nrml_id)=con(iel,i)
                        !node nrml relationship nnormn(j) = i  nrml id j has base node i 
                    end if

                end do
                if(iel.eq.nelemf) nnf = next_node_id - 1
            end do
        end if    
        !
        ! -----------------------------------------
        !
        do ieb=1,  nelemb
            iel=ieb+nelemf
            do i=1, ncn(iel)

                p=xyze(:,i,iel)
                dp=dxyze(:,i,iel)
                id = get_nid(p,xyztp,next_node_id-1)
                id2 = get_nid(dp,dxyztp,next_nrml_id-1)
                !C
                !C ** if the point (XPM, YPM) does not concide with any other
                !C    points had been checked,  give a new code

                if(id.ne.-1) then
                    ncon(iel, i)=next_node_id
                    xyztp(1:3,next_node_id) = p
                    next_node_id = next_node_id + 1
                else
                    ncon(iel,i) = id
                end if

                if(id2.ne.-1) then
                    ncond(iel, i)=next_nrml_id
                    dxyztp(1:3,next_nrml_id) = dp
                    nnormn(next_nrml_id)=ncon(iel,i)
                    !node nrml relationship nnormn(j) = i  nrml id j has base node i 
                    next_nrml_id = next_nrml_id + 1
                else
                    ncond(iel,i) = id2
                    nnormn(next_nrml_id)=con(iel,i)
                    !node nrml relationship nnormn(j) = i  nrml id j has base node i 
                end if
            end do
        end do

        nnode = next_node_id - 1
        nnoded = next_nrml_id - 1
        nnb=nnode-nnf

    end subroutine

    subroutine output_mesh(filename)
        implicit none

        integer :: n,iel,i
        character(len=*) :: filename
        open(102,file=filename)


        write(102,'(i6)') isys 
        write(102,'(4i6)')  nelem, nnode,nnoded,1 
        write(102,'(a)') ' 1  0  0.0  0.0  0.0' 

        do n=1, nnode 
            write(102,1010) n, 1, xyz(1,n),xyz(2,n),xyz(3,n) 
        enddo 

        do n=1, nnoded 
            write(102,1010) n, 1, dxyz(1,n),dxyz(2,n),dxyz(3,n) 
        enddo 

        do iel=1, nelem 
            write(102,1001) iel, ncn(iel) 
            write(102,1005) (ncon(iel, i), i=1, ncn(iel)) 
        enddo 

        do iel=1, nelem 
            write(102,1001) iel, ncn(iel) 
            write(102,1005) (ncond(iel, i), i=1, ncn(iel)) 
        enddo 
        close(102)
        1001   format(1x,2i6,3f14.6) 
        1010   format(1x,i6,i4,3f14.6) 
        1005   format(8(1x,i6)) 
        end subroutine

    subroutine output_surface(filename)
        implicit none

        integer :: n,iel,i
        character(len=*) :: filename
        open(102,file=filename)


        write(102,'(i6)') isys 
        write(102,'(4i6)')  nelemf, nnf,nnf,1 
        write(102,'(a)') ' 1  0  0.0  0.0  0.0' 

        do n=1, nnf 
            write(102,1010) n, 1, xyz(1,n),xyz(2,n),xyz(3,n) 
        enddo 

        do n=1, nnf 
            write(102,1010) n, 1, dxyz(1,n),dxyz(2,n),dxyz(3,n) 
        enddo 

        do iel=1, nelemf
            write(102,1001) iel, ncn(iel) 
            write(102,1005) (ncon(iel, i), i=1, ncn(iel)) 
        enddo 

        do iel=1, nelemf
            write(102,1001) iel, ncn(iel) 
            write(102,1005) (ncond(iel, i), i=1, ncn(iel)) 
        enddo 
        close(102)
        1001   format(1x,2i6,3f14.6) 
        1010   format(1x,i6,i4,3f14.6) 
        1005   format(8(1x,i6)) 
        end subroutine

    subroutine output_damp_info(filename)
        implicit none

        integer :: n,iel,i
        character(len=*) :: filename
        open(102,file=filename)

        do n=1, nnf 
            write(102,'(i6,6x,f20.15)') n, dampf(n)
        enddo 
        close(102)

        end subroutine

    subroutine comp_link(ielem,inode,ii) 
        !use mvar_mod
        ! return how many nodes in this element 
        implicit none
        integer,intent(in) :: inode,ielem
        integer,intent(out) :: ii

        integer :: i

        ii = 0 
        do i=1, nodnoe(inode)!!list of linked element by inode (how manys times inode appear in element)
        if(ielem .eq. nodele(inode,i)) then!
          ii=ii+1
        endif
        enddo
    end subroutine 

    logical function is_connected(ielem,inode) 
        !use mvar_mod
        ! return how many nodes in this element 
        implicit none
        integer,intent(in) :: inode,ielem
        !logical :: is_connected

        integer :: i,ii

        ii = 0 
        do i=1, nodnoe(inode)!!list of linked element by inode (how manys times inode appear in element)
        if(ielem .eq. nodele(inode,i)) then!
          ii=ii+1
        endif
        enddo
        if (ii==0) then
            is_connected=.false.
        else
            is_connected=.true.
        endif
    end function 

    subroutine topology_analysis()
        !use mvar_mod
        implicit none
        integer :: inode,ielem,j,l
        do 50 inode=1, nnode 
        l=0
        do 40 ielem=1,  nelem
        do 30 j=1,      ncn(ielem)
        if(inode.eq.ncon(ielem,j)) then
        l=l+1
        nodele(inode,l)=ielem! elem num linked to inode
        nodelj(inode,l)=j!in node-linked-element, inode appear as j-th node 
        endif
30      continue
40      continue
              nodnoe(inode)=l !total number of links
!        below related to symmetry
        nodqua(inode)=0
        if( nsys .ge. 2) then
          if( abs(xyz(2,inode)).lt.1.0e-06 ) then
          nodqua(inode)=2
          end if
        end if
!
        if( nsys .eq. 4) then
          if( abs(xyz(1,inode)).lt.1.0e-06.and.&
     &        abs(xyz(2,inode)).lt.1.0e-06) then
           nodqua(inode)=5
          else if( abs(xyz(1,inode)).lt.1.0e-06 ) then
           nodqua(inode)=4
          endif
        end if
!
50      continue
        PRINT *,"topology analysis finished"
    end subroutine
   end module
