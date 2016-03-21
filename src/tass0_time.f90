!C  TASSB0
!C *********************************************************
!C *                                                       *
!C * Calculate the element contribution and assembly the   *
!C * coefficients of the corresponding system of equationn *
!C *                                                       *
!C *********************************************************
    !include './add_on/common_block.f90'

    subroutine tassb0
        use mvar_mod
        use free_term,only:fterm,output_fterms
        use mfunc_mod,only:rludcmp
        use proj_cnst,only:rsn

        implicit   none  
        integer  inode,ielem,j,ip
        integer :: i,is,ii
        real(8)  xp,yp,zp
        real(8)  bmatrix(4,8),amatrix(4,8)

        real(8)  s_angle
        real(8) :: dsign

        call topology_analysis

        amata = 0.0d0
        cmata = 0.0d0

        ! =======================================================================
        !call output_fterms()
        pause
        !print *,"finished fterm output"
        do  500   inode=1,  nnf
                !print *,inode
            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode)

            !fterm_coef = 0
            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
             &                        h,xyz,dxyze,s_angle)    

            angle(inode)=1.0d0- s_angle
            cmata(inode,inode,1:nsys)= -1.0d0+s_angle!angle(inode)
            !amata(inode,inode,1:nsys)= 1.0d0-s_angle!angle(inode)
            !  ---------------------------
            !  Integration on the free surface
                
            do   ielem=1,  nelemf
                call comp_link(ielem,inode,ii)
                if (ii .eq. 0)   then! if src node not on element 
                    call norm_ele0(ielem,xp,yp,zp,amatrix,bmatrix)

                else if (ii .ne. 0)   then 
                    call sing_ele0(inode,ielem,nodqua(inode),xp,yp,zp,&
                        &                   amatrix,bmatrix)
                end if 
                !call testopt(1,1,ielem)
                !stop
                call common_block(0,ielem,inode,amatrix,bmatrix,1)!,fterm_coef)
            end do
            !  Integration on the body surface
            do    ielem=nelemf+1,  nelem

                call comp_link(ielem,inode,ii)
                if (ii .eq. 0)   then 
                    call norm_ele0(ielem,xp,yp,zp,amatrix,bmatrix)
                else if (ii .ne. 0)   then 
                    call sing_ele0(inode,ielem,nodqua(inode),xp,yp,zp,&
                     &                   amatrix,bmatrix)
                end if
                call common_block(1,ielem,inode,amatrix,bmatrix,1)

            end do

            
            do ip = 1,nsys
                    fra3(inode,ip) = fterm(inode,ip,1)!
                    frc31(inode,ip)=fterm(inode,ip,2)-fterm(inode,ip,1)*xp!
                    frc32(inode,ip)=fterm(inode,ip,3)-fterm(inode,ip,1)*yp!
                    frc33(inode,ip)=fterm(inode,ip,4)-fterm(inode,ip,1)*zp!
            write(*,5000) inode,fra3(inode,1),frc31(inode,1),frc32(inode,1),frc33(inode,1)
            end do

            !TODO only works for ip=1
            !write(404,5001) xp,yp,fterm(inode,1,1:4),angle(inode)
            write(405,5000) inode,fra3(inode,1),frc31(inode,1),frc32(inode,1),frc33(inode,1)
            call output_fterms()
            5000 format(I6,4f14.6)
            5001 format(7f14.8)
            


        500     continue
!
        
! =======================================================================
!    Source point is on the body surface
!
        do  1000   inode=nnf+1, nnode   
            !print *,inode

            xp=xyz(1,inode)
            yp=xyz(2,inode)
            zp=xyz(3,inode) 
            
            !fterm_coef = 0
            call solidangle(inode,nnode,nelem,ncn,ncon,nodqua,&
             &                    h,xyz,dxyze,s_angle) 

            angle(inode)=1.0d0-s_angle

            amata(inode,inode,1:nsys)= 1.0d0-s_angle! angle(inode)

            do   ielem=1,  nelemf

                ii=0   
                call norm_ele0(ielem,xp,yp,zp,amatrix,bmatrix)
                call common_block(0,ielem,inode,amatrix,bmatrix,0)

            end do

            do  ielem=1+nelemf, nelem

                call comp_link(ielem,inode,ii)!
                if (ii .eq. 0)   then 
                    call norm_ele0(ielem,xp,yp,zp,amatrix,bmatrix)
                else if (ii .ne. 0)   then 
                    call sing_ele0(inode,ielem,nodqua(inode),xp,yp,zp,&
                     &                   amatrix,bmatrix)
                end if

                !print *,"before common"

                call common_block(1,ielem,inode,amatrix,bmatrix,0)
                !print *,"after common"
            end do
            write(404,5001) xp,yp,angle(inode)
1000     continue
!
! =============================================

!        if( nsys .eq. 2) then
            !do inode=1, nnf
            !if(nodqua(inode) .eq. 2) then
                !amata(inode,inode,2)=1.0e20         
            !endif
            !enddo
        !else if( nsys .eq. 4) then
            !do inode=1, nnf
            !if(nodqua(inode) .eq. 2) then
                !amata(inode,inode,2)=1.0e20
                !amata(inode,inode,4)=1.0e20            
            !else if(nodqua(inode) .eq. 4) then
                !amata(inode,inode,3)=1.0e20
                !amata(inode,inode,2)=1.0e20
            !else if(nodqua(inode) .eq. 5) then
                !amata(inode,inode,2)=1.0e20
                !amata(inode,inode,3)=1.0e20            
                !amata(inode,inode,4)=1.0e20            
            !endif
            !enddo
        !endif

!! =============================================
        print *,"I am here!"
!
        !write(102, *) '  =========== before rludcmp =============='

        do i = 1,nnode
            do j = 1,nnode
                write(400,*) amata(i,j,1:nsys)
        end do;end do

  !      do i = 1,nnode
                !write(401,*) amata(i,i,1:nsys)
        !end do
        write(*,*) "Begin inversing LHS matrix............"
        do ip=1, nsys
            call rludcmp(ip,amata,nnode,nnode,nsys,indx,dsign)  
        enddo
        write(*,*) "Finished inversing LHS matrix............"
        stop

                write(102, *) 
        write(102, *)
        write(102, *) '  =========== after rludcmp =============='

      end subroutine

    subroutine comp_link(ielem,inode,ii) 
        use mvar_mod,only: nodnoe,nodele
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

    subroutine common_block(flag1,ielem,inode,amatrix,bmatrix,is_ft)!,fterm_coef)
        use mvar_mod,only :ncn,ncon,ncond,nsys,nnf,amata,cmata,xyz,dxyz
        use mfunc_mod,only: dinp0
        use free_term,only:fterm
        use proj_cnst,only:rsn,ex,ey

        implicit none
        integer,intent(in) :: ielem,inode,flag1
        real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)
        integer :: is_ft

        integer :: jncon,jnrml,j,ip,is,i
        real(8) :: dpdn,phi,prefix(3),np_j(3),p_j(3)
        real(8) :: dpdx(3)
        
        do     j=1,  ncn(ielem) 
            jncon=ncon(ielem,j)!jth node in ielem
            jnrml=ncond(ielem,j)  !jth nrml in ielem
            do     ip=1, nsys          
                prefix=(/ex(ip),ey(ip),1.0d0/)
                p_j=dot_product(prefix,xyz(:,jncon))
                np_j=dot_product(prefix,dxyz(1:3,jnrml))

                if (flag1.eq.0) then
                    do   is=1, nsys    
                        amata(inode,jncon,ip)=amata(inode,jncon,ip)+rsn(is,ip)*bmatrix(is,j)

                        cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)+rsn(is,ip)*amatrix(is,j)
                    enddo
                else 
                    do  is=1, nsys    
                        !waterline node need special treatment
                        if(jncon .gt. nnf)  then
                            amata(inode,jncon,ip)=amata(inode,jncon,ip)-rsn(is,ip)*amatrix(is,j)
                        else
                            cmata(inode,jncon,ip)=cmata(inode,jncon,ip)+rsn(is,ip)*amatrix(is,j)!*phi2
                        endif
                        cmata(inode,jnrml,ip)=cmata(inode,jnrml,ip)-rsn(is,ip)*bmatrix(is,j)!*dpdn
                    enddo
                end if

                if (is_ft.eq.1) then!
                    do i = 1,4
                        !call dinp0(i,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)       
                        !dpdn=dpox*nx+dpoy*ny+dpoz*nz
                        call dinp0(i,p_j,phi,dpdx)
                        dpdn=dot_product(dpdx,np_j)
                        do    is=1, nsys             
                            !fterm_coef(i,ip)=fterm_coef(i,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                            !fterm_coef(i,ip)=fterm_coef(i,ip)+rsn(is,ip)*amatrix(is,j)*phi
                            fterm(inode,ip,i)=fterm(inode,ip,i)-rsn(is,ip)*bmatrix(is,j)*dpdn        
                            fterm(inode,ip,i)=fterm(inode,ip,i)+rsn(is,ip)*amatrix(is,j)*phi
                        enddo
                    end do
                else 
                    !call dinp0(0,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)      
                    !phi = 1.0d0
                    !do    is=1, nsys             
                    !!fterm_coef(0,ip) = fterm_coef(0,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn
                    !fterm_coef(0,ip) = fterm_coef(0,ip)+rsn(is,ip)*amatrix(is,j)*phi
                    !end do
 
                end if
                end do;end do !j,ip

     end subroutine


     subroutine topology_analysis
         use mvar_mod,only:nnode,ncon,ncn,nodele,nodelj,nodqua,nodnoe,xyz
         use pvar_mod

         integer inode,j,ielem,l

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
