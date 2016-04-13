!  TASSBT+ TINBODT+ NSBDT
!
!C *********************************************************
!C *                                                       *
!C * Calculate the element contribution and assembly the   *
!C * coefficients of the corresponding system of equationn *
!C *                                                       *
!C *********************************************************
!C
    subroutine tassbt
  
      use data_all
      use motion
      use linalg,only:rlubksb
      use wave_funcs_simple,only:dinp,poxy
      use gradient,only:eval_gradient
      use free_term,only:fterm

      implicit none  

      integer  inode,ind,is,ip,j,i,n1,n2
      real(8) xp,yp,zp,dpox,dpoy,dpoz,dpdn,nx,ny,nz
      real(8) rsn(4,4),ex(4),ey(4)
      real(8) bmat(nsys),cmat(nnoded,nsys)
      real(8) :: dpoxyz_save(2,4,nnoded),tmp(8),tmp2(2,1)

      data rsn /1.,  1.,  1.,  1.,& 
     &            1., -1.,  1., -1.,&
     &            1.,  1., -1., -1.,& 
     &            1., -1., -1.,  1./ 
      data ex / 1.,  1., -1., -1./                                     
      data ey / 1., -1., -1.,  1./

       print *,"Entering tassbt",timerk


      ! establish nrml to node relationship
      do i=1,nelem
          do j=1,ncn(i)
              n1 = ncon(i,j)
              n2 = ncond(i,j)
              nrml_2_node(n2) = n1
      enddo;enddo


      cmat(:,:)=0.0
      dpoxyz_save = 0.0

      do inode =1,nnf
          xp=ex(ip)*xyz(1,inode)
          yp=ey(ip)*xyz(2,inode)
          zp=       xyz(3,inode)

          !bkn(inode,1) =poxy(xp,yp,zp) 
          call dinp(xp,yp,zp,dpox,dpoy,dpoz)
          dpoxyz_save(1,1,inode) = dpox 
          dpoxyz_save(2,1,inode) = dpoy 
      end do

      ! assign boundary value has potentials--------
      ! bk is the potential from time stepping----on surface node
      do  inode=1,  nnf  
          do is=1, nsys; do ip=1, nsys
              cmat(inode,is)=cmat(inode,is)+rsn(is,ip)*bkn(inode,ip)
          enddo;enddo

          !      do j=1,ncn(nodele(inode,1))
          !tmp(j) = bkn(ncon(nodele(inode,1),j),1)!get surface value
          !end do
          !call eval_gradient(inode,tmp,tmp2)
          !dpoxyz_save(1,1,inode) = tmp2(1,1)
          !dpoxyz_save(2,1,inode) = tmp2(2,1)
      end do

      1021 format(8f10.6)

      ! assign boudary value has dp/dn ---------------:
      do 40 inode=nnf+1, nnoded
          do 40 ip=1, nsys 

              !FIXME inode is wrong
              n2 = nrml_2_node(inode)
              xp=ex(ip)*xyz(1,n2)
              yp=ey(ip)*xyz(2,n2)
              zp=       xyz(3,n2)
              call dinp(xp,yp,zp,dpox,dpoy,dpoz)   !get initial condition    
              !================================================================
              !timerk at each rugga kutta time step
              !call dpoxyz(h,g,ampn,phi_w,beta,wkn,freq,timerk,rampf,xp,yp,zp,&
              !&              nfreq,nwave,iorder,dpox,dpoy,dpoz)

              dpdn=(dpox*ex(ip)*dxyz(1,inode)+&
                  &           dpoy*ey(ip)*dxyz(2,inode)+&
                  &           dpoz       *dxyz(3,inode) )
              !       
              !         WRITE(102,*) 'DPDN=',DPDN,DPOX,DPOY,DPOZ
              ! consider body velocity contribution to dpdn,commented for
              ! testing 
              !dpdn=dpdn-dsdt(1)*ex(ip)*dxyz(1,inode)
              !dpdn=dpdn-dsdt(2)*ey(ip)*dxyz(2,inode)
              !dpdn=dpdn-dsdt(3)*       dxyz(3,inode)
              !dpdn=dpdn-dsdt(4)*ey(ip)*dxyz(4,inode)
              !dpdn=dpdn-dsdt(5)*ex(ip)*dxyz(5,inode)
              !dpdn=dpdn-dsdt(6)*ex(ip)*ey(ip)*dxyz(6,inode)
              !

              ! cmat[nff+1:nnoded] store dpdn
              do  is=1, nsys 
                  cmat(inode,is)=cmat(inode,is)+rsn(is,ip)*dpdn
              enddo

              !this is only incident wave, missing the diffraction wave
        40      continue   


        ! ==============================================
        !
        ! Generate bmata matrix---------
        !|||----------------------------------

        bmata(:,:)=0.0d0
        do is=1, nsys; do ind=1,   nnode
            !---------------loop-body--------------------------
            bmata(ind,is)=dot_product(cmata(ind,:,is),cmat(:,is))

            if (ind<=nnf) then!potential only
                bmata(ind,is) = bmata(ind,is)-fterm(ind,is,1)*cmat(ind,is)&
                    &-fterm(ind,is,2)*dpoxyz_save(1,is,ind) &
                    &-fterm(ind,is,3)*dpoxyz_save(2,is,ind)
            end if
            !---------------------------------------------------
        enddo;enddo 

        do  is=1, nsys   
            call rlubksb(is,amata,nnode,nnode, 1,nsys, 1,indx,bmata)
        enddo


        !do i = 1,nnode
        !write(401,*) bmata(i,1:nsys)
        !end do
        !do i=1,nnoded
        !!print *,cmat(i,1)
        !write(403,*) cmat(i,1)
        !enddo

        do ip=1, nsys; do ind=1, nnode 
            bmat(ip)=(0.0d0, 0.0d0)
            do  is=1, nsys
                bmat(ip)=bmat(ip)+bmata(ind,is)*rsn(ip,is)
            enddo
            bmat(ip)=bmat(ip)/nsys
            unkn(ind,ip)=bmat(ip)        
        enddo;enddo

    end subroutine
