!  TASSBT+ TINBODT+ NSBDT
!
!C *********************************************************
!C *                                                       *
!C * Calculate the element contribution and assembly the   *
!C * coefficients of the corresponding system of equationn *
!C *                                                       *
!C *********************************************************
!C
        SUBROUTINE TASSBT
  
      USE MVAR_MOD
      USE PVAR_MOD
      use wave2,only:dpoxyz
      use mfunc_mod,only:rlubksb
      use gradient,only:eval_gradient

      implicit none  

      integer  inode,ind,is,ip,j
      real(8) xp,yp,zp,dpox,dpoy,dpoz,dpdn
      real(8) rsn(4,4),ex(4),ey(4)
      real(8) bmat(nsys),cmat(nnoded,nsys)
      real(8) :: dpoxyz_save(2,4,nnoded),tmp(8),tmp2(2,1)

      data rsn /1.,  1.,  1.,  1.,& 
     &            1., -1.,  1., -1.,&
     &            1.,  1., -1., -1.,& 
     &            1., -1., -1.,  1./ 
      data ex / 1.,  1., -1., -1./                                     
      data ey / 1., -1., -1.,  1./

       print *,"Entering tassbt"
!
!        WRITE(102,*)
!        WRITE(102,*)   'T=', TimeRK                 

!C
!C ***************************************
!C
      cmat(:,:)=0.0
      dpoxyz_save = 0.0
        ! assign boundary value has potentials--------
        do  20   inode=1,  nnf  
            do is=1, nsys
                do ip=1, nsys
                  cmat(inode,is)=cmat(inode,is)+rsn(is,ip)*bkn(inode,ip)
            enddo;enddo
                !call dpoxyz(h,g,ampn,phi_w,beta,wkn,freq,timerk,rampf,xp,yp,zp,&
     !&              nfreq,nwave,iorder,dpox,dpoy,dpoz)
          !dpoxyz_save(1,ip,inode) = dpox
        !dpoxyz_save(2,ip,inode) = dpoy
        !print *,inode,"surface list getting"
        !print *,"nodele",nodele(inode,1)
        !print *,ncon(nodele(inode,1),:)
        do j=1,ncn(nodele(inode,1))
             tmp(j) = bkn(ncon(nodele(inode,1),j),1)!get surface value
        end do
        !print *,"surface list finished"
                !write(*,1021) tmp
                call eval_gradient(inode,tmp,tmp2)
                dpoxyz_save(1,1,inode) = tmp2(1,1)
                dpoxyz_save(2,1,inode) = tmp2(2,1)
                !write(*,*) tmp2(1,1),tmp2(2,1)
                
                !print *,"finished free surface node ",inode
        
 
    20   continue
    1021 format(8f10.6)
    !print *,"cmat",cmat
!
!        WRITE(102,*)  'INODE,    CMAT(INODE,1) --20'
!        DO    INODE=1,  NNF
!        WRITE(102, 620)  INODE,CMAT(INODE,1)
!        ENDDO


! ==============================================
!
!         WRITE(102,*) 'DSDT=',DSDT

        ! assign boudary value has dp/dn ---------------:
        do 40 inode=nnf+1, nnoded
        do 40 ip=1, nsys 
        xp=ex(ip)*xyz(1,inode)
        yp=ey(ip)*xyz(2,inode)
        zp=       xyz(3,inode)
 !      call dinp(h,amp,beta,wk,w1,g,timerk,rampf,
     !1                   xp,yp,zp,dpox,dpoy,dpoz)
        !timerk at each rugga kutta time step
        call dpoxyz(h,g,ampn,phi_w,beta,wkn,freq,timerk,rampf,xp,yp,zp,&
     &              nfreq,nwave,iorder,dpox,dpoy,dpoz)

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

        ! cmat[nff+1:nnoded] store dpdn
        do  is=1, nsys 
            cmat(inode,is)=cmat(inode,is)+rsn(is,ip)*dpdn
        enddo

             !this is only incident wave, missing the diffraction wave
40      continue   
!

!        WRITE(102,*)  'INODE,    CMAT(INODE,1) --40'
!        DO    INODE=NNF+1,  NNODED
!        WRITE(102, 620)  INODE,CMAT(INODE,1)
!        ENDDO

! ==============================================
!
        ! Generate bmata matrix---------
        !|||----------------------------------

        bmata(:,:)=0.0d0
        do 200 is=1, nsys
            do 200 ind=1,   nnode
                do inode=1, nnoded
        !---------------loop-body--------------------------
          bmata(ind,is)=bmata(ind,is)+cmata(ind,inode,is)*cmat(inode,is)
                end do
                if (ind<=nnf) then!potential only
                        bmata(ind,is) = bmata(ind,is)-fra3(ind,is)*cmat(ind,is)&
                         &-frc31(ind,is)*dpoxyz_save(1,is,ind) &
                         &-frc32(ind,is)*dpoxyz_save(2,is,ind)
                end if
        !---------------------------------------------------
         200     continue   
        
       ! 
       !        WRITE(102,*)  'INODE,XP,YP,ZP,BMATA(INODE,1) '
!        DO    INODE=1,  NNODE
!         XP=XYZ(1,INODE)
!         YP=XYZ(2,INODE)
!         ZP=XYZ(3,INODE) 
!        WRITE(102, 620)  INODE,XP,YP,ZP,BMATA(INODE,1)
!        ENDDO 
!

          do 300 is=1, nsys   
            call rlubksb(is,amata,nnode,nnode, 1,nsys, 1,indx,bmata)
300       continue


!C                 
!C ** output the results, compute unkn[1:NNF] is dpdn,unkn[NNF+1:nnode]
! is potential
!C      
        do 500 ip=1, nsys 
        do 360 ind=1, nnode 
        bmat(ip)=(0.0d0, 0.0d0)
        do 350 is=1, nsys
350     bmat(ip)=bmat(ip)+bmata(ind,is)*rsn(ip,is)
        bmat(ip)=bmat(ip)/nsys
        unkn(ind,ip)=bmat(ip)        
360     continue
500     continue

!        WRITE(102,*) ' INODE,XP,YP,ZP,BMATA(INODE,1),UNKN(INODE,1' 
!        DO    INODE=1,  NNODE
!        XP=XYZ(1,INODE)
!        YP=XYZ(2,INODE)
!        ZP=XYZ(3,INODE) 
!        WRITE(102, 620)  INODE,XP,YP,ZP,BMATA(INODE,1),UNKN(INODE,1)
!        ENDDO
  
 620   FORMAT(1X,I4,5(1X,E13.6))                            

            
       RETURN
       END