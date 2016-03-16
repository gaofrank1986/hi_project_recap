program
end program
    subroutine eval_regular_elem()
        implicit none


        !call choose_gaussian_pts()
        if(sum(nsub(1:3).eq.3) then
            do igl=1,ngsp(2)
                !if(ndim.eq.3) then
                ksi(2)=gpl(igl)
                wl=gwl(igl)
                !else
                    !ksi(2) = 0
                    !wl= 1
                !endif
                do igr=1,ngsp(1)
                    ksi(1)=gpr(1)
                    wr=gwr(igr)
                    call integrate_elem(beta,cp,ck,xi,wr*wl,1.0d0,cdl,nf,v1e)
            enddo;enddo
            return
        endif

        !ncor=2**2
        !do id=1,ncor
            !csub(2,id)=-1.+(cdl(2*id)+1.)/nsub(2)
        !enddo
        !do 70 isub3=1,nsub(2)
            !do id=1,ncor
                !csub(1,id)=-1.+(cdl(2*id-1)+1.)/nsub(1)
            !enddo
            !do 60 isub2=1,nsub(1)
                !call mini_distance
                !call shapef
                !call shapef
                !call choose_gaussian_pts
                !do 40 igl=1,ngsp(2)
                    !if(ndim.eq.3) then
                        !xic(2)=gpl(igl)
                        !wl=gwl(igl)
                    !else
                        !xic(2)=0.
                        !wl=1.
                    !endif
                    !do 40 igr=1,ngsp(1)
                    !xic(1)=gpr(igr)
                    !wr=gwr(igr)
                    !call shapef
                    !call dsahpe
                    !call integrate_elem
            !40  continue
        !60  csub(1,1:ncor)=csub(1,1:ncor)+dgs(1)
    !70  csub(2,1:ncor)=csub(2,1:ncor)+dgs(2)
    end subroutine

    subroutine integrate_elem(beta,xp,ck,xi,wxy,fjcbl,cdl,nf,v1e)        
        !beta order of r on denominator
        !xi local position
        !xp : src point global pos
        !ck : global corner matrix
        !cdl :local corner matrix
        !nf :number of f_bar component
        !fcbl :was set to 1

        call shapef(3,8,ksi,ck,xp,ri,cdl,shap)
        call dshape(3,2,8,ksi,ck,cosn,fjcb,cdl,gd)
        rq=norm2(ri)
        xq=xp+ri!xq is local point's global position
        drdx=ri/rq!ri is normalized
        drdn=dot_product(cosn,drdx)
        !fjcb is jacobian determinant
        comt=fjcb*wxy*fjcbl/rq**beta
        call f_bar
        v1e=v1e+comt*fq
    end subroutine

                    

!    subroutine choose_gaussian_pts(3,2,node,cp,ck,cdl,ngr,ngl,tolgp,beta,dgs,xic,al,disl,1,nsub,ngsp,gpr,gwr,gpl,gwl)
        !implicit none
        !integer :: nsub(3),ngsp(3)

        !ng(1)=ngr
        !ng(2)=ngl
        !pb=dsqrt(order*2./3.+0.4)
        
        !if(nflg.eq.0) then
            !!goto 30
            !print *,"choose gussian pts nflg is set to 0"
            !pause
        !endif
        
        !nsub(1:2)=1!number of sub along each direction
        !ngsp(1:2)=iabs(ng(1:2))
        
        !do i=2+1,3
            !dgs(i)=2
            !nsub(i)=1!direction 3
            !ngsp(i)=1!number of gaussion point along i
        !enddo
        
        !if(ng(1).gt.0.and.ng(2).gt.0) goto 90
        !! -10 and -10

        !call setup_gaussian(8,1,gpr,gwr,gpl,gwl,total_gpts)
        !! 8 in x,1 in y

        !avl = 0.
        
        !do 20 isid = 1,2!ndimb changing index 1=ksi,2=eta
            !al(isid) = 0
            !do i=1,2!ndimb
                !if(i.ne.isid) xi(i) = 0.0d0
            !enddo
            !do 10 ig=1,total_gpts!8 points
                !xi(isid) = gpr(ig)! 8 points 
                !call dshape---
                !fjcb=norm2(gd(1:3,isid)
                !al(isid)=al(isid)+fjcb*gwr(ig)
      !10    continue     
        !avl=avl+al(isid)/num_dim
   !20   continue
        !call mini_distance()

    !30  wfa=pb*dlog(dabs(tolgp)/2.)
        !do 80 i=1,num_dim - 1
            !if(ng(i) > 0) goto 80
            !if(tolgp<0)  goto 40
            !alm=al(i)
            !if (alm>3.9*disl) alm=3.9*disl
            !ngsp(i) = 0.5*wfa/dlog(alm/disl/4.)
            !goto 50
    !40      ngsp(i)=-wfa/10.*(8./3.*al(i)/disl)**(3.d0/4.)+1.)
    !50      if(ngsp(i)<2) ngsp=2
            !if((nflag.eq.1).and.(ng(i)<0) goto 60
            !if(****** goto 80
    !60      if  *****
            !if(tolgp<0) ali=3./8/*disl(-10*ngsp(i)/wfa-1.)**(4.d0/3)
            !if(tolgp>0) li=4.*disl*(tolgp/2.)**(0.5*pb/ngsp(i))
            !nsub(i)=al(i)/ali+0.95d0 
            !al(i)=al(i)/nsub(i)
    !70      dgs(i)=2.0d0/nsub(i)
    !80  continue
    !90  call setup_gaussian(ngsp(1),ngsp(2),gpr,gwr,gpl,gwl,ngss)
    !end subroutine choose_gaussian_pts

    subroutine setup_gaussian(n1,n2,pts_x,weigt_x,pts_y,weight_y,total_gpts)
        implicit none
        call gaussv(n1,pts_x,weight_x)
        if(n2.ne.1) call gaussv(n2,pts_y,weight_y)
        total_gpts=n1*n2
    end subroutine



!    subroutine mini_distance()
        !!xic given local
        !!ksi is searched iteratively
        !told=1.0d-11
        !tolb=told*10
        !num_iter=50
        !ksi(1:2) = 0.0d0
        !disl=1.0d12
        !do 60 iter=1,num+iter
            !call dshape()
            !if(dabs(fjcb)>told) then 
                !ksi(1:2)=0.382d0*xic(1:2)++0.618d0*ksi(1:2)
                !call dshape
            !endif
            !call shapef()!dist from ksi,to cordl
            !rm=norm2(ri)
            !if(dabs(rm-disl)/avl < tolb) return
            !! if rm converged exit
            !if(rm<disl) then
                !disl = rm
                !xic(1:2) =ksi(1:2)
            !else
                !ksi(1:2) = 0.382d0*xic(1:2)+0.618d0*ksi(1:2)
                !!golden ratio?
                !goto 60
            !endif
            !!if(2.eq.3) then
                !!goto 40
            !!endif
            !!if((3-1).eq.1) then
                !!goto 20
            !!endif
            !do 10 i=1,2
                !do 10 j=1,2
                    !gdt(i,j)=0
                    !do k=1,3
                        !gdt(i,j)=gdt(i,j)+gd(k,i)*gd(k,j)
                    !end do
        !10  continue
            !fjcb=gdt(1,1)*gdt(2,2)-gdt(1,2)*gdt(2,1)
            !call ivsnr123d()
        !20  do 30 i=1,2
                !do 30 j=1,3
                    !gdt(i,j)=0.0d
                    !do k=1,2
                        !gdt= rev*gd
                    !enddo
        !30  continue
        !! gdt(2,3)
        !40  do 50 i=1,2
            !term=0.0d0
            !do k=1,3
                !term=term+gdt(i,k)*ri(k)
            !enddo
        !50  ksi(i)=ksi(i)-term
            !do i=1,2
                !if(dabs(ksi(i))>1.0d0) ksi(i)=ksi(i)/dabs(ksi(i))
            !enddo
        !60 continue
    !end subroutine


