      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:),nF(:,:)
      integer ipars(2),d(3)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:),V(:,:),beta(:),alpha(:)
      real *8, allocatable :: dV(:,:,:)
      real *8 dpars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: ffform(:,:,:),ffformex(:,:,:)
      real *8, allocatable :: ffforminv(:,:,:),ffformexinv(:,:,:)

      complex *16, allocatable :: rhs(:)
      real *8, allocatable :: sigma(:), pot(:),pot1(:),rrhs(:)
      real *8, allocatable :: errs(:),rrhs2(:),u(:),rrhs2t(:)

      real *8, allocatable :: wnear(:),rrhs1(:),nFt(:,:)
      real *8, allocatable :: targs(:,:),F(:,:),Ft(:,:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8, allocatable :: srcover(:,:),wover(:),H(:,:)
      real *8, allocatable :: xmat(:,:),xtmp(:,:),slmat(:,:)
      real *8, allocatable :: slapmat(:,:), s0mat(:,:),nsgbeta(:,:)
      real *8, allocatable :: pdis(:)

      real *8, allocatable :: rfds(:),sgalpha(:,:),sgbeta(:,:),gu(:,:)
      integer, allocatable :: ifds(:)
      complex *16, allocatable :: zfds(:)


      integer, allocatable :: col_ptr(:),row_ind(:)
      integer, allocatable :: ixyzso(:),novers(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      integer, allocatable :: seed(:),seed_old(:)


      real *8 thet,phi,eps_gmres,Wg,Wu,rnd,Hnorm,Fnorm,dipvec(3)
      complex * 16 zpars(3)
      integer numit,niter,sizer,n_min,n_max,ran_len,ran_int,clock
      character *100 title,dirname
      character *300 fname,fname1,fname2,fname3,fname4,fname5

      real *8, allocatable :: w(:,:)

      logical isout0,isout1

      complex *16 ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      call random_seed(sizer)
      print *, sizer
      allocate(seed(sizer+1),seed_old(sizer+1))
 
!  illustrate use of put and get

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, sizer) /)
      CALL RANDOM_SEED(PUT = seed)



!  confirm value of seed
      print *, "seed = ", seed
      ran_len = 3         ! length of sequence
      n_min = 1
      n_max = 5
      do i = 1,ran_len
        call random_number(rnd)
        d(i) = (n_max - n_min + 1)*rnd + n_min
        print *,d(i)
      end do


         
            
      zk = 4.4d0+ima*0.0d0
      zpars(1) = zk 
      zpars(2) = -ima*zk
      zpars(3) = 2.0d0

      
      xyz_in(1) = 0.11d0
      xyz_in(2) = 0.0d-5
      xyz_in(3) = 0.37d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      igeomtype = 4

      if(igeomtype.eq.1) then
        ipars(1) = 3
        npatches=12*(4**ipars(1))
        fname='sphere.vtk'
      endif

      if(igeomtype.eq.2) then
        ipars(1) = 40
        ipars(2) = ipars(1)*3
        npatches = 2*ipars(1)*ipars(2)
        fname='stell.vtk'
      endif

      if(igeomtype.eq.3) then
        ipars(1) = 40
        ipars(2) = 20
        npatches = 2*ipars(1)*ipars(2)
        fname='wtorus.vtk'
      endif

      if(igeomtype.eq.4) then
        ipars(1) = 4*8
        ipars(2) = 4*8
        npatches = 2*ipars(1)*ipars(2)
        fname = 'torus.vtk'
        xyz_out(1) = 0.1d0
        xyz_out(2) = 0.2d0
        xyz_out(3) = 2.1d0
      endif

      dipvec(1) = 0.37d0
      dipvec(2) = 0.48d0
      dipvec(3) = -0.8d0
      

      norder = 8 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
      allocate(wts(npts))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      print *, 'npts=',npts

      ixyzs(npatches+1) = 1+npols*npatches



c      fname = '../../geometries/Genus_10_files/Genus_10_o08_r03.go3'
c      
c      call open_gov3_geometry_mem(fname,npatches,npts)
c
c      call prinf('npatches=*',npatches,1)
c      call prinf('npts=*',npts,1)
c
c      allocate(srcvals(12,npts),srccoefs(9,npts))
c      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
c      allocate(wts(npts))
c
c      call open_gov3_geometry(fname,npatches,norders,ixyzs,
c     1   iptype,npts,srcvals,srccoefs,wts)

       fname1='tmp1.vtk'
       fname2='tmp2.vtk'
       fname3='tmp3.vtk'
       fname4='tmp4.vtk'
       fname5='tmp5.vtk'
c      fname1='../../../lbres/res/g10f/hodge_tanproj_g10_o08_r03_'
c     1             //CHAR(d(1)+48)//
c     1        '_'//CHAR(d(2)+48)//'_'//CHAR(d(3)+48)//'.vtk'
c      fname2='../../../lbres/res/g10f/hodge_dfree_g10_o08_r03_'
c     1            //CHAR(d(1)+48)//
c     1        '_'//CHAR(d(2)+48)//'_'//CHAR(d(3)+48)//'.vtk'
c
c      fname3='../../../lbres/res/g10f/hodge_cfree_g10_o08_r03_'
c     1            //CHAR(d(1)+48)//
c     1        '_'//CHAR(d(2)+48)//'_'//CHAR(d(3)+48)//'.vtk'
c
c      fname4='../../../lbres/res/g10f/hodge_harm_g10_o08_r03_'
c     1               //CHAR(d(1)+48)//
c     1        '_'//CHAR(d(2)+48)//'_'//CHAR(d(3)+48)//'.vtk'
c
c      fname5='../../../lbres/res/g10f/hodge_geo_g10_o08_r03_'
c     1             //CHAR(d(1)+48)//
c     1        '_'//CHAR(d(2)+48)//'_'//CHAR(d(3)+48)//'.vtk'


      print *, 'F1:', fname1
      print *, 'F2:', fname2
      print *, 'F3:', fname3
      print *, 'F4:', fname4
      print *, 'F5:', fname5

      
      allocate(rrhs2t(npts),H(3,npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      allocate(pdis(npatches))
      do i=1,npatches
        pdis(i) = 0
      enddo
      call get_patch_distortion(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,wts,pdis)
      pmax = maxval(pdis)
      call prin2('max patch distortion=*',pmax,1)
      pavg = 0
      do i=1,npatches
        pavg = pavg+log(pdis(i))
      enddo
      pavg = pavg/npatches
      pavg = exp(pavg)
      call prin2('geometric mean of patches=*',pavg,1)



      allocate(sigma(npts),rhs(npts),pot(npts),rrhs(npts))
      allocate(ffform(2,2,npts),rrhs2(npts),u(npts))
      allocate(rrhs1(npts),Ft(3,npts),nFt(3,npts))
      allocate(V(3,npts),F(3,npts),nF(3,npts))
      allocate(alpha(npts),beta(npts),nsgbeta(3,npts))
      allocate(sgalpha(3,npts),sgbeta(3,npts))
      allocate(gu(3,npts))
c
c       define rhs to be one of the ynm's
c
      nn = 2
      mm = 1
      nmax = nn
      allocate(w(0:nmax,0:nmax))
      call l3getsph(nmax,mm,nn,12,srcvals,rhs,npts,w)
      
      do i=1,npts
        rrhs(i) = real(rhs(i))
      enddo

      eps = 0.51d-10


c
c       precompute near quadrature correction
c
c
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO     

      ntarg = npts
      allocate(targs(3,npts))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO      

c
c    find near quadrature correction interactions
c
      call findnearmem(cms,npatches,rad_near,3,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,3,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,
     1         iquad)

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)
      


c
c    estimate oversampling for far-field, and oversample geometry
c

      ikerorder = 0
      allocate(novers(npatches),ixyzso(npatches+1))

      zpars = 0
      ndtarg = 3
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      allocate(wnear(4*nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,4*nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    

      call prinf('finished generating near field info*',i,0)
      call prinf('finished generating far field orders*',i,0)
      call prinf('npts_over=*',npts_over,1)
      call prin2('eps=*',eps,1)

      

      iquadtype = 1

      call getnearquad_lap_bel2fast(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
      call prinf('finished generating near quadrature correction*',i,0)


      
c     Evaluate DL on surface to check geometry info
      dpars(2) = 1.0d0
      dpars(1) = 0
 
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        sigma(i) = 1
      enddo
C$OMP END PARALLEL DO      


      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2  dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(1*nquad+1),
     3  sigma,novers,npts_over,ixyzso,srcover,wover,beta)
 
      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=  erra + ((beta(i)+0.5)**2)*wts(i)
        ra= ra + (0.5)**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in DL =*',erra,1)

      



c      call prin2('pot=*',pot,24)
c      call prin2('rrhs=*',rrhs,24)
     
      call prin2('starting FAST iterative solve*',i,0)
      numit = 50
      allocate(errs(numit+1))
     
      Wg = 0 
      do i=1,npts
        V(1,i) = 1.0d0
        V(2,i) = 0.0d0
        V(3,i) = 0.0d0
        Ft(1,i)=1.0d0-srcvals(1,i)**2
        Ft(2,i)=-1.0d0*srcvals(1,i)*srcvals(2,i)
        Ft(3,i)=-1.0d0*srcvals(1,i)*srcvals(3,i) 
        nFt(1,i)=0
        nFt(2,i)=srcvals(3,i)
        nFt(3,i)=-1.0d0*srcvals(2,i)

        rrhs2t(i) = -2.0d0*srcvals(1,i)
        V(1,i) = (srcvals(1,i))**2
        V(2,i) = (5*srcvals(1,i)**2+6*srcvals(2,i)**3)
        V(3,i) = (srcvals(3,i)**4)
        Ft(1,i)=-srcvals(1,i)*srcvals(3,i)**2
        Ft(2,i)=-srcvals(2,i)*srcvals(3,i)**2
        Ft(3,i)=srcvals(3,i)-srcvals(3,i)**3 
        nFt(1,i)=srcvals(2,i)*srcvals(3,i)
        nFt(2,i)=-srcvals(1,i)*srcvals(3,i)
        nFt(3,i)=0.0d0
        rrhs2t(i)=1.0d0-3*srcvals(3,i)**2    

        u(i) = srcvals(1,i)**2+srcvals(2,i)**2
        gu(1,i)=2.0d0*srcvals(1,i)*(srcvals(3,i)**2)
        gu(2,i)=2.0d0*srcvals(2,i)*(srcvals(3,i)**2)
        gu(3,i)=2.0d0*srcvals(3,i)*(srcvals(3,i)**2-1.0d0)
        sigma(i) = 0 
        alpha(i) = 0
        beta(i) = 0
      enddo

      allocate(dV(3,3,npts))
 
      call biot_savart(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,xyz_out,dipvec,V,dV)


c      call poly_field(npatches,norders,ixyzs,iptype,npts,
c     1   srccoefs,srcvals,d,V)


      title = 'u(x)'
      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,u,fname5,title)

      do i=1,npts
        Wg = Wg + (1.0d0)**2*wts(i) 
      enddo
      Wg = (Wg )
      call prin2('area of surface=*',Wg,1)

 

      call surf_grad(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,u,sgalpha)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=erra+((sgalpha(1,i))**2)*wts(i)
        erra=erra+((sgalpha(2,i))**2)*wts(i)
        erra=erra+((sgalpha(3,i))**2)*wts(i) 
        ra=ra+((gu(1,i))**2)*wts(i)
        ra=ra+((gu(2,i))**2)*wts(i)
        ra=ra+((gu(3,i))**2)*wts(i) 


      enddo
      erra = sqrt(erra)
      call prin2('norm of surface grad=*',erra,1)


      print *, 'calling tangential projection'
   

      call tangential_projection(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,V,F)

      print *, 'F calculated'

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=erra+((F(1,i))**2)*wts(i)
        erra=erra+((F(2,i))**2)*wts(i)
        erra=erra+((F(3,i))**2)*wts(i) 

      enddo
      erra = sqrt(erra)
      Fnorm = erra
      rr = 1.0d0/Fnorm !scale factor
      call prin2('unscaled L2 norm of tan proj F=*',erra,1)

c     scale F to norm 1
      do i=1,npts
        F(1,i)=rr*F(1,i)
        F(2,i)=rr*F(2,i)
        F(3,i)=rr*F(3,i) 
      enddo

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=erra+((F(1,i))**2)*wts(i)
        erra=erra+((F(2,i))**2)*wts(i)
        erra=erra+((F(3,i))**2)*wts(i) 

      enddo
      erra = sqrt(erra)
      Fnorm = erra
      call prin2('scaled L2 norm of tan proj F=*',erra,1)



      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,F,fname1,'F') 

       

      call ncross(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,F,nF)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=  erra + ((nF(1,i)-nFt(1,i))**2)*wts(i)
        erra=  erra + ((nF(2,i)-nFt(2,i))**2)*wts(i)
        erra=  erra + ((nF(3,i)-nFt(3,i))**2)*wts(i) 
        ra = ra + (nFt(1,i))**2*wts(i)
        ra = ra + (nFt(2,i))**2*wts(i)
        ra = ra + (nFt(3,i))**2*wts(i)
 
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in nxF =*',erra,1)

      


      call surf_div(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,F,rrhs1)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=  erra + ((rrhs1(i)-rrhs2t(i))**2)*wts(i)
        ra = ra + (rrhs2t(i))**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in rrhs1 =*',erra,1)

      

      call surf_div(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,nF,rrhs2)

       do i=1,npts
         rrhs2(i) = -1.0d0*rrhs2(i) 
      enddo

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=  erra + ((rrhs2(i))**2)*wts(i)
      enddo
      call prin2('error in rrhs2 =*',erra,1)



 
c    Surface integral should be zero
      Wu = 0 
      do i=1,npts
        Wu = Wu + rrhs1(i)*wts(i) 
      enddo
      call prin2('surface integral of rrhs1=*',Wu,1)


 
c    Surface integral should be zero
      Wu = 0 
      do i=1,npts
        Wu = Wu + rrhs2(i)*wts(i) 
      enddo
      call prin2('surface integral of rrhs2=*',Wu,1)

       

      eps_gmres = 1.0d-14
      call lap_bel_solver2fast(npatches,norders,ixyzs,iptype,
     1        npts,srccoefs,
     1  srcvals,eps,numit,rrhs1,eps_gmres,niter,errs,rres,sigma) 

      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)

      dpars(1) = 1.0d0
      dpars(2) = 0
 

      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2  dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(0*nquad+1),
     3  sigma,novers,npts_over,ixyzso,srcover,wover,alpha)
 
c    Surface integral should be zero
      Wu = 0 
      do i=1,npts
        Wu = Wu + alpha(i)*wts(i) 
      enddo
      call prin2('surface integral of alpha=*',Wu,1)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=erra+((alpha(i))**2)*wts(i)
      enddo
      erra = sqrt(erra)
      call prin2('norm of surface alpha=*',erra,1)



      eps_gmres = 1.0d-14
      call lap_bel_solver2fast(npatches,norders,ixyzs,iptype,
     1      npts,srccoefs,
     1  srcvals,eps,numit,rrhs2,eps_gmres,niter,errs,rres,sigma) 

      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)

      dpars(1) = 1.0d0
      dpars(2) = 0
 

      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2  dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(0*nquad+1),
     3  sigma,novers,npts_over,ixyzso,srcover,wover,beta)
 
c    Surface integral should be zero
      Wu = 0 
      do i=1,npts
        Wu = Wu + beta(i)*wts(i) 
      enddo
      call prin2('surface integral of beta=*',Wu,1)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=erra+((beta(i))**2)*wts(i)
      enddo
      erra = sqrt(erra)
      call prin2('norm of surface beta=*',erra,1)




      call surf_grad(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,alpha,sgalpha)


      call surf_grad(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,beta,sgbeta)


      call ncross(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,sgbeta,nsgbeta)

      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,sgalpha,fname2,
     1  'sg_alpha') 

      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,nsgbeta,fname3,'nsg_beta') 



C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        H(1,i) = F(1,i)-sgalpha(1,i)-nsgbeta(1,i)
        H(2,i) = F(2,i)-sgalpha(2,i)-nsgbeta(2,i)
        H(3,i) = F(3,i)-sgalpha(3,i)-nsgbeta(3,i)
      enddo
C$OMP END PARALLEL DO      

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=erra+((H(1,i))**2)*wts(i)
        erra=erra+((H(2,i))**2)*wts(i)
        erra=erra+((H(3,i))**2)*wts(i) 

      enddo
      erra = sqrt(erra)
      Hnorm = erra
      call prin2('L2 norm of harmonic H=*',erra,1)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=erra+((H(1,i))**2)
        erra=erra+((H(2,i))**2)
        erra=erra+((H(3,i))**2) 
        erra = sqrt(erra)
        ra = max(ra, erra)
        erra = 0.0d0
      enddo
      call prin2('L_inf norm of harmonic H=*',ra,1)

      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,H,fname4,'H') 



      call surf_div(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,H,rrhs1)





      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=  erra + ((rrhs1(i))**2)*wts(i)
      enddo
      erra = sqrt(erra)
      call prin2('harmonic comp error in div =*',erra,1)
      call prin2('harmonic comp error in div rel H=*',erra/Hnorm,1)
      call prin2('harmonic comp error in div rel F=*',erra/Fnorm,1)



      call ncross(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,H,nsgbeta)


      call surf_div(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,nsgbeta,rrhs2)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=  erra + ((rrhs2(i))**2)*wts(i)
      enddo
      erra = sqrt(erra)
      call prin2('harmonic comp error in curl =*',erra,1)
      call prin2('harmonic comp error in curl rel H=*',erra/Hnorm,1)
      call prin2('harmonic comp error in curl rel F=*',erra/Fnorm,1)



      
      stop
      end





      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval,xtri_wtorus_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.3) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)
        call prinf('npatches=*',npatches,1)
         
        p1(1) = 1
        p1(2) = 2
        p1(3) = 0.25d0

        p2(1) = 1.2d0
        p2(2) = 1.0d0
        p2(3) = 1.7d0

c
c         numberof oscillations
c
        p4(1) = 5.0d0


        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_eval
        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,
     2         'Triangulated surface of the wtorus')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      if(igeomtype.eq.4) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)
        call prinf('npatches=*',npatches,1)
         
        p1(1) = 1.0d0
        p1(2) = 3.0d0
        p1(3) = 0.25d0

        p2(1) = 1.0d0
        p2(2) = 1.0d0
        p2(3) = 1.0d0

c
c         number of oscillations
c
        p4(1) = 0.0d0


        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_eval
        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,
     2         'Triangulated surface of the torus')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end


      subroutine test_exterior_pt(npatches,norder,npts,srcvals,
     1   srccoefs,wts,xyzout,isout)
c
c
c  this subroutine tests whether the pt xyzin, is
c  in the exterior of a surface, and also estimates the error
c  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
c  centered at the interior point. Whether a point 
c  is in the interior or not is tested using Gauss' 
c  identity for the flux due to a point charge
c
c
c  input:
c    npatches - integer
c       number of patches
c    norder - integer
c       order of discretization
c    npts - integer
c       total number of discretization points on the surface
c    srccoefs - real *8 (9,npts)
c       koornwinder expansion coefficients of geometry info
c    xyzout -  real *8 (3)
c       point to be tested
c
c  output: 
c    isout - boolean
c      whether the target is in the interior or not
c

      implicit none
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2


      zk = 0

      ra = 0



      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

      return
      end

   



      subroutine l3getsph(nmax,mm,nn,ndx,xyzs,ynms,npts,ynm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(npts),ima
      real *8 rat1(10000),rat2(10000)
      real *8 ynm(0:nmax,0:nmax)
      data ima/(0.0d0,1.0d0)/
  
      call ylgndrini(nmax, rat1, rat2)
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        r=sqrt(x**2+y**2+z**2)
        call cart2polar(xyzs(1,i),r,theta,phi)
        ctheta = cos(theta)
        call ylgndrf(nmax, ctheta, ynm, rat1, rat2)
        ynms(i) = ynm(nn,abs(mm))*exp(ima*mm*phi)        
      enddo
       
      return
      end


      subroutine zero(npts,F)
      implicit none
      integer  npts
      real *8 :: F(3,npts)
      integer i

      do i=1,npts
        F(1,i) = 0 
        F(2,i) = 0 
        F(3,i) = 0 
      enddo

      return
      end


      subroutine ncross(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,F,nF) 
      implicit none
      integer  npatches,norders(npatches)
      integer  ixyzs(npatches+1),iptype(npatches)
      integer   npts
      real *8  srccoefs(9,npts),srcvals(12,npts)
      real *8   F(3,npts)
      real *8   nF(3,npts)
      integer i

      do i=1,npts
        nF(1,i) = srcvals(11,i)*F(3,i)-srcvals(12,i)*F(2,i) 
        nF(2,i) = srcvals(12,i)*F(1,i)-srcvals(10,i)*F(3,i) 
        nF(3,i) = srcvals(10,i)*F(2,i)-srcvals(11,i)*F(1,i) 
      enddo

      return
      end


      subroutine tangential_projection(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,V,F) 
      implicit none
      integer  npatches,norders(npatches)
      integer  ixyzs(npatches+1),iptype(npatches)
      integer  npts
      real *8  srccoefs(9,npts),srcvals(12,npts)
      real *8  V(3,npts)
      real *8  F(3,npts)
      integer i
      real *8 ndotV

      ndotV = 0
      do i=1,npts
        ndotV = srcvals(10,i)*V(1,i) + srcvals(11,i)*V(2,i) 
     1         + srcvals(12,i)*V(3,i)

        F(1,i) = V(1,i)-srcvals(10,i)*ndotV
        F(2,i) = V(2,i)-srcvals(11,i)*ndotV
        F(3,i) = V(3,i)-srcvals(12,i)*ndotV 
      enddo

      return
      end


      subroutine biot_savart(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,xyz0,dipvec,V,dV) 
      implicit none
      integer  npatches,norders(npatches)
      integer  ixyzs(npatches+1),iptype(npatches)
      integer  npts
      real *8  srccoefs(9,npts),srcvals(12,npts)
      real *8  V(3,npts),xyz0(3),dV(3,3,npts)
      integer i,j,l
      real *8 lcross(3,npts),dipvec(3),nrm
      real *8 vtmp1(3),vtmp2(3),dxyz(3),r,wtmp1(3),ru,rv
      real *8 rinv3,rinv5

      do i=1,npts
        dxyz(1) = srcvals(1,i) - xyz0(1)
        dxyz(2) = srcvals(2,i) - xyz0(2)
        dxyz(3) = srcvals(3,i) - xyz0(3)
        r  = sqrt(dxyz(1)**2 + dxyz(2)**2 + dxyz(3)**2)
        rinv3 = 1.0d0/r**3
        rinv5 = 1.0d0/r**5
        wtmp1(1:3) = dxyz(1:3)/r**3
        call cross_prod3d(dipvec,wtmp1,V(1,i))
        do j=1,3
          do l=1,3
            dV(j,l,i) = -3*dxyz(j)*dxyz(l)*rinv5/2
          enddo
          dV(j,j,i) = dV(j,j,i) + 1.0d0/2*rinv3
        enddo
      enddo

cc      do i=1,npts
cc        lcross(1,i) = l(2)*(srcvals(3,i)-xyz0(3))-
cc     1     l(3)*(srcvals(2,i)-xyz0(2)) 
cc        lcross(2,i) = l(3)*(srcvals(1,i)-xyz0(1))-
cc     1      l(1)*(srcvals(3,i)-xyz0(3)) 
cc        lcross(3,i) = l(1)*(srcvals(2,i)-xyz0(2))-
cc     1     l(2)*(srcvals(1,i)-xyz0(1)) 
cc      enddo
cc
cc
cc      do i=1,npts
cc        nrm = (srcvals(1,i)-xyz0(1))**2+(srcvals(2,i)-xyz0(2))**2
cc     1           +(srcvals(3,i)-xyz0(3))**2
cc        nrm = nrm**(1.5)
cc        V(1,i) = lcross(1,i)/nrm
cc        V(2,i) = lcross(2,i)/nrm
cc        V(3,i) = lcross(3,i)/nrm 
cc      enddo

      return
      end


      subroutine poly_field(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,d,V) 
      implicit none
      integer  npatches,norders(npatches)
      integer  ixyzs(npatches+1),iptype(npatches)
      integer  npts,d(3)
      real *8  srccoefs(9,npts),srcvals(12,npts)
      real *8  V(3,npts)
      integer i,npoly1,npoly2,npoly3
      real *8 lcross(3,npts),l(3),nrm
c      real *8 val11(npts,0:d(1)),val12(npts,0:d(2)),val13(npts,0:d(3))
c      real *8 val21(npts,0:d(1)),val22(npts,0:d(2)),val23(npts,0:d(3))
c      real *8 val31(npts,0:d(1)),val32(npts,0:d(2)),val33(npts,0:d(3))
      real *8, allocatable :: pols1(:),pols2(:),pols3(:)
     
      print *, 'Points assigned polyfield'
      npoly1 = (d(1)+1)*(d(1)+2)*(d(1)+3)/6
      npoly2 = (d(2)+1)*(d(2)+2)*(d(2)+3)/6
      npoly3 = (d(3)+1)*(d(3)+2)*(d(3)+3)/6
 
      allocate(pols1(npoly1),pols2(npoly2),pols3(npoly3))
          
c      call j_polynomial(npts,d(1),0,0,xx,val11)
c      call prin2('legen poly =*',val11,24)
c
c      call j_polynomial(npts,d(1),0,0,xx,val11)
c      print *, 'first j poly called'
c      call j_polynomial(npts,d(1),0,0,yy,val12)
c      call j_polynomial(npts,d(1),0,0,zz,val13)
c      call j_polynomial(npts,d(2),0,0,xx,val21)
c      call j_polynomial(npts,d(2),0,0,yy,val22)
c      call j_polynomial(npts,d(2),0,0,zz,val23)
c      call j_polynomial(npts,d(3),0,0,xx,val31)
c      call j_polynomial(npts,d(3),0,0,yy,val32)
c      call j_polynomial(npts,d(3),0,0,zz,val33)
 
      print *,'all j polys called'
      do i=1,npts
c        V(1,i) = val11(i,d(1))*val12(i,d(1))*val13(i,d(1))
c        V(2,i) = val21(i,d(2))*val22(i,d(2))*val23(i,d(2))
c        V(3,i) = val31(i,d(3))*val32(i,d(3))*val33(i,d(3))
c        V(1,i) = 1.0d0
c        V(2,i) = 1.0d0
c        V(3,i) = 1.0d0
        l(1) = 0.7*(0.2*srcvals(1,i) -1.0d0)
        l(2) = 0.7*(1.0d0/7)*srcvals(2,i)
        l(3) = 0.7*(2*srcvals(3,i) - 1.0d0)
        call legetens_pols_3d(l,d(1),'T',pols1)
        call legetens_pols_3d(l,d(2),'T',pols2)
        call legetens_pols_3d(l,d(3),'T',pols3)
        V(1,i) = pols1(npoly1/2)
        V(2,i) = pols2(npoly2/2)
        V(3,i) = pols3(npoly3/2)
      enddo
      print *, 'V assigned'
      return
      end


