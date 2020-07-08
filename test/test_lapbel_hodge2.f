      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:),nF(:,:)
      integer ipars(2)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:),V(:,:),beta(:),alpha(:)
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

      real *8, allocatable :: rfds(:),sgalpha(:,:),sgbeta(:,:)
      integer, allocatable :: ifds(:)
      complex *16, allocatable :: zfds(:)


      integer, allocatable :: col_ptr(:),row_ind(:)
      integer, allocatable :: ixyzso(:),novers(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)



      real *8 thet,phi,eps_gmres,Wg,Wu
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname

      real *8, allocatable :: w(:,:)

      logical isout0,isout1

      complex *16 ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


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

      igeomtype = 1
      ipars(1) = 3
      npatches=12*(4**ipars(1))

      norder = 4 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
c      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0

c      call setup_geom(igeomtype,norder,npatches,ipars, 
c     1       srcvals,srccoefs,ifplot,fname)


      fname = '../../geometries/Genus_10_o04_r03.go3'
      
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,
     1   iptype,npts,srcvals,srccoefs,wts)
 

c      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

c      do i=1,npatches
c        norders(i) = norder
c        ixyzs(i) = 1 +(i-1)*npols
c        iptype(i) = 1
c      enddo

c      print *, 'npts=',npts

c      ixyzs(npatches+1) = 1+npols*npatches
      allocate(rrhs2t(npts),H(3,npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      allocate(sigma(npts),rhs(npts),pot(npts),rrhs(npts))
      allocate(ffform(2,2,npts),rrhs2(npts),u(npts))
      allocate(rrhs1(npts),Ft(3,npts),nFt(3,npts))
      allocate(V(3,npts),F(3,npts),nF(3,npts))
      allocate(alpha(npts),beta(npts),nsgbeta(3,npts))
      allocate(sgalpha(3,npts),sgbeta(3,npts))
      do i=1,npts
        u(i) = srcvals(1,i)
      enddo

 
      fname = 'g10_lapbel.vtk'
      title = 'f(x)'
      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,u,fname,title)

      

c
c       define rhs to be one of the ynm's
c
      nn = 4
      mm = 1
      nmax = nn
      allocate(w(0:nmax,0:nmax))
      call l3getsph(nmax,mm,nn,12,srcvals,rhs,npts,w)
      
      do i=1,npts
        rrhs(i) = real(rhs(i))
      enddo

      eps = 0.51d-8


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
      call findnearmem(cms,npatches,rad_near,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,targs,npts,row_ptr, 
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

      call getnearquad_lap_bel2(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
      call prinf('finished generating near quadrature correction*',i,0)


      call prinf('entering layer potential eval*',i,0)
      call prinf('npts=*',npts,1)

c      call lpcomp_lap_bel_addsub2(npatches,norders,ixyzs,iptype,npts,
c     1  srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,
c     2  rrhs,novers,npts_over,ixyzso,srcover,wover,pot)

c     Evaluate DL on surface to check geometry info
      dpars(2) = 1.0d0/4/pi
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
        erra=  erra + ((beta(i)+0.5d0)**2)*wts(i)
        ra= ra + (0.5)**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in DL =*',erra,1)

 

      call prin2('starting iterative solve*',i,0)
      numit = 50
      allocate(errs(numit+1))
     
      Wg = 0 
      do i=1,npts
        V(1,i) = srcvals(1,i)**2
        V(2,i) = srcvals(2,i)
        V(3,i) = srcvals(3,i)
        Ft(1,i)=1.0d0-srcvals(1,i)**2
        Ft(2,i)=-1.0d0*srcvals(1,i)*srcvals(2,i)
        Ft(3,i)=-1.0d0*srcvals(1,i)*srcvals(3,i) 
        nFt(1,i)=0
        nFt(2,i)=srcvals(3,i)
        nFt(3,i)=-1.0d0*srcvals(2,i)
        rrhs2t(i) = -2.0d0*srcvals(1,i)
        u(i) = srcvals(1,i)
        sigma(i) = 0 
        alpha(i) = 0
        beta(i) = 0
      enddo

    

      call tangential_projection(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,V,F)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=  erra + ((F(1,i)-Ft(1,i))**2)*wts(i)
        erra=  erra + ((F(2,i)-Ft(2,i))**2)*wts(i)
        erra=  erra + ((F(3,i)-Ft(3,i))**2)*wts(i) 
        ra = ra + (Ft(1,i))**2*wts(i)
        ra = ra + (Ft(2,i))**2*wts(i)
        ra = ra + (Ft(3,i))**2*wts(i)
 
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in tangential proj=*',erra,1)


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

       

200      eps_gmres = 1.0d-14
      call lap_bel_solver2(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,eps,numit,rrhs1,eps_gmres,niter,errs,rres,sigma) 

      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)

300      dpars(1) = 1.0d0/4/pi
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


      eps_gmres = 1.0d-14
      call lap_bel_solver2(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,eps,numit,rrhs2,eps_gmres,niter,errs,rres,sigma) 

      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)

      dpars(1) = 1.0d0/4/pi
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


      call surf_grad(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,alpha,sgalpha)


      call surf_grad(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,beta,sgbeta)


      call ncross(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,sgbeta,nsgbeta)

C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        H(1,i) = F(1,i)-sgalpha(1,i)-nsgbeta(1,i)
        H(2,i) = F(2,i)-sgalpha(2,i)-nsgbeta(2,i)
        H(3,i) = F(3,i)-sgalpha(3,i)-nsgbeta(3,i)
      enddo
C$OMP END PARALLEL DO      


      call surf_div(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,H,rrhs1)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=  erra + ((rrhs1(i))**2)*wts(i)
      enddo
      call prin2('harmonic comp error in div =*',erra,1)


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
      call prin2('harmonic comp error in curl =*',erra,1)



      
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


      external xtri_stell_eval,xtri_sphere_eval
      
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
        vmin = 0
        vmax = 2*pi
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
