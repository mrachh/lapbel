      subroutine lpcomp_curl_s_tan_grad(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,eps,sigma,irep,jvec)
c
c  This subroutine evaluates the following layer potential
c  reprensetation
c
c     \nabla \times S [n \times \nabla_{\Gamma} S^{j} \sigma]
c
c  j = 1 or irep, and j=2 for irep=3.
c
c
c  Input arguments:
c 
c    - npatches: integer
c        number of patches
c    - norders: integer(npatches)
c        order of discretization on each patch 
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer(npatches)
c        type of patch
c        iptype = 1, triangular patch discretized using RV nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (9,npts)
c        koornwinder expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. 
c        For each point 
c          * srccoefs(1:3,i) is xyz info
c          * srccoefs(4:6,i) is dxyz/du info
c          * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: real *8 (12,npts)
c        xyz(u,v) and derivative info sampled at the 
c        discretization nodes on the surface
c          * srcvals(1:3,i) - xyz info
c          * srcvals(4:6,i) - dxyz/du info
c          * srcvals(7:9,i) - dxyz/dv info
c          * srcvals(10:12,i) - normals info
c    - eps: real *8
c        precision requested
c    - sigma: real *8 (npts)
c        the density sigma
c    - irep: integer
c        representation for solving the Laplace Beltrami problem
c        * irep = 1, s_lap_s
c        * irep = 3, lap_s2
c
c  Output arguments:
c    - jvec: real *8 (3,npts)
c        the potential pot above
c        
c
      implicit real *8 (a-h,o-z)
      
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      integer, intent(in) :: irep
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: sigma(npts)
      real *8, intent(out) :: jvec(3,npts)
      

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg
      

      integer norder,npols
      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:,:)
      real *8, allocatable :: s_one(:),dens_one(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      real *8 dpars(2),erra,ra
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      real *8 did,dtmp
      integer it,iind,it1,k,l,ndquad
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: wts(:)
      real *8, allocatable :: sigma_use(:)

      complex *16 ztmp


      done = 1
      pi = atan(done)*4

      if(irep.lt.1.or.irep.gt.3) then
        print *, "invalid argument for irep, returning"
        return
      endif



c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
C$OMP END PARALLEL DO   


c
c    initialize patch_id and uv_targ for on surface targets
c
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)

c
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

c
c    find near quadrature correction interactions
c
      print *, "entering find near mem"
      call findnearmem(cms,npatches,rad_near,3,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,3,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)


c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      ztmp = 0
      ikerorder = 0
      print *, eps

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      allocate(wts(npts))

      call get_qwts(npatches,norders,ixyzs,iptype,npts,
     1        srcvals,wts)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)

c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(nquad,4))
      ndquad = 4
      
      do j=1,ndquad
C$OMP PARALLEL DO DEFAULT(SHARED)      
        do i=1,nquad
          wnear(i,j) = 0
        enddo
C$OMP END PARALLEL DO    
      enddo


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call getnearquad_lap_bel_s_lap_s_noc(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
      print *, "done generating near quadrature"

      allocate(sigma_use(npts))

      print *, "Starting to generate right hand side for linear system"


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        sigma_use(i) = sigma(i)
      enddo
C$OMP END PARALLEL DO

      if(irep.eq.3) then
        dpars(1) = 1.0d0
        dpars(2) = 0
        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1      iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2      dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(1,1),
     3      sigma,novers,npts_over,ixyzso,srcover,wover,sigma_use)
       endif


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        jvec(1:3,i) = 0
      enddo
C$OMP END PARALLEL DO

c
c  now sigma_use holds S[\sigma] for irep=3, and \sigma
c  for irep=1, so our job is to evaluate for both
c  reps the following quantity
c
c  \nabla \times S [n \times \nabla_{\Gamma} S[\sigma_use]]
c
c
c  Compute \nabla_{\Gamma} S[\sigma_use]
c
      call lpcomp_curl_s_tan_grad_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2  nquad,wnear,sigma_use,novers,npts_over,ixyzso,
     3  srcover,wover,jvec)


      return
      end
c
c
c
c
c
c
      subroutine lpcomp_curl_s_tan_grad_addsub(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,
     2   iquad,nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover,
     3   jvec)
c
c
c
c  This subroutine evaluates the layer potential for the integral
c  representation:
c
c     \nabla \times S [n \times \nabla_{\Gamma} S \sigma]
c
c  where the near field is precomputed and stored in the row
c  sparse compressed format.
c
c  The subroutine also returns S[\sigma] which is the representation
c  for the solution
c
c  The fmm is used to accelerate the far-field with two calls to the
c  vector fmm, and the near-field interactions are 
c  handled via precomputed quadrature.
c
c  Using add and subtract: no need to call tree and set fmm parameters,
c  can call existing fmm library directly
c  
c
c  Input arguments:
c 
c    - npatches: integer
c        number of patches
c    - norders: integer(npatches)
c        order of discretization on each patch 
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer(npatches)
c        type of patch
c        iptype = 1, triangular patch discretized using RV nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (9,npts)
c        koornwinder expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. 
c        For each point 
c          * srccoefs(1:3,i) is xyz info
c          * srccoefs(4:6,i) is dxyz/du info
c          * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: real *8 (12,npts)
c        xyz(u,v) and derivative info sampled at the 
c        discretization nodes on the surface
c          * srcvals(1:3,i) - xyz info
c          * srcvals(4:6,i) - dxyz/du info
c          * srcvals(7:9,i) - dxyz/dv info
c          * srcvals(10:12,i) - normals info
c    - eps: real *8
c        precision requested
c    - iquadtype: integer
c        quadrature type
c          * iquadtype = 1, use ggq for self + adaptive integration
c            for rest
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(npts+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear_ij array where quadrature for col_ind(i)
c        starts for a single kernel. In this case the different kernels
c        are matrix entries are located at (m-1)*nquad+iquad(i), where
c        m is the kernel number
c    - nquad: integer
c        number of near field entries corresponding to each source target
c        pair. The size of wnear is (nquad,4) since there are 4 kernels
c        per source target pair
c    - wnear: real *8(nquad,4)
c        The desired near field quadrature
c        * The first kernel is S
c        * The second kernel is \nabla_{1} S
c        * The third kernel is \nabla_{2} S
c        * The fourth kernel is \nabla_{3} S 
c    - sigma: real *8(npts)
c        density for the layer potential
c    - novers: integer(npatches)
c        order of discretization for oversampled sources and
c        density
c    - ixyzso: integer(npatches+1)
c        ixyzso(i) denotes the starting location in srcover,
c        corresponding to patch i
c    - nptso: integer
c        total number of oversampled points
c    - srcover: real *8 (12,nptso)
c        oversampled set of source information
c    - whtsover: real *8 (nptso)
c        smooth quadrature weights at oversampled nodes
c    - s_one: real *8 (npts)
c        S[1] on surface
c
c  Output arguments
c    - pot: real *8 (npts)
c        (S \Delta_{\Gamma} S + SWS) \sigma applied using Calderon
c        identities
c------------------------

      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: wnear(nquad,4),sigma(npts)
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      real *8, intent(out) :: jvec(3,npts)

      integer norder,npols,nover,npolso

      integer ndtarg,ntarg

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)

      real *8, allocatable :: curv(:)
      real *8, allocatable :: abc(:,:),pot_tmp(:,:),grad_tmp(:,:,:)
      real *8, allocatable :: hess_tmp(:,:,:),charges_tmp(:,:)
      real *8, allocatable :: dipvec_tmp(:,:,:),abcover(:,:)
      real *8, allocatable :: ffforminv(:,:,:)
      real *8, allocatable :: wtmp1(:,:),wtmp2(:,:)
      real *8, allocatable :: wtmp3(:,:),wtmp4(:,:)
      integer ns,nt
      real *8 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 tmp(10),val(3),grad(3),nvgrad,nvhess
      real *8 hess(6),vgrad(3,3),vhess(3,6)
      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize
      integer i,j,jpatch,jquadstart,jstart

      integer ifaddsub

      integer ntj
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp2(:,:),dtmp2(:,:,:)
      real *8, allocatable :: wts(:)
      real *8 radexp,epsfmm

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime
      real *8 dpars(2)
      real *8 w1,w2,w3,u1,u2,uf,vf

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin,rint,rsurf
      integer nss,ii,l,npover

      integer nd,ntarg0,nmax
      integer ier,iper

      real *8 ttot,done,pi,over4pi

      parameter (ntarg0=1)



      ns = nptso
      done = 1.0d0
      pi = atan(done)*4.0d0
      over4pi = 1.0d0/4.0d0/pi

      nmax = 0
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,
     1   ixyzso,nmax)
      allocate(srctmp2(3,nmax))

      nd = 1
      allocate(abc(3,npts),charges_tmp(nd,ns),dipvec_tmp(nd,3,ns))
      allocate(abcover(3,ns),sigmaover(ns),curv(npts))
      allocate(wtmp1(3,npts),wtmp2(3,npts),ffforminv(2,2,npts))
      allocate(wtmp3(3,npts),wtmp4(3,npts))

      allocate(ctmp2(nd,nmax),dtmp2(nd,3,nmax))
      allocate(pot_tmp(nd,npts),grad_tmp(nd,3,npts),hess_tmp(nd,6,npts))

c
c   Set source info and targinfo for the fmm
c
c
      allocate(sources(3,ns),targvals(3,npts))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        targvals(1,i) = srcvals(1,i)
        targvals(2,i) = srcvals(2,i)
        targvals(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO

c
c  get mean curvature
c
c

      call get_mean_curvature(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,curv)
      call get_inv_first_fundamental_form(npatches,norders,ixyzs, 
     1   iptype,npts,srccoefs,srcvals,ffforminv)
      

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
         wtmp1(1:3,i) = ffforminv(1,1,i)*srcvals(4:6,i) + 
     1      ffforminv(1,2,i)*srcvals(7:9,i)
         wtmp2(1:3,i) = ffforminv(2,1,i)*srcvals(4:6,i) + 
     1      ffforminv(2,2,i)*srcvals(7:9,i)
         call cross_prod3d(srcvals(10,i),wtmp1(1,i),wtmp3(1,i))
         call cross_prod3d(srcvals(10,i),wtmp2(1,i),wtmp4(1,i))
      enddo
C$OMP END PARALLEL DO     
      

c
c
c  oversample the density sigma
c

      call oversample_fun_surf(1,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)
c
c  Set the charge and dipole densities for evaluating
c   
c
c     abc(1:3,:) = \nabla_{\Gamma} S[\sigma]
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        charges_tmp(1,i) = sigmaover(i)*whtsover(i)*over4pi

        dipvec_tmp(1,1,i) = 0
        dipvec_tmp(1,2,i) = 0
        dipvec_tmp(1,3,i) = 0

      enddo
C$OMP END PARALLEL DO


c
c   call the FMM 
c
c

      iper = 0
      ier = 0
      ifpgh = 0
      ifpghtarg = 2
      ifcharge = 1
      ifdipole = 0


      call lfmm3d(nd,eps,ns,sources,ifcharge,charges_tmp,
     1  ifdipole,dipvec_tmp,iper,ifpgh,tmp,tmp,tmp,npts,targvals,
     1  ifpghtarg,pot_tmp,grad_tmp,hess_tmp,ier)

c
c  Now start assembling components
c  abc(1:3,:) will hold n\times \nabla_{\Gamma} S[\sigma]
c 
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,u1,u2)       
      do i=1,npts
         call dot_prod3d(grad_tmp(1,1:3,i),srcvals(4,i),u1)
         call dot_prod3d(grad_tmp(1,1:3,i),srcvals(7,i),u2)
         abc(1:3,i) = u1*wtmp3(1:3,i) + u2*wtmp4(1:3,i)
      enddo
C$OMP END PARALLEL DO     
c
c  Fix quadrature corrections
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,npols,l,w1,w2,w3,uf,vf)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            w1 = wnear(jquadstart+l-1,2)
            w2 = wnear(jquadstart+l-1,3)
            w3 = wnear(jquadstart+l-1,4)
            uf = w1*srcvals(4,i) + w2*srcvals(5,i) + w3*srcvals(6,i)
            vf = w1*srcvals(7,i) + w2*srcvals(8,i) + w3*srcvals(9,i)
            abc(1:3,i) = abc(1:3,i) + 
     1         (uf*wtmp3(1:3,i)+vf*wtmp4(1:3,i))*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c      print *, "Near quad addition done"

      thresh = 0.0d0
      call get_fmm_thresh(3,ns,sources,3,npts,targvals,thresh)

c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,grad,u1,u2,npover)
      do i=1,npts

        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
            ii = ii+1
            srctmp2(1,ii) = srcover(1,jstart+l)
            srctmp2(2,ii) = srcover(2,jstart+l)
            srctmp2(3,ii) = srcover(3,jstart+l)
            

            ctmp2(1,ii) = charges_tmp(1,jstart+l)

          enddo
        enddo
        nss = ii

        val(1) = 0
        grad(1:3) = 0


        call l3ddirectcg(nd,srctmp2,ctmp2,
     1        nss,targvals(1,i),ntarg0,val,grad,thresh)
        call dot_prod3d(grad,srcvals(4,i),u1)
        call dot_prod3d(grad,srcvals(7,i),u2)

       
        abc(1:3,i) = abc(1:3,i) - (u1*wtmp3(1:3,i) + u2*wtmp4(1:3,i)) 
      enddo
C$OMP END PARALLEL DO     

c      print *, "Subtraction done"
c
c    abc(1:3,:) = n\times\nabla_{\Gamma} S[\sigma]
c

      nd = 3
      deallocate(ctmp2,dtmp2)
      deallocate(pot_tmp,grad_tmp,hess_tmp)
      deallocate(charges_tmp,dipvec_tmp)

      allocate(ctmp2(nd,nmax),dtmp2(nd,3,nmax))
      allocate(charges_tmp(nd,ns),dipvec_tmp(nd,3,ns))
      allocate(pot_tmp(nd,npts),grad_tmp(nd,3,npts),hess_tmp(nd,6,npts))


      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype, 
     1    npts,abc,novers,ixyzso,ns,abcover)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        charges_tmp(1:3,i) = abcover(1:3,i)*whtsover(i)*over4pi
        dipvec_tmp(1:3,1:3,i) = 0
      enddo
C$OMP END PARALLEL DO

c
c   call the FMM 
c
c

      iper = 0
      ier = 0
      ifpgh = 0
      ifpghtarg = 2
      ifcharge = 1
      ifdipole = 0
      call lfmm3d(nd,eps,ns,sources,ifcharge,charges_tmp,
     1  ifdipole,dipvec_tmp,iper,ifpgh,tmp,tmp,tmp,npts,targvals,
     1  ifpghtarg,pot_tmp,grad_tmp,hess_tmp,ier)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)       
      do i=1,npts
         jvec(1,i) = grad_tmp(3,2,i) - grad_tmp(2,3,i)
         jvec(2,i) = grad_tmp(1,3,i) - grad_tmp(3,1,i)
         jvec(3,i) = grad_tmp(2,1,i) - grad_tmp(1,2,i)
      enddo
C$OMP END PARALLEL DO      

c
c  Fix quadrature corrections
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,l,npols,w1,w2,w3)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(jquadstart+l-1,2)
            w2 = wnear(jquadstart+l-1,3)
            w3 = wnear(jquadstart+l-1,4)
            jvec(1,i) = jvec(1,i) + w2*abc(3,jstart+l-1) - 
     1           w3*abc(2,jstart+l-1)
            jvec(2,i) = jvec(2,i) + w3*abc(1,jstart+l-1) - 
     1           w1*abc(3,jstart+l-1)
            jvec(3,i) = jvec(3,i) + w1*abc(2,jstart+l-1) - 
     1           w2*abc(1,jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c      print *, "Near quad addition done"



c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,vgrad,npover)
      do i=1,npts

        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
            ii = ii+1
            srctmp2(1,ii) = srcover(1,jstart+l)
            srctmp2(2,ii) = srcover(2,jstart+l)
            srctmp2(3,ii) = srcover(3,jstart+l)
            
            ctmp2(1,ii) = charges_tmp(1,jstart+l)
            ctmp2(2,ii) = charges_tmp(2,jstart+l)
            ctmp2(3,ii) = charges_tmp(3,jstart+l)

          enddo
        enddo
        nss = ii

        val(1) = 0
        val(2) = 0
        val(3) = 0
        vgrad(1:3,1:3) = 0

        call l3ddirectcg(nd,srctmp2,ctmp2,
     1        nss,targvals(1,i),ntarg0,val,vgrad,thresh)
         jvec(1,i) = jvec(1,i) - (vgrad(3,2)-vgrad(2,3))
         jvec(2,i) = jvec(2,i) - (vgrad(1,3)-vgrad(3,1))
         jvec(3,i) = jvec(3,i) - (vgrad(2,1)-vgrad(1,2))
      enddo
C$OMP END PARALLEL DO

      return
      end

c
c
c
c
c
c
