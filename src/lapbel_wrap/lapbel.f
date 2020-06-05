c
c     Laplace Beltrami solver, based on the representation
c     S_{0} (\Delta_{\Gamma}+ \int_{\Gamma} )S_{0} \sigma = f
c
c     it stores five sets of quadratures in wnear, 
c       \nabla_{t} S_{0} \cdot x_{u}, \nabla_{t} S_{0} \cdot x_{v},
c       \nabla_{s} S_{0} \cdot x_{u}, \nabla_{s} S_{0} \cdot x_{v},
c       and finally S_{0}
c    
c
c    This way there is only loss of one order of accuracy
c    (hopefully)
c
c     This file contains the following user callable
c     routines: 
c 
c       getnearquad_lap_bel - generates the near
c        field quadrature for applying S_{0} \Delta_{\Gamma} S_{0},
c        see note above
c
c       lpcomp_lap_bel 
c          simpler version of Laplace layer potential evaluator
c          only geometry, targets, representation parameters (alpha,beta)
c          and density sampled at discretization required on input,
c          output is the application of S_{0} \Delta S_{0} 
c          (note that the identity term is not included for targets on
c           surface)
c
c       lap_bel_solver - Solves the laplace beltrami integral
c         equation 
c
c
c    Advanced user interfaces: 
c*****************************
c       Note for developers: One of the two advanced user interfaces
c         is necessary for the easy user interface and
c         efficient iterative solvers. It seems that the add and subtract
c         version tends to have better CPU-time performance, but we expect
c         the setsub version to be more numerically stable
c**************************************
c       lpcomp_lap_bel_addsub
c         Apply S_{0} \Delta S_{0} using add and subtract
c



      subroutine getnearquad_lap_bel(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,
     2   ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
c
c       this subroutine generates the near field quadrature
c       for the representation u = S_{0} \Delta_{\Gamma} S_{0} \, ,
c       by storing (\nabla_{t} S_{0} [\sigma]) \cdot x_{u}),
c       (\nabla_{t} S_{0} [\sigma]) \cdot x_{v}), 
c       (\nabla_{s} S_{0} \cdot x_{u} [\sigma]), 
c       (\nabla_{s} S_{0} \cdot x_{v} [\sigma]) and that for S_{0}
c       where the near field is specified by the user 
c       in row sparse compressed format.
c
c
c       The quadrature is computed by the following strategy
c        targets within a sphere of radius rfac0*rs
c        of a chunk centroid is handled using adaptive integration
c        where rs is the radius of the bounding sphere
c        for the patch
c  
c       All other targets in the near field are handled via
c        oversampled quadrature
c
c       The recommended parameter for rfac0 is 1.25d0
c
c       Note that all the parameters are real and the quadrature
c       correction returned is real. Currently implemented in
c       a hacky manner by calling the complex routine and then
c       converting to real routines
c        
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders - integer(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer(npatches+1)
c            starting location of data on patch i
c  
c         iptype - integer(npatches)
c           type of patch
c           iptype = 1 -> triangular patch discretized with RV nodes
c
c         npts - integer
c            total number of discretization points on the boundary
c
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du,
c            and dxyz/dv on each patch. 
c            For each point srccoefs(1:3,i) is xyz info
c                           srccoefs(4:6,i) is dxyz/du info
c                           srccoefs(7:9,i) is dxyz/dv info
c
c          srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c             srcvals(10:12,i) - normals info
c 
c         ipatch_id - integer(npts)
c            id of patch of target i, id = -1, if target is off-surface
c
c         uvs_targ - real *8 (2,npts)
c            local uv coordinates on patch if on surface, otherwise
c            set to 0 by default
c            
c         eps - real *8
c             precision requested
c
c         iquadtype - integer
c           quadrature type
c           iquadtype = 1, use ggq for self + adaptive integration
c                 for rest
c 
c
c           nnz - integer
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer(npts+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           rfac0 - integer
c               radius parameter for near field
c
c           nquad - integer
c               number of entries in wnear for one kernel
c
c        output
c            wnear - real *8(nquad*5)
c               the desired near field quadrature
c               
c

      implicit none 
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: iquadtype
      integer, intent(in) :: ipatch_id(npts)
      real *8, intent(in) :: uvs_targ(2,npts)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(5*nquad)


      integer ipars
      integer ndd,ndz,ndi
      complex *16 zpars
      real *8 dpars

      integer ndtarg

      real *8 alpha,beta,done,pi
      integer i,j
      integer ipv

      procedure (), pointer :: fker
      external l3d_slp, l3d_sgradu, l3d_sgradv,l3d_dgradu,l3d_dgradv

      done = 1
      pi = atan(done)*4

c
c
c        initialize the appropriate kernel function
c


      ndd = 0
      ndi = 0
      ndz = 0
      ndtarg = 12
      if(iquadtype.eq.1) then
        ipv = 1
        fker=>l3d_sgradu
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1))
        print *, "done with kernel 1"

        fker=>l3d_sgradv
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(nquad+1))
        print *, "done with kernel 2"

        fker=>l3d_dgradu
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(2*nquad+1))
        print *, "done with kernel 3"

        fker=>l3d_dgradv
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(3*nquad+1))
        print *, "done with kernel 4"

        fker=> l3d_slp
        ipv = 0
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(4*nquad+1))
        print *, "done with kernel 5"
      endif

      dpars = 1.0d0/4/pi

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,5*nquad
        wnear(i) = wnear(i)*dpars
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
      subroutine lpcomp_lap_bel(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,eps,sigma,pot)
c
cf2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
cf2py intent(in) eps
cf2py intent(in) sigma
cf2py intent(out) pot
c
c
c------------------------------
c  This subroutine evaluates the layer potential for the representation 
c
c
c  .. math ::
c  
c      u = S_{0} \Delta_{\Gamma} S_{0}
c
c  
c  Note: The subroutine only computes 
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential (I/4) is not included in the layer potential.
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
c        and srcvals array where information for patch i begins
c    - iptype: integer(npatches)
c        type of patch
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        koornwinder expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$. 
c    - srcvals: double precision (12,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - eps: double precision
c        precision requested
c     - sigma: double precision(npts)
c         density for layer potential
c
c  Output arguments
c    - pot: double precision(ntarg)
c        layer potential evaluated at the target points
c
c-----------------------------------
c
      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: sigma(npts)

      real *8, intent(out) :: pot(npts)

      real *8, allocatable :: targs(:,:),uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)
      integer nptso,nnz,nquad


      integer nover,npolso
      integer norder,npols
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      complex *16 zpars
      integer ipars
      real *8 timeinfo(10),t1,t2,omp_get_wtime
      integer ntarg,ndtarg


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c        this might need fixing
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
      allocate(wnear(5*nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,5*nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    

      call prinf('finished generating near field info*',i,0)
      call prinf('finished generating far field orders*',i,0)
      call prinf('npts_over=*',npts_over,1)
      call prin2('eps=*',eps,1)

      iquadtype = 1

      call getnearquad_lap_bel(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
      call prinf('finished generating near quadrature correction*',i,0)


c
c
c   compute layer potential
c
      call lpcomp_lap_bel_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,
     2  eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3  sigma,novers,npts_over,ixyzso,srcover,wover,pot)



      return
      end
c
c
c
c
c
      subroutine lpcomp_lap_bel_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,
     2   eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,
     3   nptso,ixyzso,srcover,whtsover,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the representation u = S_{0} \Delta_{\Gamma} S_{0} 
c      where the near field is precomputed and stored
c      in the row sparse compressed format.
c
c     The fmm is used to accelerate the far-field and 
c     near-field interactions are handled via precomputed quadrature
c
c     Using add and subtract - no need to call tree and set fmm parameters
c      can directly call existing fmm library
c
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders- integer(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer(npatches+1)
c            ixyzs(i) denotes the starting location in srccoefs,
c               and srcvals array corresponding to patch i
c   
c         iptype - integer(npatches)
c            type of patch
c             iptype = 1, triangular patch discretized using RV nodes
c
c         npts - integer
c            total number of discretization points on the boundary
c 
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du,
c            and dxyz/dv on each patch. 
c            For each point srccoefs(1:3,i) is xyz info
c                           srccoefs(4:6,i) is dxyz/du info
c                           srccoefs(7:9,i) is dxyz/dv info
c
c         srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c             srcvals(10:12,i) - normals info
c 
c          eps - real *8
c             precision requested
c
c           nnz - integer *8
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer(ntarg+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           nquad - integer
c               number of entries in wnear
c
c           wnear - real *8(5*nquad)
c               the near field quadrature correction
c
c           sigma - real *8(npts)
c               density for layer potential
c
c           novers - integer(npatches)
c              order of discretization for oversampled sources and
c               density
c
c         ixyzso - integer(npatches+1)
c            ixyzso(i) denotes the starting location in srcover,
c               corresponding to patch i
c   
c           nptso - integer
c              total number of oversampled points
c
c           srcover - real *8 (12,nptso)
c              oversampled set of source information
c
c           whtsover - real *8 (nptso)
c             smooth quadrature weights at oversampled nodes
c
c
c         output
c           pot - real *8(npts)
c              layer potential evaluated at the target points
c
c           
c               
c
      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: wnear(5*nquad),sigma(npts)
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      real *8, intent(out) :: pot(npts)

      integer norder,npols,nover,npolso

      integer ndtarg,ntarg

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)

      real *8, allocatable :: pot1(:),grad1(:,:),gradsurf(:,:)
      real *8, allocatable :: gradsurf2(:,:)
      real *8, allocatable :: gradsurf2over(:,:)
      real *8, allocatable :: ffforminv(:,:,:)
      
      integer ns,nt
      real *8 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 tmp(10),val,vgrad(3)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize
      integer i,j,jpatch,jquadstart,jstart

      integer ifaddsub

      integer ntj
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp2(:),dtmp2(:,:)
      real *8, allocatable :: wts(:)
      real *8 radexp,epsfmm

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime
      real *8 dpars(2)

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin,rint
      integer nss,ii,l,npover

      integer nd,ntarg0

      real *8 ttot,done,pi,over4pi

      parameter (nd=1,ntarg0=1)


c
c
c       first compute \nabla_{\Gamma} S_{0}
c

      ns = nptso
      done = 1
      pi = atan(done)*4
      over4pi = 1.0d0/4/pi

           
      ifpgh = 0
      ifpghtarg = 2
      ntarg = npts
      ndtarg = 3
      allocate(sources(3,ns),targvals(3,npts))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))

      allocate(pot1(npts),grad1(3,npts),gradsurf(2,npts))

c 
c       oversample density
c

      call oversample_fun_surf(1,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)

      ra = 0


c
c       set relevatn parameters for the fmm
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        charges(i) = sigmaover(i)*whtsover(i)*over4pi
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        targvals(1,i) = srcvals(1,i)
        targvals(2,i) = srcvals(2,i)
        targvals(3,i) = srcvals(3,i)

        grad1(1,i) = 0
        grad1(2,i) = 0
        grad1(3,i) = 0
      enddo
C$OMP END PARALLEL DO      

      ifcharge = 1
      ifdipole = 0

c
c
c       call the fmm
c

      call lfmm3d(nd,eps,ns,sources,ifcharge,charges,
     1  ifdipole,dipvec,ifpgh,tmp,tmp,tmp,npts,targvals,ifpghtarg,
     1  pot1,grad1,tmp)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        gradsurf(1,i) = grad1(1,i)*srcvals(4,i) + 
     1                  grad1(2,i)*srcvals(5,i) + 
     2                  grad1(3,i)*srcvals(6,i)

        gradsurf(2,i) = grad1(1,i)*srcvals(7,i) + 
     1                  grad1(2,i)*srcvals(8,i) + 
     2                  grad1(3,i)*srcvals(9,i)
      enddo
C$OMP END PARALLEL DO


c
c        compute threshold for ignoring local computation
c
      
      xmin = sources(1,1)
      xmax = sources(1,1)
      ymin = sources(2,1)
      ymax = sources(2,1)
      zmin = sources(3,1)
      zmax = sources(3,1)

      do i=1,ns
        if(sources(1,i).lt.xmin) xmin = sources(1,i)
        if(sources(1,i).gt.xmax) xmax = sources(1,i)
        if(sources(2,i).lt.ymin) ymin = sources(2,i)
        if(sources(2,i).gt.ymax) ymax = sources(2,i)
        if(sources(3,i).lt.zmin) zmin = sources(3,i)
        if(sources(3,i).gt.zmax) zmax = sources(3,i)
      enddo

      do i=1,npts
        if(targvals(1,i).lt.xmin) xmin = targvals(1,i)
        if(targvals(1,i).gt.xmax) xmax = targvals(1,i)
        if(targvals(2,i).lt.ymin) ymin = targvals(2,i)
        if(targvals(2,i).gt.ymax) ymax = targvals(2,i)
        if(targvals(3,i).lt.zmin) zmin = targvals(3,i)
        if(targvals(3,i).gt.zmax) zmax = targvals(3,i)
      enddo
      
      boxsize = xmax-xmin
      sizey = ymax - ymin
      sizez = zmax - zmin

      if(sizey.gt.boxsize) boxsize = sizey
      if(sizez.gt.boxsize) boxsize = sizez

      thresh = 2.0d0**(-51)*boxsize
c
c
c       add in precomputed quadrature
c



      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            gradsurf(1,i) = gradsurf(1,i) + 
     1          wnear(jquadstart+l-1)*sigma(jstart+l-1)
            gradsurf(2,i) = gradsurf(2,i) + 
     1          wnear(nquad+jquadstart+l-1)*sigma(jstart+l-1)
            pot1(i) = pot1(i) + 
     1          wnear(4*nquad+jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO


c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,nss,l,jstart,ii,val,npover,vgrad)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(3,nss),ctmp2(nss))

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

            ctmp2(ii) = charges(jstart+l)
          enddo
        enddo

        val = 0
        vgrad(1) = 0
        vgrad(2) = 0
        vgrad(3) = 0

        call l3ddirectcg(nd,srctmp2,ctmp2,
     1        nss,targvals(1,i),ntarg0,val,vgrad,thresh)
        gradsurf(1,i) = gradsurf(1,i) - vgrad(1)*srcvals(4,i)
     1        - vgrad(2)*srcvals(5,i) - vgrad(3)*srcvals(6,i)
        gradsurf(2,i) = gradsurf(2,i) - vgrad(1)*srcvals(7,i)
     1        - vgrad(2)*srcvals(8,i) - vgrad(3)*srcvals(9,i)
        pot1(i) = pot1(i) - val
        deallocate(srctmp2,ctmp2)
      enddo

c
c      convert to proper surface gradient by multiplying by 
c       inverse of first fundamental form
c            
      allocate(ffforminv(2,2,npts),gradsurf2(2,npts))
      call get_inv_first_fundamental_form(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ffforminv)

      do i=1,npts
        gradsurf2(1,i) = ffforminv(1,1,i)*gradsurf(1,i) + 
     1     ffforminv(1,2,i)*gradsurf(2,i)
        gradsurf2(2,i) = ffforminv(2,1,i)*gradsurf(1,i) + 
     1     ffforminv(2,2,i)*gradsurf(2,i)
      enddo

      rint = 0
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,
     1        srcvals,wts)
      do i=1,npts
        rint = rint + pot1(i)*wts(i)
      enddo

      do i=1,npts
        pot1(i) = rint
      enddo

      dpars(1) = over4pi
      dpars(2) = 0

      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targvals,eps,
     2  dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(4*nquad+1),
     3  pot1,novers,nptso,ixyzso,srcover,whtsover,pot)

      allocate(gradsurf2over(2,ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype, 
     1    npts,gradsurf2,novers,ixyzso,ns,gradsurf2over)

      do i=1,npts
        pot1(i) = 0
      enddo

      do i=1,ns
        dipvec(1,i) = (gradsurf2over(1,i)*srcover(4,i)+
     1      gradsurf2over(2,i)*srcover(7,i))*whtsover(i)*over4pi
        dipvec(2,i) = (gradsurf2over(1,i)*srcover(5,i)+
     1      gradsurf2over(2,i)*srcover(8,i))*whtsover(i)*over4pi
        dipvec(3,i) = (gradsurf2over(1,i)*srcover(6,i)+
     1      gradsurf2over(2,i)*srcover(9,i))*whtsover(i)*over4pi
        charges(i) = 0
      enddo

      ifcharge = 0
      ifdipole = 1

      ifpghtarg = 1
c
c
c       call the fmm
c

      call lfmm3d(nd,eps,ns,sources,ifcharge,charges,
     1  ifdipole,dipvec,ifpgh,tmp,tmp,tmp,npts,targvals,ifpghtarg,
     1  pot1,grad1,tmp)


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot1(i) = pot1(i) + 
     1          wnear(2*nquad+jquadstart+l-1)*gradsurf2(1,jstart+l-1)
            pot1(i) = pot1(i) + 
     1          wnear(3*nquad+jquadstart+l-1)*gradsurf2(2,jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO

c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(dtmp2,nss,l,jstart,ii,val,npover)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(3,nss),dtmp2(3,nss))

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

            dtmp2(1,ii) = dipvec(1,jstart+l)
            dtmp2(2,ii) = dipvec(2,jstart+l)
            dtmp2(3,ii) = dipvec(3,jstart+l)
          enddo
        enddo

        val = 0
        call l3ddirectdp(nd,srctmp2,dtmp2,
     1        nss,targvals(1,i),ntarg0,val,thresh)

        pot1(i) = pot1(i) - val

        deallocate(srctmp2,dtmp2)
      enddo
      
      do i=1,npts
        pot(i) = pot(i) - pot1(i)
        pot(i) = -pot1(i)
      enddo


      return
      end

c
c
c
c
c
c
c
c        
      subroutine lap_bel_solver(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,numit,
     2    rhs,eps_gmres,niter,errs,rres,soln)
c
c
c        this subroutine solves the Laplace beltrami equation
c        using the representation S_{0} \Delta_{\Gamma} S_{0} \sigma 
c         = S f
c
c
c     The linear system is solved iteratively using GMRES
c
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders- integer(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer(npatches+1)
c            ixyzs(i) denotes the starting location in srccoefs,
c               and srcvals array corresponding to patch i
c   
c         iptype - integer(npatches)
c            type of patch
c             iptype = 1, triangular patch discretized using RV nodes
c
c         npts - integer
c            total number of discretization points on the boundary
c 
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du,
c            and dxyz/dv on each patch. 
c            For each point srccoefs(1:3,i) is xyz info
c                           srccoefs(4:6,i) is dxyz/du info
c                           srccoefs(7:9,i) is dxyz/dv info
c
c         srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c 
c          eps - real *8
c             precision requested for computing quadrature and fmm
c             tolerance
c
c           rhs - real *8(npts)
c              right hand side
c
c           eps_gmres - real *8
c                gmres tolerance requested
c
c           numit - integer
c              max number of gmres iterations
c
c         output
c           niter - integer
c              number of gmres iterations required for relative residual
c          
c           errs(1:iter) - relative residual as a function of iteration
c              number
c 
c           rres - real *8
c              relative residual for computed solution
c              
c           soln - real *8(npts)
c              density which solves the dirichlet problem
c
c
      implicit none
      integer npatches,norder,npols,npts
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      real *8 rhs(npts)
      real *8 soln(npts)

      real *8, allocatable :: rhs2(:)
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:)

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
      real *8 rb,wnrm2
      integer numit,it,iind,it1,k,l
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: vmat(:,:),hmat(:,:)
      real *8, allocatable :: cs(:),sn(:)
      real *8, allocatable :: svec(:),yvec(:),wtmp(:)

      complex *16 ztmp


      allocate(vmat(npts,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(npts),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4


c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))

      allocate(rhs2(npts))

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
      call findnearmem(cms,npatches,rad_near,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,targs,npts,row_ptr, 
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
      call prinf('novers=*',novers,100)

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(5*nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,5*nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call getnearquad_lap_bel(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)

c
c
c     compute the right hand side S_{0}[f], note that
c     quadrature is precomputed in wnear[2*nquad+1:3*nquad]   
c
      dpars(1) = 1.0d0/4/pi
      dpars(2) = 0
      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2  dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(4*nquad+1),
     3  rhs,novers,npts_over,ixyzso,srcover,wover,rhs2)
      erra = 0
      ra = 0
      do i=1,npts
        erra = erra + (rhs(i)/5.0d0-rhs2(i))**2
        ra = ra + rhs2(i)**2
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in rhs=*',erra,1)


      
      print *, "done generating near quadrature, now starting gmres"


c
c
c     start gmres code here
c
c     NOTE: matrix equation should be of the form (z*I + K)x = y
c       the identity scaling (z) is defined via did below,
c       and K represents the action of the principal value 
c       part of the matvec
c
      did = -0.25d0
      did = 0


      niter=0

c
c      compute norm of right hand side and initialize v
c 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo


c
      do i=1,npts
        rb = rb + abs(rhs2(i))**2
      enddo
      rb = sqrt(rb)

      do i=1,npts
        vmat(i,1) = rhs2(i)/rb
      enddo

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


        call lpcomp_lap_bel_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,
     2    eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3    vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp)

        do k=1,it
          hmat(k,it) = 0
          do j=1,npts
            hmat(k,it) = hmat(k,it) + wtmp(j)*vmat(j,k)
          enddo

          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
        enddo
          
        hmat(it,it) = hmat(it,it)+did
        wnrm2 = 0
        do j=1,npts
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
        wnrm2 = sqrt(wnrm2)

        do j=1,npts
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        dtmp = wnrm2

        call rotmat_gmres(hmat(it,it),dtmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

c
c            solve the linear system corresponding to
c            upper triangular part of hmat to obtain yvec
c
c            y = triu(H(1:it,1:it))\s(1:it);
c
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo



c
c          estimate x
c
          do j=1,npts
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo


          rres = 0
          do i=1,npts
            wtmp(i) = 0
          enddo
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


          call lpcomp_lap_bel_addsub(npatches,norders,ixyzs,
     1      iptype,npts,srccoefs,srcvals,
     2      eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3      soln,novers,npts_over,ixyzso,srcover,wover,wtmp)
            
          do i=1,npts
            rres = rres + abs(did*soln(i) + wtmp(i)-rhs(i))**2
          enddo
          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo
c
      return
      end
c
