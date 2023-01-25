      implicit real *8 (a-h,o-z)
      integer norder_start,norder_end
      integer iref,norder0,irep0

      norder_start = 4
      norder_end = 6
      do iref=1,1
        do norder0 = norder_start,norder_end,2
          do irep0=1,3,2
            call solve_sphere(iref,norder0,irep0)
          enddo
        enddo
      enddo
        
          

      stop
      end
c
c
c
c
c

      subroutine solve_sphere(iref,norder0,irep0)
      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer irep0,norder0
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: novers(:),ixyzso(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      real *8, allocatable :: rhs(:),sigma(:),u(:),pot(:)
      real *8, allocatable :: errs(:)

      real *8, allocatable :: ynm(:)
      complex *16, allocatable :: coef(:,:)
      real *8, allocatable :: sigma_ex(:),u_mv_ex(:),u_ex(:),pot_ex(:)
      complex *16, allocatable :: zynm(:,:,:)

      real *8, allocatable :: w(:,:)
      real *8 xyz_out(3),xyz_in(3)
      complex *16 zpars,ima
c
c  variables to report
c     t_mv_q: matvec time (including quadrature generation)
c     t_mv_nq: matvec time (excluding quadrature time)
c     t_solve: solve time
c     niter: number of gmres iterations
c     err_mv: error matvec (error in pot)
c     err_solve: solve error (error in u)
c
c     n_ord, eps_quad, eps_gmres 
c       4, 0.51e-6, 0.5e-6
c       6, 0.51e-8, 0.5e-8
c       8, 0.51e-10, 0.5e-10
c
      real *8 t_mv_q, t_mv_nq, t_solve, err_mv, err_solve

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      igeomtype = 1

      ipars(1) = iref  

      npatches = 12*(4**ipars(1)) 
      norder = norder0 
      npols = (norder+1)*(norder+2)/2
      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0

      if(norder.eq.4.or.norder.eq.5) eps = 0.51d-2
      if(norder.eq.6.or.norder.eq.7) eps = 0.51d-3
      if(norder.ge.8) eps = 0.51d-10
      
      eps_gmres = eps
      irep = irep0

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      ra = 0
      do i=1,npts
        ra = ra + wts(i)
      enddo
      call prin2('surface area of sphere=*',ra,1)
      print *, "npts=",npts


      nmax = 80
      allocate(w(0:nmax,0:nmax))
      allocate(coef(0:nmax,0:nmax))
      allocate(zynm(0:nmax,0:nmax,npts))
      allocate(ynm(npts),sigma_ex(npts),u_ex(npts))
      allocate(pot_ex(npts),u_mv_ex(npts))
      allocate(rhs(npts),u(npts),sigma(npts),pot(npts))


      ru = 0

      call l3getsph_all(nmax,12,srcvals,zynm,npts,w)


      open(unit=34,file='sphere-conv/coefs.dat')
      rewind(34)
      

      do i=0,nmax
        do j=0,nmax
          read(34,*) i1,j1,tmp1,tmp2
          coef(i,j) = tmp1 + ima*tmp2
        enddo
      enddo

      
      close(34)


      do i=1,npts
        pot_ex(i) = 0
        rhs(i) = 0
        sigma_ex(i) = 0
        u_mv_ex(i) = 0
        u_ex(i) = 0
        do j=1,nmax
          do l=0,j
            rtmp = real(zynm(j,l,i)*coef(j,l))
            rhs(i) = rhs(i) + rtmp 

            rpot = -(j+0.0d0)*(j+1.0d0)/(2*j+1.0d0)**2
            if(irep.eq.1.or.irep.eq.2) then
              ru_mv = 1.0d0/(2*j+1.0d0)
            endif
            if(irep.eq.3) then
              ru_mv = 1.0d0/(2*j+1.0d0)**2
            endif
             ru = -1.0d0/(j+0.0d0)/(j+1.0d0)
             if(irep.eq.1.or.irep.eq.2) rsigma = (2*j+1.0d0)*ru
             if(irep.eq.3) rsigma = (2*j+1.0d0)**2*ru

             pot_ex(i) = pot_ex(i) + rtmp*rpot 
             sigma_ex(i) = sigma_ex(i) + rtmp*rsigma
             u_mv_ex(i) = u_mv_ex(i) + rtmp*ru_mv
             u_ex(i) = u_ex(i) + rtmp*ru
          enddo
        enddo
        u(i) = 0
        sigma(i) = 0
      enddo

c
c   first test matvec
c
c
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call lpcomp_lap_bel(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,eps,rhs,irep,pot,u,t_mv_nq)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      t_mv_q = t2-t1
      
      erru = 0
      ra = 0
      errpot = 0
      do i=1,npts
        ra = ra + rhs(i)**2*wts(i)
        erru = erru + (u(i) - u_mv_ex(i))**2*wts(i)
        errpot = errpot + (pot(i)-pot_ex(i))**2*wts(i)
      enddo
      erru = sqrt(erru/ra)
      errpot = sqrt(errpot/ra)
      call prin2('error in pot=*',errpot,1)
      call prin2('error in u=*',erru,1)

      err_mv = errpot

c
c  starting solver test
c

      

      numit = 100
      allocate(errs(numit+1))
      rres = 0
      niter = 0

      call cpu_time(t1)     
C$       t1 = omp_get_wtime()      
      call lap_bel_solver(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,eps,numit,rhs,irep,eps_gmres,niter,errs,rres,sigma,u)
      call cpu_time(t2)     
C$       t2 = omp_get_wtime()      
      t_solve = t2-t1
      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)
      
      
      erru = 0
      ra = 0
      errsig = 0
      do i=1,npts
        ra = ra + rhs(i)**2*wts(i)
        erru = erru + (u(i) - u_ex(i))**2*wts(i)
        errsig = errsig + (sigma(i)-sigma_ex(i))**2*wts(i)
      enddo
      erru = sqrt(erru/ra)
      errsig = sqrt(errsig/ra)
      call prin2('error in sigma=*',errsig,1)
      call prin2('error in u=*',erru,1)

      err_solve = erru
      open(unit=33,file='sphere-conv/res_summary.dat',access='append')
      write(33,'(2x,i1,2x,i1,2x,i7,2x,i4,5(2x,e11.5))') irep,norder,
     1    npts,niter,t_mv_q,t_mv_nq,t_solve,err_mv,err_solve
      close(33)
      
      
      

      return
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
          call h3d_sprime(xyzout,12,srcvals(1,i),0,dpars,1,zk,0,ipars,
     1       val)

          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-1) isout = .false.
      if(abs(ra).le.1.0d-1) isout = .true.

      return
      end

   






      subroutine l3getsph(nmax,mm,nn,ndx,xyzs,ynms,npts,ynm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(npts),ima
      real *8 rat1(10000),rat2(10000)
      complex *16 ynm(0:nmax,0:nmax)
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
        do j=0,nmax
          do l=0,j
            ynm(j,l) = ynm(j,l)*exp(ima*l*phi)
          enddo
        enddo
        ynms(i) = ynm(nn,abs(mm))
      enddo
       
      return
      end





      subroutine l3getsph_all(nmax,ndx,xyzs,ynms,npts,ynm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(0:nmax,0:nmax,npts),ima
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
        do j=0,nmax
          do l=0,j
            ynms(j,l,i) = ynm(j,l)*exp(ima*l*phi)
          enddo
        enddo
      enddo
       
      return
      end




