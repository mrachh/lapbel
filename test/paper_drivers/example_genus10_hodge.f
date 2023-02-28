      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: novers(:),ixyzso(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      real *8, allocatable :: rhs(:,:),cfree(:,:),dfree(:,:)
      real *8, allocatable :: vharm(:,:),vharm_cross(:,:)
      real *8, allocatable :: vharm_div(:),vharm_cross_div(:)
      
      real *8, allocatable :: errs(:,:)

      integer niter(2)
      real *8 rres(2)

      real *8, allocatable :: ynm(:),unm(:,:),xnm(:,:)
      complex *16, allocatable :: zynm(:)
      real *8, allocatable :: pols1(:),pols2(:),pols3(:)

      real *8, allocatable :: w(:,:)
      real *8 xyz_out(3),xyz_in(3),wtmp(3),dl(3)
      complex *16 zpars

      character *300 dirname,fname_rhs,fname_vharm
      character *300 fname_geom,fname_geom_dir

      integer pdeg(3)


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      iref = 2
      norder = 8
      npols = (norder+1)*(norder+2)/2
      
      fname_geom_dir = 'genus10-geometries/'
      write(fname_geom,'(a,a,i2.2,a)') trim(fname_geom_dir),
     1  'Genus_10_o08_r',iref,'.go3'
      print *, trim(fname_geom)
      call open_gov3_geometry_mem(trim(fname_geom),npatches,npts)
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(trim(fname_geom),npatches,norders,ixyzs,
     1  iptype,npts,srcvals,srccoefs,wts)


      pdeg(1) = 1
      pdeg(2) = 1
      pdeg(3) = 4

      npoly1 = (pdeg(1)+1)*(pdeg(1)+2)*(pdeg(1)+3)/6
      npoly2 = (pdeg(2)+1)*(pdeg(2)+2)*(pdeg(2)+3)/6
      npoly3 = (pdeg(3)+1)*(pdeg(3)+2)*(pdeg(3)+3)/6
      allocate(rhs(3,npts),cfree(3,npts),dfree(3,npts))
      allocate(vharm(3,npts),vharm_cross(3,npts))
      allocate(vharm_div(npts),vharm_cross_div(npts))
      allocate(pols1(npoly1),pols2(npoly2),pols3(npoly3))
      
      rmax = 0
      ra = 0
      do i=1,npts
        dx = 0.7d0*(0.2d0*srcvals(1,i)-1.0d0) 
        dy = 0.7d0*(1.0d0/7.0d0)*srcvals(2,i)
        dz = 0.7d0*(2*srcvals(3,i)-1.0d0)
        
        if(abs(dx).ge.rmax) rmax = abs(dx)
        if(abs(dy).ge.rmax) rmax = abs(dy)
        if(abs(dz).ge.rmax) rmax = abs(dz)
        call legetens_pols_3d(dx,pdeg(1),'T',pols1)
        call legetens_pols_3d(dy,pdeg(2),'T',pols2)
        call legetens_pols_3d(dz,pdeg(3),'T',pols3)
        rhs(1,i) = pols1(npoly1/2)
        rhs(2,i) = pols2(npoly2/2)
        rhs(3,i) = pols3(npoly3/2)

        ra = ra + rhs(1,i)**2*wts(i)
        ra = ra + rhs(2,i)**2*wts(i)
        ra = ra + rhs(3,i)**2*wts(i)
        
        cfree(1:3,i) = 0
        dfree(1:3,i) = 0
        vharm(1:3,i) = 0
      enddo

      ra = sqrt(ra)
      print *, "norm of rhs=",ra



      eps = 0.5d-7
      eps_gmres = 0.5d-9
      irep = 3


      numit = 400
      allocate(errs(numit+1,2))
      rres(1:2) = 0
      niter(1:2) = 0

      call cpu_time(t1)
C$       t1 = omp_get_wtime() 
      call get_hodge_decomposition(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,numit,rhs,irep,eps_gmres,niter,errs,rres,
     2  cfree,dfree,vharm)
      
      call cpu_time(t2)
C$       t2 = omp_get_wtime() 

      do i=1,npts
        call cross_prod3d(srcvals(10,i),vharm(1,i),vharm_cross(1,i))
      enddo

      call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,vharm,vharm_div)
      
      call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,vharm_cross,vharm_cross_div)
      
      err1 = 0
      err2 = 0
      do i=1,npts
        err1 = err1 + vharm_div(i)**2*wts(i)
        err2 = err2 + vharm_cross_div(i)**2*wts(i) 
      enddo

      err1 = sqrt(err1)
      err2 = sqrt(err2)
      call prin2('error in div hvec',err1,1)
      call prin2('error in n times div hvec',err2,1)


      open(unit=33,file='genus10-hodge/res_summary.dat',access='append')
      write(33,'(2x,i1,2x,i1,2x,i7,2(2x,i3),3(2x,e11.5))') irep,norder,
     1   npts,niter(1),niter(2),t2-t1,err1,err2
      close(33)
      
      dirname = 'genus10-hodge/'
      write(fname_rhs,'(a,a,i1,a,i1,a)') trim(dirname),
     1   'rhs_iref',iref,'_norder',norder,'.dat'
      open(unit=78,file=trim(fname_rhs),form='unformatted')
      write(78) rhs
      close(78)
      
      write(fname_vharm,'(a,a,i1,a,i1,a,i1,a)') trim(dirname),
     1   'hvec_iref',iref,'_norder',norder,'_irep',irep,
     2   '.dat'
      open(unit=78,file=trim(fname_vharm),form='unformatted')
      write(78) vharm
      close(78)
      

      

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



