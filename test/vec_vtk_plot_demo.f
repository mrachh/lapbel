      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: sigma(:,:)
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,fname

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima

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

      fname = '../../fmm3dbie/geometries/sphere_192_o07.go3'
      
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,
     1   iptype,npts,srcvals,srccoefs,wts)


      allocate(sigma(3,npts))
      ifinout = 1

      do i=1,npts
        sigma(1,i) = srcvals(10,i)
        sigma(2,i) = srcvals(11,i)
        sigma(3,i) = srcvals(12,i)
      enddo

      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,sigma,'vecplot.vtk','a') 


      stop
      end

