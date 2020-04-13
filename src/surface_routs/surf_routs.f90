!
!
!   This file contains the following user callable routines
!   
!   get_first_fundamental_form - get the first fundamental form at all
!      points
!
!   get_inv_first_fundamental_form - get inverse of first fundamental 
!      form
!
!   get_surf_grad - compute surface gradient of a scalar function
!   get_surf_grad_fast - compute surface gradient of a scalar function 
!      (With precomputed inverse of first fundamental form)
!
!   get_surf_uv_grad - compute the uv gradient of a function defined
!      on a surface
!
!   get_surf_div - compute the surface divergence of a vector field
!       defined on a surface
!

subroutine get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)
  
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts
  real *8 srccoefs(9,npts),srcvals(12,npts),ffform(2,2,npts)
  integer i
  real *8 xuu,xuv,xuv

  do i=1,npts
    xuu = srcvals(4,i)**2 + srcvals(5,i)**2 + srcvals(6,i)**2
    xuv = srcvals(7,i)**2 + srcvals(8,i)**2 + srcvals(9,i)**2
    xuv = srcvals(4,i)*srcvals(7,i) + srcvals(5,i)*srcvals(8,i)+ &
       srcvals(6,i)*srcvals(9,i)
    ffform(1,1,i) = xuu
    ffform(2,1,i) = xuv
    ffform(1,2,i) = xuv
    ffform(2,2,i) = xvv
  enddo

  return
end subroutine get_first_fundamental_form
!
!
!
!
!

subroutine get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)
  
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts
  real *8 srccoefs(9,npts),srcvals(12,npts),ffforminv(2,2,npts)
  integer i
  real *8 xuu,xuv,xuv,d

  do i=1,npts
    xuu = srcvals(4,i)**2 + srcvals(5,i)**2 + srcvals(6,i)**2
    xuv = srcvals(7,i)**2 + srcvals(8,i)**2 + srcvals(9,i)**2
    xuv = srcvals(4,i)*srcvals(7,i) + srcvals(5,i)*srcvals(8,i)+ &
       srcvals(6,i)*srcvals(9,i)
    d = xuu*xvv - xuv*xuv
    ffforminv(1,1,i) = xvv/d
    ffforminv(2,1,i) = -xuv/d
    ffforminv(1,2,i) = -xuv/d
    ffforminv(2,2,i) = xuu/d
  enddo

  return
end subroutine get_inv_first_fundamental_form
!
!
!
!

subroutine get_surf_grad(nd,npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,f,df)
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts
  real *8 srccoefs(9,npts),srcvals(12,npts),f(nd,npts),df(nd,2,npts)
  real *8, allocatable :: ffforminv(:,:,:)

  allocate(ffforminv(2,2,npts))

  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  call get_surf_grad_fast(nd,npatches,norders,ixyzs,iptype,npts, &
    ffforminv,f,df)

  return
end subroutine get_surf_grad
!
!
!
!
!


subroutine get_surf_grad_fast(nd,npatches,norders,ixyzs,iptype,npts, &
  ffforminv,f,df)
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
  integer npts
  real *8 ffforminv(2,2,npts),f(nd,npts),df(nd,2,npts)
  real *8, allocatable :: fcoefs(:,:),dfuv(:,:)

  allocate(fcoefs(nd,npts),dfuv(nd,2,npts))
  call get_surf_uv_grad(nd,npatches,norders,ixyzs,iptype,npts,f,dfuv)

  return
end subroutine get_surf_grad_fast
!
!
!
!
!
subroutine get_surf_uv_grad(nd,npatches,norders,ixyzs,iptype,npts,f, &
   dfuv)
 
  implicit none
  integer nd,npatches,norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches),npts
  real *8 f(nd,npts),dfuv(nd,2,npts)

  !continue from here



