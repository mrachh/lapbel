!
!  This file contains the following user callable routines:
!    
!    form_surf_lap_mat - form the matrix for applying the
!      surface laplacian operator
!
!
!
!
!
subroutine form_surf_lap_mat(npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,xmat)
!
!-----------------------------
!  Compute the surface laplacian matrix  
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - xmat: double precision(npts,npts)
!        matrix for applying surface laplacian
!        
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: xmat(npts,npts)
  real *8, allocatable :: ffforminv(:,:,:),xtmp(:,:)

  integer i,istart,npols,j,l

  allocate(ffforminv(2,2,npts))

  do i=1,npts
    do j=1,npts
      xmat(j,i) = 0
    enddo
  enddo


  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  do i=1,npatches
    istart = ixyzs(i)
    npols = ixyzs(i+1)-ixyzs(i)
    allocate(xtmp(npols,npols))

    if(iptype(i).eq.1) then
      call form_surf_lap_mat_tri(norders(i),npols, &
        ffforminv(1,1,istart),xtmp)
    endif
    do j=1,npols
      do l=1,npols
        xmat(istart+l-1,istart+j-1) = xtmp(l,j)
      enddo
    enddo
    deallocate(xtmp)

  enddo

  return
end subroutine form_surf_lap_mat
!
!
!
!
!
!
!
!
!
subroutine form_surf_lap_mat_tri(norder,n,ffforminv,xmat)
!
!-----------------------------
!  Compute the surface laplacian matrix for a triangle
!
!  Input arguments:
!
!    - norder: integer
!        order of discretization 
!    - n: integer 
!        number of points on patch 
!    - ffforminv: real *8 (2,2,n)
!        inverse of first fundamental form
!
!  Output arguments:
!
!    - xmat: double precision(n,n)
!        matrix for applying surface laplacian
!        
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: norder,n
  real *8, intent(in) :: ffforminv(2,2,n) 
  real *8, intent(out) :: xmat(n,n)
  real *8, allocatable :: umat(:,:),dumat(:,:),xutmp(:,:)
  real *8, allocatable :: xvtmp(:,:),dvmat(:,:)
  real *8, allocatable :: xutmp2(:,:),xvtmp2(:,:)
  real *8, allocatable :: xutmp3(:,:),xvtmp3(:,:)
  real *8 uvs(2,n),pols(n)
  real *8 ders(2,n)
  integer i,istart,npols,j,l
  real *8 done,dzero,gg(n),gginv(n)

  allocate(umat(n,n),dumat(n,n),dvmat(n,n),xutmp(n,n),xvtmp(n,n))
  allocate(xutmp2(n,n),xvtmp2(n,n),xutmp3(n,n),xvtmp3(n,n))

  do i=1,n
    gginv(i) = sqrt(ffforminv(1,1,i)*ffforminv(2,2,i)- &
       ffforminv(1,2,i)*ffforminv(2,1,i))
    gg(i) = 1.0d0/gginv(i) 
  enddo

  
  call get_vioreanu_nodes(norder,n,uvs)

  call koorn_vals2coefs(norder,n,uvs,umat)

  do i=1,n
    call koorn_ders(uvs(1,i),norder,n,pols,ders)
    do j=1,n
      dumat(j,i) = ders(1,j)
      dvmat(j,i) = ders(2,j)
    enddo
  enddo

  done =1
  dzero = 0
  call dgemm('t','n',n,n,n,done,dumat,n,umat,n,dzero,xutmp,n)
  call dgemm('t','n',n,n,n,done,dvmat,n,umat,n,dzero,xvtmp,n)

!
!    Now apply the diagonal scaling 
!
  do j=1,n
    do i=1,n
      xutmp2(i,j) = gg(i)*(ffforminv(1,1,i)*xutmp(i,j) + &
        ffforminv(1,2,i)*xvtmp(i,j))
      xvtmp2(i,j) = gg(i)*(ffforminv(2,1,i)*xutmp(i,j) + &
        ffforminv(2,2,i)*xvtmp(i,j))
    enddo
  enddo

  call dgemm('n','n',n,n,n,done,xutmp,n,xutmp2,n,dzero,xutmp3,n)
  call dgemm('n','n',n,n,n,done,xvtmp,n,xvtmp2,n,dzero,xvtmp3,n)

  do j=1,n
    do i=1,n
      xmat(i,j) = (xutmp3(i,j) + xvtmp3(i,j))*gginv(i)
    enddo
  enddo


  return
end subroutine form_surf_lap_mat_tri
!
!
!
!
!
!
!
!
!
subroutine get_second_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sfform)
!
!------------------------
!  This subroutine computes the second fundamental form at
!  the discretization nodes on the surface.
!
!  The second fundamental form is
!  
!  .. math::
!    
!    \begin{bmatrix} x_{uu} \cdot n & x_{uv} \cdot n \\
!    x_{uv} \cdot n & x_{vv} \cdot n \end{bmatrix}
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - sfform: double precision(2,2,npts)
!        second fundamental form at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: sfform(2,2,npts)
  integer i
  real *8 L,M,N
  real *8, allocatable :: dxuv(:,:)
  real *8, allocatable :: dxuv2(:,:,:)


  allocate(dxuv(6,npts))
! Calculating x_{uu}, x_{uv}, x_{uv} stored in dxuv2 
 
  do i=1,npts
    dxuv(1,i) = srcvals(4,i)  
    dxuv(2,i) = srcvals(5,i)  
    dxuv(3,i) = srcvals(6,i)  
    dxuv(4,i) = srcvals(7,i)  
    dxuv(5,i) = srcvals(8,i)  
    dxuv(6,i) = srcvals(9,i)   
  enddo

  allocate(dxuv2(6,2,npts))

  call get_surf_uv_grad(6,npatches,norders,ixyzs,iptype,npts,dxuv,dxuv2)



  do i=1,npts
    L = dxuv2(1,1,i)*srcvals(10,i) + dxuv2(2,1,i)*srcvals(11,i) + dxuv2(3,1,i)*srcvals(12,i) ! Calculation of L, M, N. L = x_uu \cdot n 
    M = dxuv2(1,2,i)*srcvals(10,i) + dxuv2(2,2,i)*srcvals(11,i) + dxuv2(3,2,i)*srcvals(12,i)  ! M = x_uv \cdot n
    N = dxuv2(4,2,i)*srcvals(10,i) + dxuv2(5,2,i)*srcvals(11,i) + dxuv2(6,2,i)*srcvals(12,i)  ! N = x_vv \cdot n
    sfform(1,1,i) = L
    sfform(2,1,i) = M
    sfform(1,2,i) = M
    sfform(2,2,i) = N
  enddo

  return
end subroutine get_second_fundamental_form
!
!
!
subroutine get_mean_curvature(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,mean_curv)
!
!------------------------
!  This subroutine computes the mean curvature at
!  the discretization nodes on the surface.
!
!  The mean curvature is
!  
!  .. math::
!    
!    0.5*Trace(II \cdot I^{-1}) \\
!    
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - mean_curv: double precision(npts)
!        mean curvature at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: mean_curv(npts)
  integer i
  real *8, allocatable :: ffform(:,:,:),ffforminv(:,:,:)
  real *8, allocatable :: sfform(:,:,:)

  allocate(ffform(2,2,npts))

  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)



  allocate(ffforminv(2,2,npts))

  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  allocate(sfform(2,2,npts))

  call get_second_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sfform)
  
  print *,"Point on surface:", srcvals(1,3),srcvals(2,3), srcvals(3,3) 
  print *,"First fundamental form=", ffform(:,:, 3) 
  print *,"Inverse first fundamental form=", ffforminv(:,:, 3)
  print *,"Second fundamental form=", sfform(:,:, 3)
 


! Calculating mean curvature 
 
  do i=1,npts
    mean_curv(i) = -0.5*(sfform(1,1,i)*ffforminv(1,1,i) + &
                     sfform(1,2,i)*ffforminv(2,1,i) + &
                     sfform(2,1,i)*ffforminv(1,2,i) + &
                     sfform(2,2,i)*ffforminv(2,2,i))
  enddo
  print *,"Mean=", mean_curv(3)
 
  return
end subroutine get_mean_curvature

