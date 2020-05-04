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
