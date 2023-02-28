      subroutine get_lapbel_matrices(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ipatch_id,uvs_targ,eps,
     2   smat,sgradumat,sgradvmat,dgradumat,dgradvmat,spmat,dmat,
     3   diffmat)
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: npatches,norders(npatches),npts
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: ipatch_id(npts)
      real *8, intent(in) :: uvs_targ(2,npts)
      real *8, intent(out) :: smat(npts,npts),sgradumat(npts,npts)
      real *8, intent(out) :: sgradvmat(npts,npts)
      real *8, intent(out) :: dgradumat(npts,npts),dgradvmat(npts,npts)
      real *8, intent(out) :: spmat(npts,npts),dmat(npts,npts)
      real *8, intent(out) :: diffmat(npts,npts)
      
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: xmattmp(:,:)

      integer iquadtype
      real *8 rfac0
      integer ipars
      real *8 dpars
      complex *16 zpars
      procedure (), pointer :: fker

      external l3d_slp,l3d_sgradu,l3d_sgradv,l3d_dgradu,l3d_dgradv
      external l3d_dlp,l3d_sprime,l3d_spp_sum_dp

      done = 1
      pi = atan(done)*4
      over4pi = 1.0d0/4/pi

      nnz = npatches*npts
      print *, nnz
      allocate(row_ptr(npts+1),col_ind(nnz),iquad(nnz+1))
      do i=1,npts+1
        row_ptr(i) = (i-1)*npatches+1
      enddo


      do i=1,npts
        do j=1,npatches
          ii = (i-1)*npatches+j
          col_ind(ii) = j
        enddo
      enddo
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,iquad)
      nquad = iquad(nnz+1)-1

      print *, "nquad=",nquad
      print *, "npts*npts=",npts*npts
      print *, "error=",npts*npts-nquad

      allocate(xmattmp(npts,npts))
      

      ndd = 0
      ndi = 0
      ndz = 0
      rfac0 = 1.25d0
      ndtarg = 12
      ipv = 1
      fker => l3d_slp

      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,rfac0,nquad,xmattmp)
      do i=1,npts
        do j=1,npts
          smat(j,i) = xmattmp(i,j)
        enddo
      enddo

      print *, "done with first matrix"

      fker => l3d_sgradu
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,rfac0,nquad,xmattmp)
      do i=1,npts
        do j=1,npts
          sgradumat(j,i) = xmattmp(i,j)
        enddo
      enddo

      print *, "done with second matrix"

      fker => l3d_sgradv
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,rfac0,nquad,xmattmp)
      do i=1,npts
        do j=1,npts
          sgradvmat(j,i) = xmattmp(i,j)
        enddo
      enddo

      print *, "done with third matrix"

      fker => l3d_dgradu
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,rfac0,nquad,xmattmp)
      do i=1,npts
        do j=1,npts
          dgradumat(j,i) = xmattmp(i,j)
        enddo
      enddo
      print *, "done with fourth matrix"


      fker => l3d_dgradv
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,rfac0,nquad,xmattmp)
      do i=1,npts
        do j=1,npts
          dgradvmat(j,i) = xmattmp(i,j)
        enddo
      enddo

      print *, "done with fifth matrix"


      fker => l3d_dlp
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,rfac0,nquad,xmattmp)
      do i=1,npts
        do j=1,npts
          dmat(j,i) = xmattmp(i,j)
        enddo
      enddo

      print *, "done with sixth matrix"


      fker => l3d_sprime
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,rfac0,nquad,xmattmp)
      do i=1,npts
        do j=1,npts
          spmat(j,i) = xmattmp(i,j)
        enddo
      enddo

      print *, "done with seventh matrix"


      fker => l3d_spp_sum_dp
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,rfac0,nquad,xmattmp)
      do i=1,npts
        do j=1,npts
          diffmat(j,i) = xmattmp(i,j)
        enddo
      enddo

      print *, "done with eighth matrix"

      return
      end
