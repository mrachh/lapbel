!
!      surf_vtk_plot_vec - generate a vtk file to plot the surface,
!         along with prescribed vector field
!
!


subroutine surf_vtk_plot_vec(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sigma,fname,title)
!
!   This subroutine writes a vtk to plot the surface along
!   with a vector field. Currently only supports triangular patches
!
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs
!f2py intent(in) srcvals,sigma,fname,title
!
  implicit none
  integer npatches,norders(npatches),ixyzs(npatches+1),npts
  integer iptype(npatches)
  real *8 srccoefs(9,npts),srcvals(12,npts),sigma(3,npts)
  real *8, allocatable :: sigma_coefs(:,:),pols(:),rtmp(:)
  character (len=*) fname,title

  real *8, allocatable :: xyzs(:,:),uvs(:,:,:),splot(:,:)
  integer, allocatable :: kovers(:),nps(:),ipstart(:)

  integer i,j,k,l,ipatch,npout,kover,npols
  integer itrip,itric1,nb,nlmax,nuv,istart,iend,nd
  integer ilstart,itri,iunit1,m,ncell,ncsize,norder,nuvl,i1
  integer idim

  real *8 ra,erra

!
!  get the coefs of the density
!
   nd = 3
   allocate(sigma_coefs(nd,npts))
   call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts, &
     sigma,sigma_coefs)


 
!
!   estimate kovers, nps 
!

  npout = 0
  allocate(kovers(npatches),nps(npatches),ipstart(npatches+1))
  do i=1,npatches
    npols = ixyzs(i+1)-ixyzs(i)
    kover = 0
    if(npols.gt.4**0) kover = 1
    if(npols.gt.4**1) kover = 2
    if(npols.gt.4**2) kover = 3
    if(npols.gt.4**3) kover = 4
    kovers(i) = kover
    nps(i) = 4**kover
    if(iptype(i).eq.1) nps(i) = nps(i)*3
    npout = npout + nps(i) 
  enddo


  ipstart(1) = 1
  nps(1) = nps(1) + 1
  call cumsum(npatches,nps,ipstart(2))
  nps(1) = nps(1) - 1

  allocate(xyzs(3,npout),splot(nd,npout))

!
!   get uvs of all patches of type = 1
!
  
  nlmax = 4
  nuv = (4**(nlmax+1)-1)/3
  allocate(uvs(2,3,nuv))

  uvs(1,1,1) = 0
  uvs(2,1,1) = 0
  uvs(1,2,1) = 1
  uvs(2,2,1) = 0
  uvs(1,3,1) = 0
  uvs(2,3,1) = 1

  do i=0,nlmax-1
    istart = (4**(i)-1)/3+1
    nb = 4**i
    iend = istart + nb-1
    do itrip = istart,iend
      itric1 = (itrip-istart)*4 + iend
      call gettrichildren(uvs(1,1,itrip),uvs(1,1,itric1+1), &
       uvs(1,1,itric1+2),uvs(1,1,itric1+3),uvs(1,1,itric1+4))   
    enddo
  enddo


  do ipatch=1,npatches
    istart = ipstart(ipatch)
    npols = ixyzs(ipatch+1)-ixyzs(ipatch)
    norder = norders(ipatch)
    allocate(pols(npols))
    if(iptype(ipatch).eq.1) then

      nuvl = ipstart(ipatch+1)-ipstart(ipatch)
      ilstart = 4**(kovers(ipatch)-1)/3+1
      nb = 4**(kovers(ipatch))
      do i=1,nb
        itri = i+ilstart-1
        do j=1,3
          call koorn_pols(uvs(1,j,itri),norder,npols,pols)
          
          do m=1,3
            xyzs(m,istart+3*(i-1)+j-1) = 0
          enddo
          do idim=1,nd
            splot(idim,istart+3*(i-1)+j-1) = 0
          enddo

          do l=1,npols
            do m=1,3
              xyzs(m,istart+3*(i-1)+j-1) = & 
                xyzs(m,istart+3*(i-1)+j-1) + &
                pols(l)*srccoefs(m,ixyzs(ipatch)+l-1)
            enddo
            do idim=1,nd
              splot(idim,istart+3*(i-1)+j-1) = &
               splot(idim,istart+3*(i-1)+j-1)+ &
               pols(l)*sigma_coefs(idim,ixyzs(ipatch)+l-1)
            enddo
          enddo
        enddo
      enddo
    endif
    deallocate(pols)
  enddo
  
  iunit1 = 877
  open(unit = iunit1, file=trim(fname))

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') trim(title)
  write(iunit1,'(a)') "ASCII"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i9,a)') "POINTS ", npout, " float"

  do i = 1,npout
    write(iunit1,"(E11.5,2(2x,e11.5))") xyzs(1,i), xyzs(2,i), xyzs(3,i)
  end do

  ncell = 0
  ncsize = 0
  do i=1,npatches
    ncell = ncell + 4**kovers(i)
    if(iptype(i).eq.1) ncsize = ncsize + 4*(4**kovers(i))
  enddo

  write(iunit1,'(a,i9,i9)') "CELLS ", ncell, ncsize

  do ipatch=1,npatches
    nb = 4**kovers(ipatch)
    if(iptype(ipatch).eq.1) then
      istart = ipstart(ipatch) 
      do i = 1,nb
        i1 = istart + 3*(i-1) 
        write(iunit1,'(a,i9,i9,i9)') "3 ", i1-1, i1, i1+1
      enddo
    endif
  end do

  write(iunit1,'(a,i9)') "CELL_TYPES ", ncell
  do ipatch = 1,npatches
    nb = 4**kovers(ipatch)
    if(iptype(ipatch).eq.1) then
      do i=1,nb
        write(iunit1,'(a)') "5"
      enddo
    endif
  end do

  write(iunit1,'(a)') ""
  write(iunit1,'(a,i9)') "POINT_DATA ", npout
  write(iunit1,'(a,i4)') "SCALARS vec_comps float ", 3
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,npout
    write(iunit1,'(E11.5,2x,E11.5,2x,e11.5)') &
      splot(1,i),splot(2,i),splot(3,i)
  end do

  close(iunit1)



end subroutine surf_vtk_plot_vec
!
!
