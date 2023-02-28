      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 dpars(2)
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      complex *16, allocatable :: rhs(:)
      real *8, allocatable :: rrhs(:),pot(:),potex(:)
      real *8, allocatable :: sigma1(:),sigma2(:),sigma3(:)
      real *8, allocatable :: errs1(:),errs2(:),errs3(:)
      
      real *8, allocatable :: smat(:,:),sgradumat(:,:),sgradvmat(:,:)
      real *8, allocatable :: spmat(:,:),dmat(:,:),diffmat(:,:)
      real *8, allocatable :: dgradumat(:,:),dgradvmat(:,:)
      real *8, allocatable :: xmat_rep1(:,:),xmat_rep2(:,:)
      real *8, allocatable :: xmat_rep3(:,:)
      real *8, allocatable :: xmat1(:,:),xmat2(:,:),xmat3(:,:)
      real *8, allocatable :: xmat4(:,:),xmat5(:,:),xmat6(:,:)
      real *8, allocatable :: xmatu(:,:),xmatv(:,:)
      real *8, allocatable :: ffforminv(:,:,:),gginv(:),gg(:)
      complex *16, allocatable :: coefs(:,:),coefs_pot(:,:)
      real *8, allocatable :: dumat(:,:),dvmat(:,:),xutmp(:,:),
     1   xvtmp(:,:),umat(:,:),uvs(:,:),meancrvs(:)
      real *8, allocatable :: pols(:),ders(:,:)
      
      real *8, allocatable :: w(:,:),work(:)



      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname
      real *8, allocatable :: s1(:),s2(:),s3(:)


      logical isout0,isout1

      complex *16 ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


      
      igeomtype = 1
      ipars(1) = 2
      npatches=12*(4**ipars(1))

      norder = 3 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      print *, 'npts=',npts

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      allocate(rhs(npts),rrhs(npts),pot(npts),potex(npts))

c
c       define rhs to be one of the ynm's
c
      nmax = 2
      allocate(coefs(0:nmax,0:nmax),coefs_pot(0:nmax,0:nmax))
      coefs = 0
      coefs_pot = 0

      fname='rep-comp-res/boundary-data.dat'
      open(unit=33,file=trim(fname))

 1121 format(2(2x,i3),2(2x,e22.16))
      do i=0,nmax
        do j=0,i
          coefs(i,j) = (hkrand(0)-0.5d0) + ima*(hkrand(0)-0.5d0)
          if(i.eq.0.and.j.eq.0) coefs(i,j) = 0
          coefs_pot(i,j) = -coefs(i,j)*i*(i+1.0d0)/(2*i+1.0d0)**2
          write(33,1121) i,j,coefs(i,j),coefs_pot(i,j)
        enddo
      enddo

      close(33)

      call prin2('coefs=*',coefs,(nmax+1)**2)
      call prin2('coefs_pot=*',coefs_pot,(nmax+1)**2)

      allocate(w(0:nmax,0:nmax))
      call l3dget_multynm_rhs(nmax,coefs,12,srcvals,rrhs,npts,w)
      call l3dget_multynm_rhs(nmax,coefs_pot,12,srcvals,potex,npts,w)
      call prin2('rrhs=*',rrhs,12)
      call prin2('potex=*',potex,12)
      
      write(fname,'(a,i3.3,a)') 'rep-comp-res/ressum-',npatches,
     1  '.dat'
      open(unit=34,file=trim(fname))

      write(34,'(2x,i5)') npatches
      write(34,'(2x,i5)') norder
      write(34,'(2x,i5)') npts


      
      allocate(smat(npts,npts),sgradumat(npts,npts))
      allocate(sgradvmat(npts,npts),dgradumat(npts,npts))
      allocate(dgradvmat(npts,npts),xmat_rep1(npts,npts))
      allocate(xmat_rep2(npts,npts))
      allocate(xmat1(npts,npts),xmat2(npts,npts),xmat5(npts,npts))
      allocate(xmatu(npts,npts),xmatv(npts,npts))

      eps = 0.51d-9

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)

      allocate(spmat(npts,npts),dmat(npts,npts),diffmat(npts,npts))

      call get_lapbel_matrices(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,ipatch_id,uvs_targ,eps,smat,sgradumat,
     2  sgradvmat,dgradumat,dgradvmat,spmat,dmat,diffmat)

      alpha = 1.0d0
      beta = 0.0d0
      allocate(pols(npols),uvs(2,npols),umat(npols,npols),
     1   ders(2,npols),ffforminv(2,2,npts))
      allocate(dumat(npols,npols),dvmat(npols,npols))
      allocate(xutmp(npols,npols),xvtmp(npols,npols))
      
      call get_inv_first_fundamental_form(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ffforminv)
      call prinf('done computing ffforminv*',i,0)


      allocate(gg(npts),gginv(npts))

      do i=1,npts
        gginv(i) = sqrt(ffforminv(1,1,i)*ffforminv(2,2,i)-
     1    ffforminv(1,2,i)*ffforminv(2,1,i))
        gg(i) = 1.0d0/gginv(i)
      enddo

      call prin2('gginv=*',gginv,12)


      do j=1,npts
        do i=1,npts
          xmatu(i,j) = (sgradumat(i,j)*ffforminv(1,1,i)+
     1      sgradvmat(i,j)*ffforminv(1,2,i))
          xmatv(i,j) = (sgradumat(i,j)*ffforminv(2,1,i)+
     1      sgradvmat(i,j)*ffforminv(2,2,i))
        enddo
      enddo


c
c
c    add in contribution of smat*w*smat
c 
      do j=1,npts
        do i=1,npts
           xmat1(i,j) = wts(j)/4/pi
        enddo
      enddo

      call dgemm('n','n',npts,npts,npts,alpha,xmat1,npts,smat,npts,
     1  beta,xmat2,npts)

      
      call dgemm('n','n',npts,npts,npts,alpha,smat,npts,xmat2,npts,
     1  beta,xmat5,npts)
      
       
c
c
c        now start constructing matrix for representation \cA_{1}
c
c
      call dgemm('n','n',npts,npts,npts,alpha,dgradumat,npts,xmatu,npts,
     1  beta,xmat_rep1,npts)
      call dgemm('n','n',npts,npts,npts,alpha,dgradvmat,npts,xmatv,npts,
     1  alpha,xmat_rep1,npts)
      do j=1,npts
        do i=1,npts 
          xmat_rep1(i,j) = -xmat_rep1(i,j) + xmat5(i,j)
        enddo
      enddo
      

c
c   test S(\Delta + W) S with representation \cA_{1}
c
      call dgemv('n',npts,npts,alpha,xmat_rep1,npts,rrhs,1,beta,pot,1)
      
      erra = 0
      ra = 0
      do i=1,npts
        erra = erra + abs(potex(i)-pot(i))**2*wts(i)
        ra = ra + abs(potex(i))**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      err_app_rep1 = erra
      call prin2('error in final operator apply for rep1=*',erra,1)

c
c   start generating matrix for rep \cA_{2}
c
      allocate(meancrvs(npts))
      call get_mean_curvature(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,meancrvs)
      do j=1,npts
        do i=1,npts
          xmat1(i,j) = spmat(i,j)*meancrvs(i)
        enddo
      enddo

      call prinf('done getting mean curvature*',i,0)

      call dgemm('n','n',npts,npts,npts,alpha,smat,npts,xmat1,npts,
     1  beta,xmat_rep2,npts)
      
      rtmp = 2.0d0
      call dgemm('n','n',npts,npts,npts,alpha,smat,npts,diffmat,npts,
     1  rtmp,xmat_rep2,npts)
      
      rtmp = -1.0d0
      call dgemm('n','n',npts,npts,npts,alpha,dmat,npts,dmat,npts,
     1  rtmp,xmat_rep2,npts)
      do i=1,npts
        xmat_rep2(i,i) = -0.25d0 + xmat_rep2(i,i)
      enddo

      do j=1,npts
        do i=1,npts
          xmat_rep2(i,j) = xmat_rep2(i,j) + xmat5(i,j)
        enddo
      enddo
        
      
      
c
c   test S(\Delta + W) S with representation \cA_{2}
c
      call dgemv('n',npts,npts,alpha,xmat_rep2,npts,rrhs,1,beta,pot,1)
      
      erra = 0
      ra = 0
      do i=1,npts
        erra = erra + abs(potex(i)-pot(i))**2*wts(i)
        ra = ra + abs(potex(i))**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      err_app_rep2 = erra
      call prin2('error in final operator apply for rep2=*',erra,1)

 1131 format(2(2x,e11.5))
 1141 format(2(2x,i4))
      write(34,*) " ==="
      write(34,*) " Error in mat apply"
      write(34,1131) err_app_rep1,err_app_rep2
      write(34,*) "===="


c
c
c   solve problem via gmres
c
      numit = 100
      allocate(errs1(numit+1),errs2(numit+1),errs3(numit+1))
      allocate(sigma1(npts),sigma2(npts),sigma3(npts))

      niter1 = 0
      niter2 = 0
      call dgmres_slow(npts,xmat_rep1,potex,sigma1,numit,niter1,errs1)
      call dgmres_slow(npts,xmat_rep2,potex,sigma2,numit,niter2,errs2)

      call prinf('niter1=*',niter1,1)
      call prin2('errs1=*',errs1,niter1)
      call prinf('niter2=*',niter2,1)
      call prin2('errs2=*',errs2,niter2)
      

      erra1 = 0
      erra2 = 0
      ra = 0

      do i=1,npts
        ra = ra + abs(rrhs(i))**2*wts(i)
        erra1 = erra1 + abs(rrhs(i)-sigma1(i))**2*wts(i)
        erra2 = erra2 + abs(rrhs(i)-sigma2(i))**2*wts(i)
      enddo
      erra1 = sqrt(erra1/ra)
      erra2 = sqrt(erra2/ra)
      
      call prin2('error in solve rep1=*',erra1,1)
      call prin2('error in solve rep2=*',erra2,1)

      write(34,*) "error in solve"
      write(34,1131) erra1,erra2
      write(34,*) "===="
      write(34,*) "niter"
      write(34,1141) niter1,niter2,niter3
      write(34,*) "===="
      write(34,*) "errors as a function of iteration number"
      do i=1,niter1
        write(34,'(2x,e11.5)') errs1(i)
      enddo
      do i=1,niter2
        write(34,'(2x,e11.5)') errs2(i)
      enddo
      write(34,*) '===='
      write(34,*) 'Singular values'


      


c
c    compute the spectra of these matrices
c
      allocate(s1(npts),s2(npts))
      lw = 100*npts
      allocate(work(lw))
      call dgesvd('n','n',npts,npts,xmat_rep1,npts,s1,xmat1,npts,
     1   xmat2,npts,work,lw,info)
      call dgesvd('n','n',npts,npts,xmat_rep2,npts,s2,xmat1,npts,
     1   xmat2,npts,work,lw,info)

      do i=1,npts
        write(34,1131) s1(i),s2(i)
      enddo
      close(34)


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


      external xtri_stell_eval,xtri_sphere_eval
      
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
        vmin = 0
        vmax = 2*pi
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
          call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

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



   



      subroutine l3dget_multynm_rhs(nmax,coefs,ndx,xyzs,ynms,npts,w)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      real *8 ynms(npts)
      complex *16 coefs(0:nmax,0:nmax),ztmp
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
        ynms(i) = 0
        do j=0,nmax
          do l=0,j
            ztmp = coefs(j,l)*ynm(j,l)*exp(ima*l*phi) 
            ynms(i) = ynms(i) + real(ztmp) 
          enddo
        enddo
      enddo
       
      return
      end
c
c
c
c
c
      subroutine dgmres_slow(npts,xmat,rhs,soln,numit,niter,errs)
      implicit real *8(a-h,o-z)
      real *8 xmat(npts,npts),rhs(npts),soln(npts),errs(numit+1)

      real *8, allocatable :: vmat(:,:),hmat(:,:)
      real *8, allocatable :: cs(:),sn(:)
      real *8, allocatable :: svec(:),yvec(:),wtmp(:)

      allocate(vmat(npts,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(npts),svec(numit+1),yvec(numit+1))


      niter=0
      alpha = 1.0d0
      beta = 0.0d0
      eps_gmres = 1.0d-14

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
        rb = rb + abs(rhs(i))**2
      enddo
      rb = sqrt(rb)

      do i=1,npts
        vmat(i,1) = rhs(i)/rb
      enddo

      svec(1) = rb


      do it=1,numit
        it1 = it + 1

        call dgemv('n',npts,npts,alpha,xmat,npts,vmat(1,it),1,beta,
     1      wtmp,1)

        do k=1,it
          hmat(k,it) = 0
          do j=1,npts
            hmat(k,it) = hmat(k,it) + wtmp(j)*vmat(j,k)
          enddo

          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
        enddo
          
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


          call dgemv('n',npts,npts,alpha,xmat,npts,soln,1,beta,
     1      wtmp,1)

          do i=1,npts
            rres = rres + abs(wtmp(i)-rhs(i))**2
          enddo
          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo

c


      return
      end
      







