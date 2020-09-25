




subroutine l3d_sgradu(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(12),dpars(ndd)
  integer ipars(ndi)
  real *8 :: val
  complex *16 :: zk

  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  d = dx*targinfo(4) + dy*targinfo(5) + dz*targinfo(6)
  r=sqrt(dx**2+dy**2+dz**2)


  val =  -d/(r**3)

  return
end subroutine l3d_sgradu






subroutine l3d_sgradv(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(12),dpars(ndd)
  integer ipars(ndi)
  real *8 :: val
  complex *16 :: zk

  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  d = dx*targinfo(7) + dy*targinfo(8) + dz*targinfo(9)
  r=sqrt(dx**2+dy**2+dz**2)

  val =  -d/(r**3)

  return
end subroutine l3d_sgradv






subroutine l3d_dgradu(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: val
  complex *16 :: zk

  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  d = dx*srcinfo(4) + dy*srcinfo(5) + dz*srcinfo(6)
  r=sqrt(dx**2+dy**2+dz**2)


  val =  d/(r**3)

  return
end subroutine l3d_dgradu






subroutine l3d_dgradv(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: val
  complex *16 :: zk

  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  d = dx*srcinfo(7) + dy*srcinfo(8) + dz*srcinfo(9)
  r=sqrt(dx**2+dy**2+dz**2)

  val =  d/(r**3)

  return
end subroutine l3d_dgradv





subroutine l3d_spp_sum_dp(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val) 

!
!---------------------
!  This subroutine evaluates the difference kernel
!   (S'' + D')*4*pi (to be consistent with fmm scaling) 
!
!  Input arguments:
!  
!    - srcinfo: real *8 (12)
!        Source information
!    - ndt: integer
!        must be at least 12, with first twelve components 
!        being the standard targ info, xyz,dxyz/du,dxyz/dv,normal
!    - targinfo: real *8 (12)
!        target information
!    - ndd: integer
!        dpars not used
!    - dpars: real *8
!        dpars not used
!    - ndz: integer
!        zpars not used
!    - zpars: double complex 
!        zpars not used
!    - ndi: integer
!        ipars not used
!    - ipars: integer
!        ipars not used
!
!  Output arugments:
!    - val: double precision
!        s'' + d'
!    
!------------------

  implicit none
  integer ndd,ndi,ndz,ndt
  real *8 :: srcinfo(*), targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8, intent(out) :: val
  complex *16 :: zk

  real *8 dx(3),dns(3),dnt(3),r,d,drns,drnt,r3inv,r5inv
  real *8 dnsnt,rinv,over4pi
  data over4pi/0.07957747154594767d0/

  dx(1) = targinfo(1) - srcinfo(1)
  dx(2) = targinfo(2) - srcinfo(2)
  dx(3) = targinfo(3) - srcinfo(3)

  dns(1) = srcinfo(10)
  dns(2) = srcinfo(11)
  dns(3) = srcinfo(12)

  dnt(1) = targinfo(10)
  dnt(2) = targinfo(11)
  dnt(3) = targinfo(12)

  r = sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
  rinv = 1.0d0/r
  r3inv = rinv**3
  r5inv = rinv**5
  drns = dx(1)*dns(1) + dx(2)*dns(2) + dx(3)*dns(3)
  drnt = dx(1)*dnt(1) + dx(2)*dnt(2) + dx(3)*dnt(3)
  dnsnt = dnt(1)*dns(1) + dnt(2)*dns(2) + dnt(3)*dns(3)
  
!
!  second derivative of single layer
!
  val = 3*drnt**2*r5inv - r3inv

!
! add in derivative of double layer
!
  val = val - 3*drns*drnt*r5inv + dnsnt*r3inv
  val = val*over4pi

end subroutine l3d_spp_sum_dp


