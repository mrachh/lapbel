




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



