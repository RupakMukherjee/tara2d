
!===================================================================================
!=================== SUBROUTINE RUNGE - KUTTA 4 ====================================
!=================== DENSITY =======================================================

subroutine rk4_nk2coh(Nx,Ny,Nh,pi,time,nk,nk_new,ukx,uky,dn_dt,dux_dt,duy_dt)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh, i, j
real ( kind = 8 ) pi,time,dt,Lx,Ly,mux,muy,ms,CS
real ( kind = 8 ) n(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny)
complex ( kind = 8 ) kNLn1(Nh,Ny),kNLn2(Nh,Ny),kNLn3(Nh,Ny),kNLn4(Nh,Ny)
complex ( kind = 8 ) kNLx1(Nh,Ny),kNLx2(Nh,Ny),kNLx3(Nh,Ny),kNLx4(Nh,Ny)
complex ( kind = 8 ) kNLy1(Nh,Ny),kNLy2(Nh,Ny),kNLy3(Nh,Ny),kNLy4(Nh,Ny)
complex ( kind = 8 ) nk(Nh,Ny),ukx(Nh,Ny),uky(Nh,Ny),nk_new(Nh,Ny),ukx_new(Nh,Ny),uky_new(Nh,Ny)
complex ( kind = 8 ) dum_nk(Nh,Ny),dum_ukx(Nh,Ny),dum_uky(Nh,Ny)
complex ( kind = 8 ) dn_dt(Nh,Ny),dux_dt(Nh,Ny),duy_dt(Nh,Ny)
complex ( kind = 8 ) dn_dt_old(Nh,Ny),dux_dt_old(Nh,Ny),duy_dt_old(Nh,Ny)
complex ( kind = 8 ) dn_dt_new(Nh,Ny),dux_dt_new(Nh,Ny),duy_dt_new(Nh,Ny)
common/comm/Lx,Ly,mux,muy,ms,CS,dt

do i = 1,Nh
  do j = 1,Ny
    kNLn1(i,j) = dn_dt(i,j)
    dum_nk(i,j) = nk(i,j) + kNLn1(i,j)*dt/2.0d0
  end do
end do

  call derive2coh(Nx,Ny,Nh,pi,time+dt/2.0d0,dum_nk,ukx,uky, &
  dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new)

do i = 1,Nh
  do j = 1,Ny
    kNLn2(i,j) = dn_dt(i,j)
    dum_nk(i,j) = nk(i,j) + kNLn2(i,j)*dt/2.0d0
  end do
end do

  call derive2coh(Nx,Ny,Nh,pi,time+dt/2.0d0,dum_nk,ukx,uky, &
  dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new)

do i = 1,Nh
  do j = 1,Ny
    kNLn3(i,j) = dn_dt(i,j)
    dum_nk(i,j) = nk(i,j) + kNLn3(i,j)*dt/2.0d0
  end do
end do

  call derive2coh(Nx,Ny,Nh,pi,time+dt/2.0d0,dum_nk,ukx,uky, &
  dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new)

do i = 1,Nh
  do j = 1,Ny
    kNLn4(i,j) = dn_dt(i,j)
    nk_new(i,j) = nk(i,j) + dt/6.0d0*(kNLn1(i,j) + 2.0d0*kNLn2(i,j) + 2.0d0*kNLn3(i,j) + kNLn4(i,j))
  end do
end do

return

end subroutine rk4_nk2coh
