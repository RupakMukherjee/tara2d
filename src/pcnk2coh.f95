
!===================================================================================
!=================== SUBROUTINE PREDICTOR - CORRECTOR ==============================
!=================== DENSITY =======================================================

subroutine pc_nk2coh(Nx,Ny,Nh,pi,time,nk,nk_new,ukx,uky,dn_dt,dux_dt,duy_dt)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh, i, j
real ( kind = 8 ) pi,time,dt,Lx,Ly,mux,muy,ms,CS
real ( kind = 8 ) n(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny)
complex ( kind = 8 ) predict(Nh,Ny),correct(Nh,Ny)
complex ( kind = 8 ) nk(Nh,Ny),ukx(Nh,Ny),uky(Nh,Ny),nk_new(Nh,Ny),ukx_new(Nh,Ny),uky_new(Nh,Ny)
complex ( kind = 8 ) dum_nk(Nh,Ny),dum_ukx(Nh,Ny),dum_uky(Nh,Ny)
complex ( kind = 8 ) dn_dt(Nh,Ny),dux_dt(Nh,Ny),duy_dt(Nh,Ny)
complex ( kind = 8 ) dn_dt_old(Nh,Ny),dux_dt_old(Nh,Ny),duy_dt_old(Nh,Ny)
complex ( kind = 8 ) dn_dt_new(Nh,Ny),dux_dt_new(Nh,Ny),duy_dt_new(Nh,Ny)
common/comm/Lx,Ly,mux,muy,ms,CS,dt

do i = 1,Nh
  do j = 1,Ny
    predict(i,j) = dn_dt(i,j)
    dum_nk(i,j) = nk(i,j) + predict(i,j)*dt   ! Predictor using Euler Method...
  end do
end do

  call derive2coh(Nx,Ny,Nh,pi,time+dt,dum_nk,ukx,uky, &
  dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new)

do i = 1,Nh
  do j = 1,Ny
    correct(i,j) = dn_dt(i,j)
    nk_new(i,j) = nk(i,j) + ( predict(i,j) + correct(i,j) )*dt/2.0d0   ! Corrector using Predicted slope...
  end do
end do

return

end subroutine pc_nk2coh
