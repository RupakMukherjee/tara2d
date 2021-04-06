
!===================================================================================
!=================== SUBROUTINE EULER ==============================================
!=================== DENSITY & VELOCITY  ===========================================

subroutine euler2coh(Nx,Ny,Nh,pi,time,nk,ukx,uky,nk_new,ukx_new,uky_new,dn_dt,dux_dt,duy_dt)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh, i, j
real ( kind = 8 ) pi,time,dt,Lx,Ly,mux,muy,ms,CS
complex ( kind = 8 ) kNLn(Nh,Ny),kNLx(Nh,Ny),kNLy(Nh,Ny)
complex ( kind = 8 ) nk(Nh,Ny),ukx(Nh,Ny),uky(Nh,Ny),nk_new(Nh,Ny),ukx_new(Nh,Ny),uky_new(Nh,Ny)
complex ( kind = 8 ) dn_dt(Nh,Ny),dux_dt(Nh,Ny),duy_dt(Nh,Ny)
common/comm/Lx,Ly,mux,muy,ms,CS,dt

do i = 1,Nh
  do j = 1,Ny
    kNLn(i,j) = dn_dt(i,j)
    nk_new(i,j) = nk(i,j) + kNLn(i,j)*dt     ! EULER SOLVER
    kNLx(i,j) = dux_dt(i,j)
    ukx_new(i,j) = ukx(i,j) + kNLx(i,j)*dt     ! EULER SOLVER
    kNLy(i,j) = duy_dt(i,j)
    uky_new(i,j) = uky(i,j) + kNLy(i,j)*dt     ! EULER SOLVER
  end do
end do

return

end subroutine euler2coh
