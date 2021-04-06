!===================================================================================
!=================== SUBROUTINE ADAMS BASHFORTH ====================================
!=================== DENSITY & VELOCITY  ===========================================

subroutine ab2coh(Nx,Ny,Nh,pi,time,nk,ukx,uky,nk_new,ukx_new,uky_new, &
dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh, i, j
real ( kind = 8 ) pi,time,dt,Lx,Ly,mux,muy,ms,CS
complex ( kind = 8 ) nk(Nh,Ny),ukx(Nh,Ny),uky(Nh,Ny),nk_new(Nh,Ny),ukx_new(Nh,Ny),uky_new(Nh,Ny)
complex ( kind = 8 ) dn_dt_old(Nh,Ny),dux_dt_old(Nh,Ny),duy_dt_old(Nh,Ny)
complex ( kind = 8 ) dn_dt_new(Nh,Ny),dux_dt_new(Nh,Ny),duy_dt_new(Nh,Ny)
common/comm/Lx,Ly,mux,muy,ms,CS,dt

!$OMP PARALLEL SHARED(dt,nk,ukx,uky,nk_new,ukx_new,uky_new,dn_dt_new,dux_dt_new,duy_dt_new),&
!$OMP & SHARED(dn_dt_old,dux_dt_old,duy_dt_old) PRIVATE(i,j)
!$OMP DO

do i = 1,Nh
  do j = 1,Ny
    ! Density Equation Evolution.
    nk_new(i,j) = nk(i,j) + ( (3.0d0/2.0d0)*dn_dt_new(i,j) - (1.0d0/2.0d0)*dn_dt_old(i,j) )*dt
    ! Momentum Equation Evolution.
    ukx_new(i,j) = ukx(i,j) + ( (3.0d0/2.0d0)*dux_dt_new(i,j) - (1.0d0/2.0d0)*dux_dt_old(i,j) )*dt
    uky_new(i,j) = uky(i,j) + ( (3.0d0/2.0d0)*duy_dt_new(i,j) - (1.0d0/2.0d0)*duy_dt_old(i,j) )*dt
  end do
end do

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine ab2coh
