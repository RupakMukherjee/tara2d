
!===================================================================================
!=================== SUBROUTINE ADAMS BASHFORTH ====================================
!===================================================================================

subroutine ab2com(Nx,Ny,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k, &
              rho_k_new,rho_ux_k_new,rho_uy_k_new,E_k_new,Bx_k_new,By_k_new, &
              d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
              d_E_k_dt_old,d_E_k_dt_new,d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh, i, j, k
real ( kind = 8 ) pi,time,dt,Lx,Ly,spheat,mu,ms,CS,mu_0,eta,kf,A

complex ( kind = 8 ) rho_k(Nh,Ny),rho_ux_k(Nh,Ny),rho_uy_k(Nh,Ny),E_k(Nh,Ny),Bx_k(Nh,Ny),By_k(Nh,Ny)
complex ( kind = 8 ) rho_k_new(Nh,Ny),rho_ux_k_new(Nh,Ny),rho_uy_k_new(Nh,Ny),E_k_new(Nh,Ny),Bx_k_new(Nh,Ny),By_k_new(Nh,Ny)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny),d_rho_ux_k_dt_old(Nh,Ny),d_rho_uy_k_dt_old(Nh,Ny),d_E_k_dt_old(Nh,Ny)
complex ( kind = 8 ) d_Bx_k_dt_old(Nh,Ny),d_By_k_dt_old(Nh,Ny)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny),d_rho_ux_k_dt_new(Nh,Ny),d_rho_uy_k_dt_new(Nh,Ny),d_E_k_dt_new(Nh,Ny)
complex ( kind = 8 ) d_Bx_k_dt_new(Nh,Ny),d_By_k_dt_new(Nh,Ny)

common/comm/Lx,Ly,spheat,mu,ms,CS,mu_0,eta,dt,kf,A

!$OMP PARALLEL SHARED(dt,rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k),&
!$OMP & SHARED(rho_k_new,rho_ux_k_new,rho_uy_k_new,E_k_new,Bx_k_new,By_k_new),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new),&
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new),&
!$OMP & SHARED(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old),&
!$OMP & SHARED(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old) PRIVATE(i,j)
!$OMP DO

do i = 1,Nh
  do j = 1,Ny
    ! Density Equation Evolution.
    rho_k_new(i,j) = rho_k(i,j) + ( (3.0d0/2.0d0)*d_rho_k_dt_new(i,j) - (1.0d0/2.0d0)*d_rho_k_dt_old(i,j) )*dt
    ! Momentum Equation Evolution.
    rho_ux_k_new(i,j) = rho_ux_k(i,j) + ( (3.0d0/2.0d0)*d_rho_ux_k_dt_new(i,j) - (1.0d0/2.0d0)*d_rho_ux_k_dt_old(i,j) )*dt
    rho_uy_k_new(i,j) = rho_uy_k(i,j) + ( (3.0d0/2.0d0)*d_rho_uy_k_dt_new(i,j) - (1.0d0/2.0d0)*d_rho_uy_k_dt_old(i,j) )*dt
    ! Energy Equation Evolution.
    E_k_new(i,j) = E_k(i,j) + ( (3.0d0/2.0d0)*d_E_k_dt_new(i,j) - (1.0d0/2.0d0)*d_E_k_dt_old(i,j) )*dt
    ! Energy Equation Evolution.
    Bx_k_new(i,j) = Bx_k(i,j) + ( (3.0d0/2.0d0)*d_Bx_k_dt_new(i,j) - (1.0d0/2.0d0)*d_Bx_k_dt_old(i,j) )*dt
    By_k_new(i,j) = By_k(i,j) + ( (3.0d0/2.0d0)*d_By_k_dt_new(i,j) - (1.0d0/2.0d0)*d_By_k_dt_old(i,j) )*dt
  end do
end do

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine ab2com
