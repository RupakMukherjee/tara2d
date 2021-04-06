
!===================================================================================
!=================== SUBROUTINE ADAMS BASHFORTH ====================================
!===================================================================================

subroutine ab3com(Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
              rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new, &
              d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
              d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
              d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nz,Nh,iret,thread_num, i, j, k
real ( kind = 8 ) pi,time,dt,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,kf,time_max

complex ( kind = 8 ) rho_k(Nh,Ny,Nz),rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) E_k(Nh,Ny,Nz),Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_new(Nh,Ny,Nz),rho_ux_k_new(Nh,Ny,Nz),rho_uy_k_new(Nh,Ny,Nz),rho_uz_k_new(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_new(Nh,Ny,Nz),Bx_k_new(Nh,Ny,Nz),By_k_new(Nh,Ny,Nz),Bz_k_new(Nh,Ny,Nz)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_ux_k_dt_old(Nh,Ny,Nz),d_rho_uy_k_dt_old(Nh,Ny,Nz),d_rho_uz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_Bx_k_dt_old(Nh,Ny,Nz),d_By_k_dt_old(Nh,Ny,Nz),d_Bz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_ux_k_dt_new(Nh,Ny,Nz),d_rho_uy_k_dt_new(Nh,Ny,Nz),d_rho_uz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_Bx_k_dt_new(Nh,Ny,Nz),d_By_k_dt_new(Nh,Ny,Nz),d_Bz_k_dt_new(Nh,Ny,Nz)

common/comm/iret,thread_num,time_max,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,dt,kf

!$OMP PARALLEL SHARED(dt,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),&
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new),&
!$OMP & SHARED(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old),&
!$OMP & SHARED(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nh
  do j = 1,Ny
    do k = 1,Nz
      ! Density Equation Evolution.
      rho_k_new(i,j,k) = rho_k(i,j,k) + ( (3.0d0/2.0d0)*d_rho_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_k_dt_old(i,j,k) )*dt

      ! Momentum Equation Evolution.
      rho_ux_k_new(i,j,k) = rho_ux_k(i,j,k) + ( (3.0d0/2.0d0)*d_rho_ux_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_ux_k_dt_old(i,j,k) )*dt
      rho_uy_k_new(i,j,k) = rho_uy_k(i,j,k) + ( (3.0d0/2.0d0)*d_rho_uy_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_uy_k_dt_old(i,j,k) )*dt
      rho_uz_k_new(i,j,k) = rho_uz_k(i,j,k) + ( (3.0d0/2.0d0)*d_rho_uz_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_rho_uz_k_dt_old(i,j,k) )*dt

      ! Energy Equation Evolution.
      E_k_new(i,j,k) = E_k(i,j,k) !+ ( (3.0d0/2.0d0)*d_E_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_E_k_dt_old(i,j,k) )*dt

      ! Energy Equation Evolution.
      Bx_k_new(i,j,k) = Bx_k(i,j,k) + ( (3.0d0/2.0d0)*d_Bx_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_Bx_k_dt_old(i,j,k) )*dt
      By_k_new(i,j,k) = By_k(i,j,k) + ( (3.0d0/2.0d0)*d_By_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_By_k_dt_old(i,j,k) )*dt
      Bz_k_new(i,j,k) = Bz_k(i,j,k) + ( (3.0d0/2.0d0)*d_Bz_k_dt_new(i,j,k) - (1.0d0/2.0d0)*d_Bz_k_dt_old(i,j,k) )*dt
    enddo
  end do
end do

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine ab3com
