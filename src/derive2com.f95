
!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive2com(Nx,Ny,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k, &
                  d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
                  d_E_k_dt_old,d_E_k_dt_new,d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new)
implicit none

include "fftw3.f"

integer ( kind = 4 ) Nx,Ny,Nh, i, j, k
real ( kind = 8 ) pi,time,dt,Lx,Ly,kx,ky,mu,ms,CS,mu_0,eta,spheat,kf,A

real ( kind = 8 ) x(Nx), y(Ny)

real ( kind = 8 ) rho(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny),E(Nx,Ny),P(Nx,Ny),Bx(Nx,Ny),By(Nx,Ny),B2(Nx,Ny),Source(Nx,Ny)
real ( kind = 8 ) rho_ux(Nx,Ny),rho_uy(Nx,Ny),Mom_x_1(Nx,Ny),Mom_x_2(Nx,Ny),Mom_y_1(Nx,Ny),Mom_y_2(Nx,Ny)

real ( kind = 8 ) d_ux_dx(Nx,Ny),d_uy_dy(Nx,Ny),Fx(Nx,Ny),Fy(Nx,Ny),ux_dum(Nx,Ny),uy_dum(Nx,Ny)
real ( kind = 8 ) Energy_x(Nx,Ny),Energy_y(Nx,Ny),E_Visc(Nx,Ny),Mag_x(Nx,Ny),Mag_y(Nx,Ny)
real ( kind = 8 ) d_Bx_dy(Nx,Ny),d_By_dx(Nx,Ny),curl_B(Nx,Ny)

complex ( kind = 8 ) rho_k(Nh,Ny),ux_k(Nh,Ny),uy_k(Nh,Ny),rho_ux_k(Nh,Ny),rho_uy_k(Nh,Ny)
complex ( kind = 8 ) rho_ux_k_dum(Nh,Ny),rho_uy_k_dum(Nh,Ny)
complex ( kind = 8 ) P_k(Nh,Ny),E_k(Nh,Ny),Fx_k(Nh,Ny),Fy_k(Nh,Ny),Bx_k(Nh,Ny),By_k(Nh,Ny)
complex ( kind = 8 ) rho_k_dum(Nh,Ny),ux_k_dum(Nh,Ny),uy_k_dum(Nh,Ny)
complex ( kind = 8 ) E_k_dum(Nh,Ny),Bx_k_dum(Nh,Ny),By_k_dum(Nh,Ny)

complex ( kind = 8 ) i_kx_rho_ux_k(Nh,Ny),i_ky_rho_uy_k(Nh,Ny)
complex ( kind = 8 ) Mom_x_1_k(Nh,Ny),Mom_x_2_k(Nh,Ny),Mom_y_1_k(Nh,Ny),Mom_y_2_k(Nh,Ny)
complex ( kind = 8 ) i_kx_Mom_x_1_k(Nh,Ny),i_ky_Mom_x_2_k(Nh,Ny),i_kx_Mom_y_1_k(Nh,Ny),i_ky_Mom_y_2_k(Nh,Ny)
complex ( kind = 8 ) kx2_ux_k(Nh,Ny),ky2_ux_k(Nh,Ny),kx2_uy_k(Nh,Ny),ky2_uy_k(Nh,Ny)
complex ( kind = 8 ) i_kx_ux_k(Nh,Ny),i_ky_uy_k(Nh,Ny)
complex ( kind = 8 ) Energy_x_k(Nh,Ny),Energy_y_k(Nh,Ny),E_Visc_k(Nh,Ny),Mag_x_k(Nh,Ny),Mag_y_k(Nh,Ny)
complex ( kind = 8 ) i_kx_Energy_x_k(Nh,Ny),i_ky_Energy_y_k(Nh,Ny),i_ky_Mag_x_k(Nh,Ny),i_kx_Mag_y_k(Nh,Ny)
complex ( kind = 8 ) i_ky_Bx_k(Nh,Ny),i_kx_By_k(Nh,Ny)
complex ( kind = 8 ) kx2_Bx_k(Nh,Ny),ky2_Bx_k(Nh,Ny),kx2_By_k(Nh,Ny),ky2_By_k(Nh,Ny)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny),d_rho_ux_k_dt_old(Nh,Ny),d_rho_uy_k_dt_old(Nh,Ny)
complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny),d_Bx_k_dt_old(Nh,Ny),d_By_k_dt_old(Nh,Ny)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny),d_rho_ux_k_dt_new(Nh,Ny),d_rho_uy_k_dt_new(Nh,Ny)
complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny),d_Bx_k_dt_new(Nh,Ny),d_By_k_dt_new(Nh,Ny)

integer ( kind = 8 ) plan_forward,plan_backward ! FFTW

common/comm/Lx,Ly,spheat,mu,ms,CS,mu_0,eta,dt,kf,A

! Keep Backup of the Arrays for FFTW.

!$OMP PARALLEL SHARED(rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k),&
!$OMP & SHARED(rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,E_k_dum,Bx_k_dum,By_k_dum) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
  rho_k_dum(i,j) = rho_k(i,j)
  rho_ux_k_dum(i,j) = rho_ux_k(i,j)
  rho_uy_k_dum(i,j) = rho_uy_k(i,j)
  E_k_dum(i,j) = E_k(i,j)
  Bx_k_dum(i,j) = Bx_k(i,j)
  By_k_dum(i,j) = By_k(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, rho_k_dum, rho, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, rho_ux_k_dum, rho_ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, rho_uy_k_dum, rho_uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, E_k_dum, E, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, Bx_k_dum, Bx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, By_k_dum, By, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

!$OMP PARALLEL SHARED(A,kf,rho,rho_ux,rho_uy,E,Fx,Fy,Bx,By,B2,ux_dum,uy_dum) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
  ! FFTW Normalisation
  rho(i,j) = rho(i,j)/(dfloat(Nx)*dfloat(Ny))

  rho_ux(i,j) = rho_ux(i,j)/(dfloat(Nx)*dfloat(Ny))
  rho_uy(i,j) = rho_uy(i,j)/(dfloat(Nx)*dfloat(Ny))

  E(i,j) = E(i,j)/(dfloat(Nx)*dfloat(Ny))

  Bx(i,j) = Bx(i,j)/(dfloat(Nx)*dfloat(Ny))
  By(i,j) = By(i,j)/(dfloat(Nx)*dfloat(Ny))

  ! Evaluate Velocity in Real Space.
  ux(i,j) = rho_ux(i,j)/rho(i,j)
  uy(i,j) = rho_uy(i,j)/rho(i,j)

  ! Evaluate Forcing.
  Fx(i,j) = - rho(i,j) * A * dsin(kf*y(j))
  Fy(i,j) = + rho(i,j) * A * dsin(kf*x(i))

  ! Evaluate Square of Magnetic Field.
  B2(i,j) = Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j)

  ! Keep Backup of the Arrays for FFTW.
  ux_dum(i,j) = ux(i,j)
  uy_dum(i,j) = uy(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Fx, Fx_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Fy, Fy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

! Evaluate derivatives of Velocity and Magnetic Field.

!$OMP PARALLEL SHARED(Lx,Ly,ux_k,uy_k,Bx_k,By_k),&
!$OMP & SHARED(i_kx_ux_k,i_ky_uy_k,i_ky_Bx_k,i_kx_By_k) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    i_kx_ux_k(i,j) = (0.0d0,1.0d0)*kx*ux_k(i,j)
    i_ky_uy_k(i,j) = (0.0d0,1.0d0)*ky*uy_k(i,j)
    i_ky_Bx_k(i,j) = (0.0d0,1.0d0)*ky*Bx_k(i,j)
    i_kx_By_k(i,j) = (0.0d0,1.0d0)*kx*By_k(i,j)
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    i_kx_ux_k(i,j) = (0.0d0,1.0d0)*kx*ux_k(i,j)
    i_ky_uy_k(i,j) = (0.0d0,1.0d0)*ky*uy_k(i,j)
    i_ky_Bx_k(i,j) = (0.0d0,1.0d0)*ky*Bx_k(i,j)
    i_kx_By_k(i,j) = (0.0d0,1.0d0)*kx*By_k(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_ux_k, d_ux_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_uy_k, d_uy_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_Bx_k, d_Bx_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_By_k, d_By_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

!$OMP PARALLEL SHARED(d_ux_dx,d_uy_dy,d_Bx_dy,d_By_dx,curl_B,P),&
!$OMP & SHARED(Mom_x_1,Mom_x_2,Mom_y_1,Mom_y_2,Energy_x,Energy_y,E_Visc,Mag_x,Mag_y) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
  ! FFTW Normalisation.
  d_ux_dx(i,j) = d_ux_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
  d_uy_dy(i,j) = d_uy_dy(i,j)/(dfloat(Nx)*dfloat(Ny))

  d_Bx_dy(i,j) = d_ux_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
  d_By_dx(i,j) = d_uy_dy(i,j)/(dfloat(Nx)*dfloat(Ny))

  ! Evaluate Curl of Magnetic Field.
  curl_B(i,j) = d_By_dx(i,j) - d_Bx_dy(i,j)

  ! Evaluate Pressure
  P(i,j) = CS*CS*rho(i,j)!( spheat - 1.0d0 ) * ( E(i,j) - 0.50d0 * ( rho_ux(i,j)*ux(i,j)+rho_uy(i,j)*uy(i,j) - B2(i,j) ) )

  ! Evaluate LHS of Momentum Equation.
  Mom_x_1(i,j) = rho_ux(i,j)*ux(i,j) + P(i,j) + B2(i,j)/2.0d0 - Bx(i,j)*Bx(i,j)
  Mom_x_2(i,j) = rho_ux(i,j)*uy(i,j) - Bx(i,j)*By(i,j)

  Mom_y_1(i,j) = rho_ux(i,j)*uy(i,j) - Bx(i,j)*By(i,j)
  Mom_y_2(i,j) = rho_uy(i,j)*uy(i,j) + P(i,j) + B2(i,j)/2.0d0 - By(i,j)*By(i,j)

  ! Evaluate LHS of Energy Equation.
  Energy_x(i,j) = ( E(i,j) + P(i,j) + B2(i,j)/2.0d0 ) * ux(i,j) - ux(i,j)*Bx(i,j)*Bx(i,j) - uy(i,j)*Bx(i,j)*By(i,j)
  Energy_x(i,j) = Energy_x(i,j) - eta * By(i,j) * curl_B(i,j)
  Energy_y(i,j) = ( E(i,j) + P(i,j) + B2(i,j)/2.0d0 ) * uy(i,j) - ux(i,j)*Bx(i,j)*By(i,j) - uy(i,j)*By(i,j)*By(i,j)
  Energy_y(i,j) = Energy_y(i,j) + eta * Bx(i,j) * curl_B(i,j)

  ! Evaluate RHS of Energy Equation.
  E_Visc(i,j) = ( d_ux_dx(i,j) + d_uy_dy(i,j) )**2

  ! Evaluate LHS of Magnetic Field Equation.
  Mag_x(i,j) = uy(i,j)*Bx(i,j) - By(i,j)*ux(i,j)
  Mag_y(i,j) = ux(i,j)*By(i,j) - Bx(i,j)*uy(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Mom_x_1, Mom_x_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Mom_x_2, Mom_x_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Mom_y_1, Mom_y_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Mom_y_2, Mom_y_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Energy_x, Energy_x_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Energy_y, Energy_y_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, E_Visc, E_Visc_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Mag_x, Mag_x_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Mag_y, Mag_y_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_destroy_plan_ (plan_backward)
  call dfftw_destroy_plan_ (plan_forward)

! Evaluate the Derivatives in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,ux_k,uy_k,rho_ux_k,rho_uy_k,Bx_k,By_k),&
!$OMP & SHARED(Mom_x_1_k,Mom_x_2_k,Mom_y_1_k,Mom_y_2_k),&
!$OMP & SHARED(Energy_x_k,Energy_y_k,Mag_x_k,Mag_y_k),&
!$OMP & SHARED(i_kx_rho_ux_k,i_ky_rho_uy_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_ky_Mag_x_k,i_kx_Mag_y_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kx2_uy_k,ky2_uy_k),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kx2_By_k,ky2_By_k) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly

    i_kx_rho_ux_k(i,j) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j)
    i_ky_rho_uy_k(i,j) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j)

    i_kx_Mom_x_1_k(i,j) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j)
    i_ky_Mom_x_2_k(i,j) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j)
    i_kx_Mom_y_1_k(i,j) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j)
    i_ky_Mom_y_2_k(i,j) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j)

    kx2_ux_k(i,j) = kx*kx*ux_k(i,j)
    ky2_ux_k(i,j) = ky*ky*ux_k(i,j)
    kx2_uy_k(i,j) = kx*kx*uy_k(i,j)
    ky2_uy_k(i,j) = ky*ky*uy_k(i,j)

    i_kx_Energy_x_k(i,j) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j)
    i_ky_Energy_y_k(i,j) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j)

    i_ky_Mag_x_k(i,j) = (0.0d0,1.0d0)*ky*Mag_x_k(i,j)
    i_kx_Mag_y_k(i,j) = (0.0d0,1.0d0)*kx*Mag_y_k(i,j)

    kx2_Bx_k(i,j) = kx*kx*Bx_k(i,j)
    ky2_Bx_k(i,j) = ky*ky*Bx_k(i,j)
    kx2_By_k(i,j) = kx*kx*By_k(i,j)
    ky2_By_k(i,j) = ky*ky*By_k(i,j)
  enddo

  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly

    i_kx_rho_ux_k(i,j) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j)
    i_ky_rho_uy_k(i,j) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j)

    i_kx_Mom_x_1_k(i,j) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j)
    i_ky_Mom_x_2_k(i,j) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j)
    i_kx_Mom_y_1_k(i,j) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j)
    i_ky_Mom_y_2_k(i,j) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j)

    kx2_ux_k(i,j) = kx*kx*ux_k(i,j)
    ky2_ux_k(i,j) = ky*ky*ux_k(i,j)
    kx2_uy_k(i,j) = kx*kx*uy_k(i,j)
    ky2_uy_k(i,j) = ky*ky*uy_k(i,j)

    i_kx_Energy_x_k(i,j) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j)
    i_ky_Energy_y_k(i,j) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j)

    i_ky_Mag_x_k(i,j) = (0.0d0,1.0d0)*ky*Mag_x_k(i,j)
    i_kx_Mag_y_k(i,j) = (0.0d0,1.0d0)*kx*Mag_y_k(i,j)

    kx2_Bx_k(i,j) = kx*kx*Bx_k(i,j)
    ky2_Bx_k(i,j) = ky*ky*Bx_k(i,j)
    kx2_By_k(i,j) = kx*kx*By_k(i,j)
    ky2_By_k(i,j) = ky*ky*By_k(i,j)
 enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

! De - Aliazing Technique With 2/3 Rule for All the Non-Linear Terms.

!$OMP PARALLEL SHARED(Lx,Ly,i_kx_rho_ux_k,i_ky_rho_uy_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,E_Visc_k,i_ky_Mag_x_k,i_kx_Mag_y_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kx2_uy_k,ky2_uy_k),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kx2_By_k,ky2_By_k) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 0) then
      i_kx_rho_ux_k(i,j) = 0.0d0
      i_ky_rho_uy_k(i,j) = 0.0d0
      i_kx_Mom_x_1_k(i,j) = 0.0d0
      i_ky_Mom_x_2_k(i,j) = 0.0d0
      i_kx_Mom_y_1_k(i,j) = 0.0d0
      i_ky_Mom_y_2_k(i,j) = 0.0d0
      kx2_ux_k(i,j) = 0.0d0
      ky2_ux_k(i,j) = 0.0d0
      kx2_uy_k(i,j) = 0.0d0
      ky2_uy_k(i,j) = 0.0d0
      i_kx_Energy_x_k(i,j) = 0.0d0
      i_ky_Energy_y_k(i,j) = 0.0d0
      E_Visc_k(i,j) = 0.0d0
      i_ky_Mag_x_k(i,j) = 0.0d0
      i_kx_Mag_y_k(i,j) = 0.0d0
      kx2_Bx_k(i,j) = 0.0d0
      ky2_Bx_k(i,j) = 0.0d0
      kx2_By_k(i,j) = 0.0d0
      ky2_By_k(i,j) = 0.0d0
      endif
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 0) then
      i_kx_rho_ux_k(i,j) = 0.0d0
      i_ky_rho_uy_k(i,j) = 0.0d0
      i_kx_Mom_x_1_k(i,j) = 0.0d0
      i_ky_Mom_x_2_k(i,j) = 0.0d0
      i_kx_Mom_y_1_k(i,j) = 0.0d0
      i_ky_Mom_y_2_k(i,j) = 0.0d0
      kx2_ux_k(i,j) = 0.0d0
      ky2_ux_k(i,j) = 0.0d0
      kx2_uy_k(i,j) = 0.0d0
      ky2_uy_k(i,j) = 0.0d0
      i_kx_Energy_x_k(i,j) = 0.0d0
      i_ky_Energy_y_k(i,j) = 0.0d0
      E_Visc_k(i,j) = 0.0d0
      i_ky_Mag_x_k(i,j) = 0.0d0
      i_kx_Mag_y_k(i,j) = 0.0d0
      kx2_Bx_k(i,j) = 0.0d0
      ky2_Bx_k(i,j) = 0.0d0
      kx2_By_k(i,j) = 0.0d0
      ky2_By_k(i,j) = 0.0d0
      endif
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(mu,eta,i_kx_rho_ux_k,i_ky_rho_uy_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,E_Visc_k,i_ky_Mag_x_k,i_kx_Mag_y_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kx2_uy_k,ky2_uy_k,Fx_k,Fy_k),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kx2_By_k,ky2_By_k),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new),&
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
    ! Density Equation.
    d_rho_k_dt_new(i,j) = - ( i_kx_rho_ux_k(i,j) + i_ky_rho_uy_k(i,j) )

    ! Momentum Equation.
    d_rho_ux_k_dt_new(i,j) = - ( i_kx_Mom_x_1_k(i,j) + i_ky_Mom_x_2_k(i,j) )
    d_rho_ux_k_dt_new(i,j) = d_rho_ux_k_dt_new(i,j) - mu * ( kx2_ux_k(i,j) + ky2_ux_k(i,j) ) + Fx_k(i,j)

    d_rho_uy_k_dt_new(i,j) = - ( i_kx_Mom_y_1_k(i,j) + i_ky_Mom_y_2_k(i,j) )
    d_rho_uy_k_dt_new(i,j) = d_rho_uy_k_dt_new(i,j) - mu * ( kx2_uy_k(i,j) + ky2_uy_k(i,j) ) + Fy_k(i,j)

    ! Energy Equation.
    d_E_k_dt_new(i,j) = - ( i_kx_Energy_x_k(i,j) + i_ky_Energy_y_k(i,j) )
    d_E_k_dt_new(i,j) = d_E_k_dt_new(i,j) + mu * E_Visc_k(i,j)

    ! Magnetic Field Equation.
    d_Bx_k_dt_new(i,j) = - i_ky_Mag_x_k(i,j) - eta * ( kx2_Bx_k(i,j) + ky2_Bx_k(i,j) )
    d_By_k_dt_new(i,j) = - i_kx_Mag_y_k(i,j) - eta * ( kx2_By_k(i,j) + ky2_By_k(i,j) )
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine derive2com
