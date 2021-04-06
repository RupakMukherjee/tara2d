
!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive3com(Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
                  d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
                  d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
                  d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
implicit none
include "fftw3.f"
integer ( kind = 4 ) Nx,Ny,Nz,Nh,iret,thread_num, i, j, k
real ( kind = 8 ) pi,time,dt,Lx,Ly,Lz,dx,dy,dz,kx,ky,kz,mu,ms,CS,mu_0,eta,spheat,kf,A,B,C,time_max

real ( kind = 8 ) x(Nx), y(Ny), z(Nz)

real ( kind = 8 ) rho(Nx,Ny,Nz),rho_dum(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz),E(Nx,Ny,Nz),P(Nx,Ny,Nz)
real ( kind = 8 ) Bx(Nx,Ny,Nz),By(Nx,Ny,Nz),Bz(Nx,Ny,Nz),B2(Nx,Ny,Nz),Source(Nx,Ny,Nz)
real ( kind = 8 ) rho_ux(Nx,Ny,Nz),rho_uy(Nx,Ny,Nz),rho_uz(Nx,Ny,Nz)
real ( kind = 8 ) ux_dum(Nx,Ny,Nz),uy_dum(Nx,Ny,Nz),uz_dum(Nx,Ny,Nz)
real ( kind = 8 ) Mom_x_1(Nx,Ny,Nz),Mom_x_2(Nx,Ny,Nz),Mom_x_3(Nx,Ny,Nz),Mom_y_1(Nx,Ny,Nz),Mom_y_2(Nx,Ny,Nz),Mom_y_3(Nx,Ny,Nz)
real ( kind = 8 ) Mom_z_1(Nx,Ny,Nz),Mom_z_2(Nx,Ny,Nz),Mom_z_3(Nx,Ny,Nz)

real ( kind = 8 ) d_ux_dx(Nx,Ny,Nz),d_uy_dy(Nx,Ny,Nz),d_uz_dz(Nx,Ny,Nz),Fx(Nx,Ny,Nz),Fy(Nx,Ny,Nz),Fz(Nx,Ny,Nz)
real ( kind = 8 ) Energy_x(Nx,Ny,Nz),Energy_y(Nx,Ny,Nz),Energy_z(Nx,Ny,Nz),E_Visc(Nx,Ny,Nz)
real ( kind = 8 ) Mag_x_1(Nx,Ny,Nz),Mag_x_2(Nx,Ny,Nz),Mag_y_1(Nx,Ny,Nz),Mag_y_2(Nx,Ny,Nz),Mag_z_1(Nx,Ny,Nz),Mag_z_2(Nx,Ny,Nz)
real ( kind = 8 ) d_Bx_dy(Nx,Ny,Nz),d_By_dx(Nx,Ny,Nz),d_Bx_dz(Nx,Ny,Nz)
real ( kind = 8 ) d_By_dz(Nx,Ny,Nz),d_Bz_dx(Nx,Ny,Nz),d_Bz_dy(Nx,Ny,Nz)
real ( kind = 8 ) curl_x_B(Nx,Ny,Nz),curl_y_B(Nx,Ny,Nz),curl_z_B(Nx,Ny,Nz)

complex ( kind = 8 ) rho_k(Nh,Ny,Nz),ux_k(Nh,Ny,Nz),uy_k(Nh,Ny,Nz),uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) nu(Nh,Ny,Nz),Fx_k(Nh,Ny,Nz),Fy_k(Nh,Ny,Nz),Fz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k(Nh,Ny,Nz),rho_uy_k(Nh,Ny,Nz),rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_ux_k_dum(Nh,Ny,Nz),rho_uy_k_dum(Nh,Ny,Nz),rho_uz_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) P_k(Nh,Ny,Nz),E_k(Nh,Ny,Nz),Bx_k(Nh,Ny,Nz),By_k(Nh,Ny,Nz),Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) rho_k_dum(Nh,Ny,Nz),ux_k_dum(Nh,Ny,Nz),uy_k_dum(Nh,Ny,Nz)
complex ( kind = 8 ) E_k_dum(Nh,Ny,Nz),Bx_k_dum(Nh,Ny,Nz),By_k_dum(Nh,Ny,Nz),Bz_k_dum(Nh,Ny,Nz)

complex ( kind = 8 ) i_kx_rho_ux_k(Nh,Ny,Nz),i_ky_rho_uy_k(Nh,Ny,Nz),i_kz_rho_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_x_1_k(Nh,Ny,Nz),Mom_x_2_k(Nh,Ny,Nz),Mom_x_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_y_1_k(Nh,Ny,Nz),Mom_y_2_k(Nh,Ny,Nz),Mom_y_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mom_z_1_k(Nh,Ny,Nz),Mom_z_2_k(Nh,Ny,Nz),Mom_z_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_x_1_k(Nh,Ny,Nz),i_ky_Mom_x_2_k(Nh,Ny,Nz),i_kz_Mom_x_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_y_1_k(Nh,Ny,Nz),i_ky_Mom_y_2_k(Nh,Ny,Nz),i_kz_Mom_y_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Mom_z_1_k(Nh,Ny,Nz),i_ky_Mom_z_2_k(Nh,Ny,Nz),i_kz_Mom_z_3_k(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_ux_k(Nh,Ny,Nz),ky2_ux_k(Nh,Ny,Nz),kz2_ux_k(Nh,Ny,Nz),kx2_uy_k(Nh,Ny,Nz),ky2_uy_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_uy_k(Nh,Ny,Nz),kx2_uz_k(Nh,Ny,Nz),ky2_uz_k(Nh,Ny,Nz),kz2_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_ux_k(Nh,Ny,Nz),i_ky_uy_k(Nh,Ny,Nz),i_kz_uz_k(Nh,Ny,Nz)
complex ( kind = 8 ) Energy_x_k(Nh,Ny,Nz),Energy_y_k(Nh,Ny,Nz),Energy_z_k(Nh,Ny,Nz),E_Visc_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mag_x_1_k(Nh,Ny,Nz),Mag_x_2_k(Nh,Ny,Nz),Mag_y_1_k(Nh,Ny,Nz),Mag_y_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) Mag_z_1_k(Nh,Ny,Nz),Mag_z_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kx_Energy_x_k(Nh,Ny,Nz),i_ky_Energy_y_k(Nh,Ny,Nz),i_kz_Energy_z_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Mag_x_1_k(Nh,Ny,Nz),i_kz_Mag_x_2_k(Nh,Ny,Nz),i_kx_Mag_y_1_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Mag_y_2_k(Nh,Ny,Nz),i_kx_Mag_z_1_k(Nh,Ny,Nz),i_ky_Mag_z_2_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_ky_Bx_k(Nh,Ny,Nz),i_kx_By_k(Nh,Ny,Nz),i_kx_Bz_k(Nh,Ny,Nz)
complex ( kind = 8 ) i_kz_Bx_k(Nh,Ny,Nz),i_ky_Bz_k(Nh,Ny,Nz),i_kz_By_k(Nh,Ny,Nz)
complex ( kind = 8 ) kx2_Bx_k(Nh,Ny,Nz),ky2_Bx_k(Nh,Ny,Nz),kz2_Bx_k(Nh,Ny,Nz),kx2_By_k(Nh,Ny,Nz),ky2_By_k(Nh,Ny,Nz)
complex ( kind = 8 ) kz2_By_k(Nh,Ny,Nz),kx2_Bz_k(Nh,Ny,Nz),ky2_Bz_k(Nh,Ny,Nz),kz2_Bz_k(Nh,Ny,Nz)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny,Nz),d_rho_ux_k_dt_old(Nh,Ny,Nz),d_rho_uy_k_dt_old(Nh,Ny,Nz),d_rho_uz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_old(Nh,Ny,Nz),d_Bx_k_dt_old(Nh,Ny,Nz),d_By_k_dt_old(Nh,Ny,Nz),d_Bz_k_dt_old(Nh,Ny,Nz)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny,Nz),d_rho_ux_k_dt_new(Nh,Ny,Nz),d_rho_uy_k_dt_new(Nh,Ny,Nz),d_rho_uz_k_dt_new(Nh,Ny,Nz)
complex ( kind = 8 ) d_E_k_dt_new(Nh,Ny,Nz),d_Bx_k_dt_new(Nh,Ny,Nz),d_By_k_dt_new(Nh,Ny,Nz),d_Bz_k_dt_new(Nh,Ny,Nz)

integer ( kind = 8 ) plan_forward,plan_backward ! FFTW

common/comm/iret,thread_num,time_max,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,dt,kf

! Keep Backup of the Arrays for FFTW.

!$OMP PARALLEL SHARED(rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum),&
!$OMP & SHARED(E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
    do k = 1,Nz
    rho_k_dum(i,j,k) = rho_k(i,j,k)
    rho_ux_k_dum(i,j,k) = rho_ux_k(i,j,k)
    rho_uy_k_dum(i,j,k) = rho_uy_k(i,j,k)
    rho_uz_k_dum(i,j,k) = rho_uz_k(i,j,k)
    E_k_dum(i,j,k) = E_k(i,j,k)
    Bx_k_dum(i,j,k) = Bx_k(i,j,k)
    By_k_dum(i,j,k) = By_k(i,j,k)
    Bz_k_dum(i,j,k) = Bz_k(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_k_dum, rho, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_ux_k_dum, rho_ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uy_k_dum, rho_uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, rho_uz_k_dum, rho_uz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, E_k_dum, E, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bx_k_dum, Bx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, By_k_dum, By, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Bz_k_dum, Bz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

!if (t <= 100  ) then
 A = 0.10d0
 B = 0.10d0
 C = 0.10d0
!else
! A = 0.0d0
! B = 0.0d0
! C = 0.0d0
!endif

!$OMP PARALLEL SHARED(dx,dy,dz,A,B,C,kf,x,y,z,rho,ux,uy,uz,rho_dum,ux_dum,uy_dum,uz_dum),&
!$OMP & SHARED(rho_ux,rho_uy,rho_uz,Fx,Fy,Fz,E,Bx,By,Bz,B2) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx
x(i)=0.0d0+real(i-1)*dx
  do j = 1,Ny
  y(j)=0.0d0+real(j-1)*dy
    do k = 1,Nz
    z(k)=0.0d0+real(k-1)*dz
    ! FFTW Normalisation
    rho(i,j,k) = rho(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

    rho_ux(i,j,k) = rho_ux(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    rho_uy(i,j,k) = rho_uy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    rho_uz(i,j,k) = rho_uz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

    E(i,j,k) = E(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

    Bx(i,j,k) = Bx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    By(i,j,k) = By(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    Bz(i,j,k) = Bz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

    ! Evaluate Velocity in Real Space.
    ux(i,j,k) = rho_ux(i,j,k)/rho(i,j,k)
    uy(i,j,k) = rho_uy(i,j,k)/rho(i,j,k)
    uz(i,j,k) = rho_uz(i,j,k)/rho(i,j,k)

    ! Evaluate Forcing.
    Fx(i,j,k) = rho(i,j,k) * ( A*dsin(kf*z(k)) + C*dcos(kf*y(j)) )
    Fy(i,j,k) = rho(i,j,k) * ( B*dsin(kf*x(i)) + A*dcos(kf*z(k)) )
    Fz(i,j,k) = rho(i,j,k) * ( C*dsin(kf*y(j)) + B*dcos(kf*x(i)) )

    ! Evaluate Square of Magnetic Field.
    B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)

    ! Keep Backup of the Arrays for FFTW.
    ux_dum(i,j,k) = ux(i,j,k)
    uy_dum(i,j,k) = uy(i,j,k)
    uz_dum(i,j,k) = uz(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, uz_dum, uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fx, Fx_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fy, Fy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Fz, Fz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

! Evaluate derivatives of Velocity and Magnetic Field.

!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,Bx_k,By_k),&
!$OMP & SHARED(i_kx_ux_k,i_ky_uy_k,i_kz_uz_k),&
!$OMP & SHARED(i_ky_Bx_k,i_kz_Bx_k,i_kx_By_k,i_kz_By_k,i_kx_Bz_k,i_ky_Bz_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*float(k-1)/Lz
      i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
      i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)

      i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
      i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)

      i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k-1)/Lz
      i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
      i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)

      i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
      i_kx_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*ux_k(i,j,k)
      i_ky_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*uy_k(i,j,k)
      i_kz_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*uz_k(i,j,k)

      i_ky_Bx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      i_kz_Bx_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k)
      i_kx_By_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k)
      i_kz_By_k(i,j,k) = (0.0d0,1.0d0)*kz*By_k(i,j,k)
      i_kx_Bz_k(i,j,k) = (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      i_ky_Bz_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_ux_k, d_ux_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_uy_k, d_uy_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_uz_k, d_uz_dz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_Bx_k, d_Bx_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_Bx_k, d_Bx_dz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_By_k, d_By_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kz_By_k, d_By_dz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_kx_Bz_k, d_Bz_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, i_ky_Bz_k, d_Bz_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

!$OMP PARALLEL SHARED(spheat,CS,eta,ux,uy,uz,rho,rho_ux,rho_uy,rho_uz,Bx,By,Bz),&
!$OMP & SHARED(d_ux_dx,d_uy_dy,d_uz_dz,P,E),&
!$OMP & SHARED(d_Bx_dy,d_Bx_dz,d_By_dx,d_By_dz,d_Bz_dx,d_Bz_dy,curl_x_B,curl_y_B,curl_z_B,B2),&
!$OMP & SHARED(Mom_x_1,Mom_x_2,Mom_x_3,Mom_y_1,Mom_y_2,Mom_y_3,Mom_z_1,Mom_z_2,Mom_z_3),&
!$OMP & SHARED(Energy_x,Energy_y,Energy_z,E_Visc),&
!$OMP & SHARED(Mag_x_1,Mag_x_2,Mag_y_1,Mag_y_2,Mag_z_1,Mag_z_2) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
    ! FFTW Normalisation.
    d_ux_dx(i,j,k) = d_ux_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_uy_dy(i,j,k) = d_uy_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_uz_dz(i,j,k) = d_uz_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

    d_Bx_dy(i,j,k) = d_Bx_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bx_dz(i,j,k) = d_Bx_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_By_dx(i,j,k) = d_By_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_By_dz(i,j,k) = d_By_dz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bz_dx(i,j,k) = d_Bz_dx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
    d_Bz_dy(i,j,k) = d_Bz_dy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

    ! Evaluate Curl of Magnetic Field.
    curl_x_B(i,j,k) = d_Bz_dy(i,j,k) - d_By_dz(i,j,k)
    curl_y_B(i,j,k) = d_Bz_dx(i,j,k) - d_Bx_dz(i,j,k)
    curl_z_B(i,j,k) = d_By_dx(i,j,k) - d_Bx_dy(i,j,k)

    ! Evaluate Pressure
    P(i,j,k) = CS*CS*rho(i,j,k)!( spheat - 1.0d0 ) * ( E(i,j,k) &
               !- 0.50d0 * ( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) &
               !- B2(i,j,k) ) )

    ! Evaluate LHS of Momentum Equation.
    Mom_x_1(i,j,k) = rho_ux(i,j,k)*ux(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - Bx(i,j,k)*Bx(i,j,k)
    Mom_x_2(i,j,k) = rho_ux(i,j,k)*uy(i,j,k) - Bx(i,j,k)*By(i,j,k)
    Mom_x_3(i,j,k) = rho_ux(i,j,k)*uz(i,j,k) - Bx(i,j,k)*Bz(i,j,k)

    Mom_y_1(i,j,k) = rho_ux(i,j,k)*uy(i,j,k) - Bx(i,j,k)*By(i,j,k)
    Mom_y_2(i,j,k) = rho_uy(i,j,k)*uy(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - By(i,j,k)*By(i,j,k)
    Mom_y_3(i,j,k) = rho_uy(i,j,k)*uz(i,j,k) - By(i,j,k)*Bz(i,j,k)

    Mom_z_1(i,j,k) = rho_uz(i,j,k)*ux(i,j,k) - Bz(i,j,k)*Bx(i,j,k)
    Mom_z_2(i,j,k) = rho_uz(i,j,k)*uy(i,j,k) - Bz(i,j,k)*By(i,j,k)
    Mom_z_3(i,j,k) = rho_uz(i,j,k)*uz(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 - Bz(i,j,k)*Bz(i,j,k)

    ! Evaluate LHS of Energy Equation.
    Energy_x(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * ux(i,j,k) &
                      - ux(i,j,k)*Bx(i,j,k)*Bx(i,j,k) - uy(i,j,k)*Bx(i,j,k)*By(i,j,k) - uz(i,j,k)*Bx(i,j,k)*Bz(i,j,k) &
                      - eta * ( By(i,j,k) * curl_z_B(i,j,k) + Bz(i,j,k) * curl_y_B(i,j,k))

    Energy_y(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * uy(i,j,k) &
                      - ux(i,j,k)*Bx(i,j,k)*By(i,j,k) - uy(i,j,k)*By(i,j,k)*By(i,j,k) - uz(i,j,k)*By(i,j,k)*Bz(i,j,k) &
                      + eta * ( Bx(i,j,k) * curl_z_B(i,j,k) - Bz(i,j,k) * curl_x_B(i,j,k))

    Energy_z(i,j,k) = ( E(i,j,k) + P(i,j,k) + B2(i,j,k)/2.0d0 ) * uz(i,j,k) &
                      - ux(i,j,k)*Bz(i,j,k)*Bx(i,j,k) - uy(i,j,k)*Bz(i,j,k)*By(i,j,k) - uz(i,j,k)*Bz(i,j,k)*Bz(i,j,k) &
                      + eta * ( Bx(i,j,k) * curl_y_B(i,j,k) + By(i,j,k) * curl_x_B(i,j,k))

    ! Evaluate RHS of Energy Equation.
    E_Visc(i,j,k) = ( d_ux_dx(i,j,k) + d_uy_dy(i,j,k) + d_uz_dz(i,j,k) )**2

    ! Evaluate LHS of Magnetic Field Equation.
    Mag_x_1(i,j,k) = ux(i,j,k)*By(i,j,k) - uy(i,j,k)*Bx(i,j,k)
    Mag_x_2(i,j,k) = ux(i,j,k)*Bz(i,j,k) - uz(i,j,k)*Bx(i,j,k)

    Mag_y_1(i,j,k) = ux(i,j,k)*By(i,j,k) - uy(i,j,k)*Bx(i,j,k)
    Mag_y_2(i,j,k) = uy(i,j,k)*Bz(i,j,k) - uz(i,j,k)*By(i,j,k)

    Mag_z_1(i,j,k) = ux(i,j,k)*Bz(i,j,k) - uz(i,j,k)*Bx(i,j,k)
    Mag_z_2(i,j,k) = uy(i,j,k)*Bz(i,j,k) - uz(i,j,k)*By(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_1, Mom_x_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_2, Mom_x_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_x_3, Mom_x_3_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_1, Mom_y_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_2, Mom_y_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_y_3, Mom_y_3_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_1, Mom_z_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_2, Mom_z_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mom_z_3, Mom_z_3_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_x, Energy_x_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_y, Energy_y_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Energy_z, Energy_z_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E_Visc, E_Visc_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_x_1, Mag_x_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_x_2, Mag_x_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_y_1, Mag_y_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_y_2, Mag_y_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_z_1, Mag_z_1_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Mag_z_2, Mag_z_2_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)


! Evaluate the Derivatives in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,rho_ux_k,rho_uy_k,rho_uz_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(Mom_x_1_k,Mom_x_2_k,Mom_x_3_k,Mom_y_1_k,Mom_y_2_k,Mom_y_3_k),&
!$OMP & SHARED(Mom_z_1_k,Mom_z_2_k,Mom_z_3_k),&
!$OMP & SHARED(Energy_x_k,Energy_y_k,Energy_z_k,Mag_x_1_k,Mag_x_2_k,Mag_y_1_k,Mag_y_2_k,Mag_z_1_k,Mag_z_2_k),&
!$OMP & SHARED(i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k),&
!$OMP & SHARED(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k),&
!$OMP & SHARED(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*float(k-1)/Lz

      i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
      i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k)
      i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k)

      i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)
      i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)
      i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k)
      i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)
      i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)
      i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)

      kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
      ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
      kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
      kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
      ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k)
      kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
      kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
      ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
      kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

      i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)

      i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)

      kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
      kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
      ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k)
      kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
      kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float(j-1)/Ly
    kz = 2.0d0*pi*float((k-1)-Nz)/Lz

      i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
      i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k)
      i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k)

      i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)
      i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)
      i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k)
      i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)
      i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)
      i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)

      kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
      ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
      kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
      kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
      ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k)
      kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
      kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
      ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
      kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

      i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)

      i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)

      kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
      kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
      ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k)
      kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
      kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
    enddo
  enddo

  do j = Ny/2+1,Ny
    do k = 1,Nz/2
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
    kz = 2.0d0*pi*float(k-1)/Lz

      i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
      i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k)
      i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k)

      i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)
      i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)
      i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k)
      i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)
      i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)
      i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)

      kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
      ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
      kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
      kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
      ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k)
      kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
      kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
      ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
      kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

      i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)

      i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)

      kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
      kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
      ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k)
      kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
      kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
    kx = 2.0d0*pi*float(i-1)/Lx
    ky = 2.0d0*pi*float((j-1)-Ny)/Ly
    kz = 2.0d0*pi*float((k-1)-Nz)/Lz

      i_kx_rho_ux_k(i,j,k) = (0.0d0,1.0d0)*kx*rho_ux_k(i,j,k)
      i_ky_rho_uy_k(i,j,k) = (0.0d0,1.0d0)*ky*rho_uy_k(i,j,k)
      i_kz_rho_uz_k(i,j,k) = (0.0d0,1.0d0)*kz*rho_uz_k(i,j,k)

      i_kx_Mom_x_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_x_1_k(i,j,k)
      i_ky_Mom_x_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_x_2_k(i,j,k)
      i_kz_Mom_x_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_x_3_k(i,j,k)
      i_kx_Mom_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_y_1_k(i,j,k)
      i_ky_Mom_y_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_y_2_k(i,j,k)
      i_kz_Mom_y_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_y_3_k(i,j,k)
      i_kx_Mom_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mom_z_1_k(i,j,k)
      i_ky_Mom_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mom_z_2_k(i,j,k)
      i_kz_Mom_z_3_k(i,j,k) = (0.0d0,1.0d0)*kz*Mom_z_3_k(i,j,k)

      kx2_ux_k(i,j,k) = kx*kx*ux_k(i,j,k)
      ky2_ux_k(i,j,k) = ky*ky*ux_k(i,j,k)
      kz2_ux_k(i,j,k) = kz*kz*ux_k(i,j,k)
      kx2_uy_k(i,j,k) = kx*kx*uy_k(i,j,k)
      ky2_uy_k(i,j,k) = ky*ky*uy_k(i,j,k)
      kz2_uy_k(i,j,k) = kz*kz*uy_k(i,j,k)
      kx2_uz_k(i,j,k) = kx*kx*uz_k(i,j,k)
      ky2_uz_k(i,j,k) = ky*ky*uz_k(i,j,k)
      kz2_uz_k(i,j,k) = kz*kz*uz_k(i,j,k)

      i_kx_Energy_x_k(i,j,k) = (0.0d0,1.0d0)*kx*Energy_x_k(i,j,k)
      i_ky_Energy_y_k(i,j,k) = (0.0d0,1.0d0)*ky*Energy_y_k(i,j,k)
      i_kz_Energy_z_k(i,j,k) = (0.0d0,1.0d0)*kz*Energy_z_k(i,j,k)

      i_ky_Mag_x_1_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_x_1_k(i,j,k)
      i_kz_Mag_x_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_x_2_k(i,j,k)
      i_kx_Mag_y_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_y_1_k(i,j,k)
      i_kz_Mag_y_2_k(i,j,k) = (0.0d0,1.0d0)*kz*Mag_y_2_k(i,j,k)
      i_kx_Mag_z_1_k(i,j,k) = (0.0d0,1.0d0)*kx*Mag_z_1_k(i,j,k)
      i_ky_Mag_z_2_k(i,j,k) = (0.0d0,1.0d0)*ky*Mag_z_2_k(i,j,k)

      kx2_Bx_k(i,j,k) = kx*kx*Bx_k(i,j,k)
      ky2_Bx_k(i,j,k) = ky*ky*Bx_k(i,j,k)
      kz2_Bx_k(i,j,k) = kz*kz*Bx_k(i,j,k)
      kx2_By_k(i,j,k) = kx*kx*By_k(i,j,k)
      ky2_By_k(i,j,k) = ky*ky*By_k(i,j,k)
      kz2_By_k(i,j,k) = kz*kz*By_k(i,j,k)
      kx2_Bz_k(i,j,k) = kx*kx*Bz_k(i,j,k)
      ky2_Bz_k(i,j,k) = ky*ky*Bz_k(i,j,k)
      kz2_Bz_k(i,j,k) = kz*kz*Bz_k(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

! De - Aliazing Technique With 2/3 Rule for All the Non-Linear Terms.

!$OMP PARALLEL SHARED(Lx,Ly,Lz,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k),&
!$OMP & SHARED(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k)&
!$OMP & SHARED(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,nu),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
        if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float(j-1)/Ly
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
        if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float(k-1)/Lz
        if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*float(i-1)/Lx
      ky = 2.0d0*pi*float((j-1)-Ny)/Ly
      kz = 2.0d0*pi*float((k-1)-Nz)/Lz
        if (dsqrt(kx*kx + ky*ky + kz*kz) .ge. (dfloat(Nx+Ny+Nz)/3.0)/3.0 + 0) then
        i_kx_rho_ux_k(i,j,k) = 0.0d0
        i_ky_rho_uy_k(i,j,k) = 0.0d0
        i_kz_rho_uz_k(i,j,k) = 0.0d0
        i_kx_Mom_x_1_k(i,j,k) = 0.0d0
        i_ky_Mom_x_2_k(i,j,k) = 0.0d0
        i_kz_Mom_x_3_k(i,j,k) = 0.0d0
        i_kx_Mom_y_1_k(i,j,k) = 0.0d0
        i_ky_Mom_y_2_k(i,j,k) = 0.0d0
        i_kz_Mom_y_3_k(i,j,k) = 0.0d0
        i_kx_Mom_z_1_k(i,j,k) = 0.0d0
        i_ky_Mom_z_2_k(i,j,k) = 0.0d0
        i_kz_Mom_z_3_k(i,j,k) = 0.0d0
        kx2_ux_k(i,j,k) = 0.0d0
        ky2_ux_k(i,j,k) = 0.0d0
        kz2_ux_k(i,j,k) = 0.0d0
        kx2_uy_k(i,j,k) = 0.0d0
        ky2_uy_k(i,j,k) = 0.0d0
        kz2_uy_k(i,j,k) = 0.0d0
        kx2_uz_k(i,j,k) = 0.0d0
        ky2_uz_k(i,j,k) = 0.0d0
        kz2_uz_k(i,j,k) = 0.0d0
        i_kx_Energy_x_k(i,j,k) = 0.0d0
        i_ky_Energy_y_k(i,j,k) = 0.0d0
        i_kz_Energy_z_k(i,j,k) = 0.0d0
        E_Visc_k(i,j,k) = 0.0d0
        i_ky_Mag_x_1_k(i,j,k) = 0.0d0
        i_kz_Mag_x_2_k(i,j,k) = 0.0d0
        i_kx_Mag_y_1_k(i,j,k) = 0.0d0
        i_kz_Mag_y_2_k(i,j,k) = 0.0d0
        i_kx_Mag_z_1_k(i,j,k) = 0.0d0
        i_ky_Mag_z_2_k(i,j,k) = 0.0d0
        kx2_Bx_k(i,j,k) = 0.0d0
        ky2_Bx_k(i,j,k) = 0.0d0
        kz2_Bx_k(i,j,k) = 0.0d0
        kx2_By_k(i,j,k) = 0.0d0
        ky2_By_k(i,j,k) = 0.0d0
        kz2_By_k(i,j,k) = 0.0d0
        kx2_Bz_k(i,j,k) = 0.0d0
        ky2_Bz_k(i,j,k) = 0.0d0
        kz2_Bz_k(i,j,k) = 0.0d0
        nu(i,j,k) = 0.0d0
        endif
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(mu,eta,i_kx_rho_ux_k,i_ky_rho_uy_k,i_kz_rho_uz_k),&
!$OMP & SHARED(i_kx_Mom_x_1_k,i_ky_Mom_x_2_k,i_kz_Mom_x_3_k,i_kx_Mom_y_1_k,i_ky_Mom_y_2_k,i_kz_Mom_y_3_k),&
!$OMP & SHARED(i_kx_Mom_z_1_k,i_ky_Mom_z_2_k,i_kz_Mom_z_3_k),&
!$OMP & SHARED(i_kx_Energy_x_k,i_ky_Energy_y_k,i_kz_Energy_z_k,E_Visc_k)&
!$OMP & SHARED(i_ky_Mag_x_1_k,i_kz_Mag_x_2_k,i_kx_Mag_y_1_k,i_kz_Mag_y_2_k,i_kx_Mag_z_1_k,i_ky_Mag_z_2_k),&
!$OMP & SHARED(kx2_ux_k,ky2_ux_k,kz2_ux_k,kx2_uy_k,ky2_uy_k,kz2_uy_k,kx2_uz_k,ky2_uz_k,kz2_uz_k,nu,Fx_k,Fy_k,Fz_k),&
!$OMP & SHARED(kx2_Bx_k,ky2_Bx_k,kz2_Bx_k,kx2_By_k,ky2_By_k,kz2_By_k,kx2_Bz_k,ky2_Bz_k,kz2_Bz_k),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),&
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
    do k = 1,Nz
      ! Density Equation.
      d_rho_k_dt_new(i,j,k) = - ( i_kx_rho_ux_k(i,j,k) + i_ky_rho_uy_k(i,j,k) + i_kz_rho_uz_k(i,j,k) )

      ! Momentum Equation.
      d_rho_ux_k_dt_new(i,j,k) = - ( i_kx_Mom_x_1_k(i,j,k) + i_ky_Mom_x_2_k(i,j,k) + i_kz_Mom_x_3_k(i,j,k) ) &
                                 - ( kx2_ux_k(i,j,k) + ky2_ux_k(i,j,k) + kz2_ux_k(i,j,k) ) / 450.0d0 !+ Fx_k(i,j,k)

      d_rho_uy_k_dt_new(i,j,k) = - ( i_kx_Mom_y_1_k(i,j,k) + i_ky_Mom_y_2_k(i,j,k) + i_kz_Mom_y_3_k(i,j,k) ) &
                                 - ( kx2_uy_k(i,j,k) + ky2_uy_k(i,j,k) + kz2_uy_k(i,j,k) ) / 450.0d0 !+ Fy_k(i,j,k)

      d_rho_uz_k_dt_new(i,j,k) = - ( i_kx_Mom_z_1_k(i,j,k) + i_ky_Mom_z_2_k(i,j,k) + i_kz_Mom_z_3_k(i,j,k) ) &
                                 - ( kx2_uz_k(i,j,k) + ky2_uz_k(i,j,k) + kz2_uz_k(i,j,k) ) / 450.0d0 !+ Fz_k(i,j,k)

      ! Energy Equation.
      d_E_k_dt_new(i,j,k) = - ( i_kx_Energy_x_k(i,j,k) + i_ky_Energy_y_k(i,j,k) + i_kz_Energy_z_k(i,j,k) ) &
                            + mu * E_Visc_k(i,j,k)

      ! Magnetic Field Equation.
      d_Bx_k_dt_new(i,j,k) = + ( i_ky_Mag_x_1_k(i,j,k) + i_kz_Mag_x_2_k(i,j,k) ) &
                             - ( kx2_Bx_k(i,j,k) + ky2_Bx_k(i,j,k) + kz2_Bx_k(i,j,k) ) / 450.0d0

      d_By_k_dt_new(i,j,k) = - ( i_kx_Mag_y_1_k(i,j,k) - i_kz_Mag_y_2_k(i,j,k) ) &
                             - ( kx2_By_k(i,j,k) + ky2_By_k(i,j,k) + kz2_By_k(i,j,k) ) / 450.0d0

      d_Bz_k_dt_new(i,j,k) = - ( i_kx_Mag_z_1_k(i,j,k) + i_ky_Mag_z_2_k(i,j,k) ) &
                             - ( kx2_Bz_k(i,j,k) + ky2_Bz_k(i,j,k) + kz2_Bz_k(i,j,k) ) / 450.0d0
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine derive3com
