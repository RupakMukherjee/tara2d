
!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive2coh(Nx,Ny,Nh,pi,time,nk,ukx,uky, &
dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new)
implicit none

include "fftw3.f"

integer ( kind = 4 ) Nx,Ny,Nh, i, j
real ( kind = 8 ) pi,time,dt,Lx,Ly,mux,muy,ms,CS,kx,ky
real ( kind = 8 ) n(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny)
real ( kind = 8 ) Px(Nx,Ny),Py(Nx,Ny),Viscxx(Nx,Ny),Viscxy(Nx,Ny),Viscyx(Nx,Ny),Viscyy(Nx,Ny),Source(Nx,Ny)
real ( kind = 8 ) dn_dx(Nx,Ny),dn_dy(Nx,Ny),dux_dx(Nx,Ny),dux_dy(Nx,Ny),duy_dx(Nx,Ny),duy_dy(Nx,Ny)
real ( kind = 8 ) d2ux_dx2(Nx,Ny),d2ux_dy2(Nx,Ny),d2uy_dx2(Nx,Ny),d2uy_dy2(Nx,Ny)
real ( kind = 8 ) NLn(Nx,Ny),NLx(Nx,Ny),NLy(Nx,Ny)
complex ( kind = 8 ) nk(Nh,Ny),nk_dum(Nh,Ny),ukx(Nh,Ny),ukx_dum(Nh,Ny),uky(Nh,Ny),uky_dum(Nh,Ny)
complex ( kind = 8 ) i_kx_nk(Nh,Ny),i_ky_nk(Nh,Ny)
complex ( kind = 8 ) i_kx_ukx(Nh,Ny),i_ky_ukx(Nh,Ny),i_kx_uky(Nh,Ny),i_ky_uky(Nh,Ny)
complex ( kind = 8 ) kx2_ukx(Nh,Ny),ky2_ukx(Nh,Ny),kx2_uky(Nh,Ny),ky2_uky(Nh,Ny)
complex ( kind = 8 ) Pkx(Nh,Ny),Pky(Nh,Ny),Visckxx(Nh,Ny),Visckxy(Nh,Ny),Visckyx(Nh,Ny),Visckyy(Nh,Ny)
complex ( kind = 8 ) NLkn(Nh,Ny),NLkx(Nh,Ny),NLky(Nh,Ny)
complex ( kind = 8 ) dn_dt_old(Nh,Ny),dn_dt_new(Nh,Ny)
complex ( kind = 8 ) dux_dt_old(Nh,Ny),duy_dt_old(Nh,Ny),dux_dt_new(Nh,Ny),duy_dt_new(Nh,Ny)

integer ( kind = 8 ) plan_forward,plan_backward

common/comm/Lx,Ly,mux,muy,ms,CS,dt

!$OMP PARALLEL SHARED(Source,NLkx,NLky,nk,ukx,uky,nk_dum,ukx_dum,uky_dum) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
    ! Declare the Source Term for Density Equation.
    Source(i,j) = 0.0d0/ms
    ! Set the variables to Zero.
    NLkx(i,j) = 0.0d0
    NLky(i,j) = 0.0d0
    ! Keep Backup of the Arrays for FFTW.
    nk_dum(i,j) = nk(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

! Evaluate the Derivatives in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,nk,ukx,uky,i_kx_nk,i_ky_nk,i_kx_ukx,i_ky_ukx,i_kx_uky,i_ky_uky),&
!$OMP & SHARED(kx2_ukx,ky2_ukx,kx2_uky,ky2_uky) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    i_kx_nk(i,j) = (0.0d0,1.0d0)*kx*nk(i,j)
    i_ky_nk(i,j) = (0.0d0,1.0d0)*ky*nk(i,j)
    i_kx_ukx(i,j) = (0.0d0,1.0d0)*kx*ukx(i,j)
    i_ky_ukx(i,j) = (0.0d0,1.0d0)*ky*ukx(i,j)
    i_kx_uky(i,j) = (0.0d0,1.0d0)*kx*uky(i,j)
    i_ky_uky(i,j) = (0.0d0,1.0d0)*ky*uky(i,j)
    kx2_ukx(i,j) = kx*kx*ukx(i,j)
    ky2_ukx(i,j) = ky*ky*ukx(i,j)
    kx2_uky(i,j) = kx*kx*uky(i,j)
    ky2_uky(i,j) = ky*ky*uky(i,j)
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    i_kx_nk(i,j) = (0.0d0,1.0d0)*kx*nk(i,j)
    i_ky_nk(i,j) = (0.0d0,1.0d0)*ky*nk(i,j)
    i_kx_ukx(i,j) = (0.0d0,1.0d0)*kx*ukx(i,j)
    i_ky_ukx(i,j) = (0.0d0,1.0d0)*ky*ukx(i,j)
    i_kx_uky(i,j) = (0.0d0,1.0d0)*kx*uky(i,j)
    i_ky_uky(i,j) = (0.0d0,1.0d0)*ky*uky(i,j)
    kx2_ukx(i,j) = kx*kx*ukx(i,j)
    ky2_ukx(i,j) = ky*ky*ukx(i,j)
    kx2_uky(i,j) = kx*kx*uky(i,j)
    ky2_uky(i,j) = ky*ky*uky(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, nk_dum, n, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_nk, dn_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_nk, dn_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_ukx, dux_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_ukx, dux_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_uky, duy_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_uky, duy_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, kx2_ukx, d2ux_dx2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ky2_ukx, d2ux_dy2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, kx2_uky, d2uy_dx2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ky2_uky, d2uy_dy2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

!$OMP PARALLEL SHARED(n,ux,uy,dn_dx,dn_dy,dux_dx,dux_dy,duy_dx,duy_dy),&
!$OMP & SHARED(d2ux_dx2,d2ux_dy2,d2uy_dx2,d2uy_dy2),&
!$OMP & SHARED(Px,Py,Viscxx,Viscxy,Viscyx,Viscyy,NLn,NLx,NLy,Source) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
    ! FFTW Normalisation.
    n(i,j) = n(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    dn_dx(i,j) = dn_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    dn_dy(i,j) = dn_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    dux_dx(i,j) = dux_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    dux_dy(i,j) = dux_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    duy_dx(i,j) = duy_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    duy_dy(i,j) = duy_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    d2ux_dx2(i,j) = d2ux_dx2(i,j)/(dfloat(Nx)*dfloat(Ny))
    d2ux_dy2(i,j) = d2ux_dy2(i,j)/(dfloat(Nx)*dfloat(Ny))
    d2uy_dx2(i,j) = d2uy_dx2(i,j)/(dfloat(Nx)*dfloat(Ny))
    d2uy_dy2(i,j) = d2uy_dy2(i,j)/(dfloat(Nx)*dfloat(Ny))
    ! Evaluate Viscous Terms in Real Space.
    Viscxx(i,j) = d2ux_dx2(i,j)/n(i,j)
    Viscxy(i,j) = d2ux_dy2(i,j)/n(i,j)
    Viscyx(i,j) = d2uy_dx2(i,j)/n(i,j)
    Viscyy(i,j) = d2uy_dy2(i,j)/n(i,j)
    ! Evaluate Pressure in Real Space.
    Px(i,j) = dn_dx(i,j)/n(i,j)
    Py(i,j) = dn_dy(i,j)/n(i,j)
    ! Evaluate Non-Linear Terms in Real Space.
    ! Density Equation.
    NLn(i,j) = ux(i,j)*dn_dx(i,j) + uy(i,j)*dn_dy(i,j) + n(i,j)*(dux_dx(i,j) + duy_dy(i,j)) !+ Source(i,j)
    ! Momentum Equation.
    NLx(i,j) = ux(i,j)*dux_dx(i,j) + uy(i,j)*dux_dy(i,j)
    NLy(i,j) = ux(i,j)*duy_dx(i,j) + uy(i,j)*duy_dy(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Px, Pkx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Py, Pky, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Viscxx, Visckxx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Viscxy, Visckxy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Viscyx, Visckyx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Viscyy, Visckyy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, NLn, NLkn, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, NLx, NLkx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, NLy, NLky, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

! De - Aliazing Technique With 2/3 Rule for All the Non-Linear Terms.

!$OMP PARALLEL SHARED(pi,Lx,Ly,NLkn,NLkx,NLky),&
!$OMP & SHARED(Pkx,Pky,Visckxx,Visckxy,Visckyx,Visckyy) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 0) then
      NLkn(i,j) = 0.0d0
      NLkx(i,j) = 0.0d0
      NLky(i,j) = 0.0d0
      Pkx(i,j) = 0.0d0
      Pky(i,j) = 0.0d0
      Visckxx(i,j) = 0.0d0
      Visckxy(i,j) = 0.0d0
      Visckyx(i,j) = 0.0d0
      Visckyy(i,j) = 0.0d0
      endif
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 0) then
      NLkn(i,j) = 0.0d0
      NLkx(i,j) = 0.0d0
      NLky(i,j) = 0.0d0
      Pkx(i,j) = 0.0d0
      Pky(i,j) = 0.0d0
      Visckxx(i,j) = 0.0d0
      Visckxy(i,j) = 0.0d0
      Visckyx(i,j) = 0.0d0
      Visckyy(i,j) = 0.0d0
      endif
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(mux,muy,CS,Pkx,Pky,Visckxx,Visckxy,Visckyx,Visckyy),&
!$OMP & SHARED(NLkn,NLkx,NLky,dn_dt_new,dux_dt_new,duy_dt_new) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
   ! Evaluate Full Non-Linear Term for Momentum Equation in Spectral Space.
    Nlkx(i,j) = NLkx(i,j) + (mux*Visckxx(i,j) + muy*Visckxy(i,j)) + (CS*CS)*Pkx(i,j)
    NLky(i,j) = NLky(i,j) + (mux*Visckyx(i,j) + muy*Visckyy(i,j)) + (CS*CS)*Pky(i,j)
    ! Generate Input for Time Solvers.
    dn_dt_new(i,j) = -NLkn(i,j)
    dux_dt_new(i,j) = -NLkx(i,j)
    duy_dt_new(i,j) = -NLky(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine derive2coh
