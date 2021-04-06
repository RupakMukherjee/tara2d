!=================== SUBROUTINE DERIVE ===================================

subroutine derive2csh (Nx, Ny, Nh, nu, time, omegak, omegak_new, dt_omegak_old, dt_omegak_new)
implicit none
include "fftw3.f"
integer (kind = 4) Nx, Ny, Nh, i, j
real (kind = 8) time, kx, ky, nu, Lx, Ly
real (kind = 8) ux(Nx,Ny), uy (Nx,Ny)
real (kind = 8) domega_dx(Nx,Ny), domega_dy(Nx,Ny)
real (kind = 8) ux_domega_dx(Nx,Ny),uy_domega_dy(Nx,Ny)
complex (kind=8) NLkx(Nh,Ny),NLky(Nh,Ny), NLk(Nh,Ny)
complex (kind=8) omegak(Nh,Ny), omegak_new(Nh,Ny), psik(Nh,Ny), ukx(Nh,Ny), uky(Nh,Ny)
complex (kind=8) i_kx_omegak(Nh,Ny),i_ky_omegak(Nh,Ny)
complex (kind=8) omegak_dum(Nh,Ny),psik_dum(Nh,Ny),ukx_dum(Nh,Ny),uky_dum(Nh,Ny)
complex (kind=8) dt_omegak_new(Nh, Ny),  dt_omegak_old(Nh, Ny)
real (kind=8), parameter :: pi = 3.1428

integer (kind=8) plan_forward, plan_backward

Lx = 2.0d0*pi
Ly = 2.0d0*pi

do i = 1, Nh
  do j = 1, Ny
  omegak_dum(i,j) = omegak(i,j)
  enddo
enddo

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak(i,j)/(kx*kx+ky*ky)
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
      endif
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/(kx*kx+ky*ky)
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
enddo

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0) * kx * omegak(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0) * ky * omegak(i,j)
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    i_kx_omegak(i,j) = (0.0d0,1.0d0) * kx * omegak(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0) * ky * omegak(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_kx_omegak, domega_dx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, i_ky_omegak, domega_dy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

do i = 1,Nx
  do j = 1,Ny
    ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
    uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_dx(i,j) = domega_dx(i,j)/(dfloat(Nx)*dfloat(Ny))
    domega_dy(i,j) = domega_dy(i,j)/(dfloat(Nx)*dfloat(Ny))
    ux_domega_dx(i,j) = ux(i,j) * domega_dx(i,j)
    uy_domega_dy(i,j) = uy(i,j) * domega_dy(i,j)
  enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_domega_dx, NLkx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_domega_dy, NLky, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

! De-Aliazing Technique - 2/3 Truncation...
do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    if ( dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 ) then
    NLkx(i,j) = 0.0d0
    NLky(i,j) = 0.0d0
    endif
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    if ( dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 ) then
    NLkx(i,j) = 0.0d0
    NLky(i,j) = 0.0d0
    endif
  enddo
enddo

do i = 1,Nh
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  do j = 1,Ny/2
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    NLk(i,j) = NLkx(i,j) + NLky(i,j)
    NLk(i,j) = NLk(i,j) + nu * (kx*kx + ky*ky) * omegak(i,j)
  enddo
  do j = Ny/2+1,Ny
    ky = 2.0d0*pi*(dfloat(j-1)-Ny)/Ly
    NLk(i,j) = NLkx(i,j) + NLky(i,j)
    NLk(i,j) = NLk(i,j) + nu * (kx*kx + ky*ky) * omegak(i,j)
  enddo
enddo

do i = 1,Nh
  do j = 1,Ny
    dt_omegak_new(i,j) = -NLk(i,j)
  enddo
enddo

return

end subroutine derive2csh
