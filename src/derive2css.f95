
!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive2css(Nx,Ny,Nh,alpha,nu,time,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
implicit none

include "fftw3.f"

integer ( kind = 4 ) Nx,Ny,Nh, i, j
real ( kind = 8 ) time,alpha,nu,lambda,kx,ky,Lx,Ly
real ( kind = 8 ) ux(Nx,Ny),uy(Nx,Ny),domega_dx(Nx,Ny),domega_dy(Nx,Ny)
real ( kind = 8 ) ux_domega_dx(Nx,Ny),uy_domega_dy(Nx,Ny)
complex ( kind = 8 ) ukx(Nh,Ny),ukx_dum(Nh,Ny),uky(Nh,Ny),uky_dum(Nh,Ny),omegak_dum(Nh,Ny),omegak_dummy(Nh,Ny)
complex ( kind = 8 ) psik(Nh,Ny),omegak(Nh,Ny),i_kx_omegak(Nh,Ny),i_ky_omegak(Nh,Ny),omegak_new(Nh,Ny)
complex ( kind = 8 ) NLkx1(Nh,Ny),NLky1(Nh,Ny),NLk1(Nh,Ny),dt_omegak_old(Nh,Ny),dt_omegak_new(Nh,Ny)
real (kind=8), parameter :: pi = 3.1428

integer ( kind = 8 ) plan_forward,plan_backward

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
     omegak_dummy(i,j) = omegak(i,j)
     if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak_dummy(i,j)/( kx*kx + ky*ky +lambda)
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
      endif
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
    omegak_dum(i,j) = (kx*kx + ky*ky) * psik(i,j)
    i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak_dum(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak_dum(i,j)
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    omegak_dummy(i,j) = omegak(i,j)
    psik(i,j) = omegak_dummy(i,j)/( kx*kx + ky*ky + lambda)
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
    omegak_dum(i,j) = (kx*kx + ky*ky) * psik(i,j)
    i_kx_omegak(i,j) = (0.0d0,1.0d0)*kx*omegak_dum(i,j)
    i_ky_omegak(i,j) = (0.0d0,1.0d0)*ky*omegak_dum(i,j)
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
    ux_domega_dx(i,j) = ux(i,j)*domega_dx(i,j)
    uy_domega_dy(i,j) = uy(i,j)*domega_dy(i,j)
  enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_domega_dx, NLkx1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_domega_dy, NLky1, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      ! De - Aliazing Technique...
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 0) then
      NLkx1(i,j) = 0.0d0
      NLky1(i,j) = 0.0d0
      endif
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      ! De - Aliazing Technique...
      if (dsqrt(kx*kx + ky*ky) .ge. (dfloat(Nx+Ny)/2.0)/3.0 + 0) then
      NLkx1(i,j) = 0.0d0
      NLky1(i,j) = 0.0d0
      endif
  enddo
enddo

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
    NLk1(i,j) = alpha * ( NLkx1(i,j) + NLky1(i,j) )
    NLk1(i,j) = NLk1(i,j) + nu * (kx*kx + ky*ky)*omegak_dum(i,j)
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    NLk1(i,j) = alpha * ( NLkx1(i,j) + NLky1(i,j) )
    NLk1(i,j) = NLk1(i,j) + nu * (kx*kx + ky*ky)*omegak_dum(i,j)
  enddo
enddo

do i = 1,Nx/2+1
  do j = 1,Ny
    dt_omegak_new(i,j) = -NLk1(i,j)
  enddo
enddo

return

end subroutine derive2css
