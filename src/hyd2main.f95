subroutine tara2dHydSer()

use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

real (kind=8), parameter :: pi = 3.14159265358979323846d0

include "fftw3.f"

integer (kind=4) Np,Nx,Ny,Nh

integer (kind=4) i,j,m,n,t
real (kind = 8) Lx, Ly, dx, dy, h, kx, ky, nu, W0, m0, d, y_Energy
real (kind = 8) time, time_min, time_max, dt

real (kind = 8), dimension (:), allocatable :: xp, yp, vxp, vyp
real (kind = 8), dimension (:), allocatable :: x, y

real (kind=8), dimension(:, :), allocatable :: omega, psi, ux, uy
complex (kind=8), dimension(:, :), allocatable :: omegak, omegak_new, psik, ukx, uky
complex (kind=8), dimension(:, :), allocatable :: omegak_dum,psik_dum,ukx_dum,uky_dum
complex (kind=8), dimension(:, :), allocatable :: dt_omegak_new,  dt_omegak_old

integer (kind=8) plan_forward, plan_backward

character (len=90) :: filename

cfg = parse_cfg("input.ini")

call cfg%get("particle","Np",Np)
call cfg%get("grid","Nx",Nx)
call cfg%get("grid","Ny",Ny)

Nh = (Nx/2) + 1

allocate(xp(Np))
allocate(yp(Np))
allocate(vxp(Np))
allocate(vyp(Np))

allocate(x(Nx))
allocate(y(Ny))

allocate(omega(Nx,Ny))
allocate(psi(Nx,Ny))
allocate(ux(Nx,Ny))
allocate(uy(Nx,Ny))

allocate(omegak(Nh,Ny))
allocate(omegak_new(Nh,Ny))
allocate(psik(Nh,Ny))
allocate(ukx(Nh,Ny))
allocate(uky(Nh,Ny))

allocate(omegak_dum(Nh,Ny))
allocate(psik_dum(Nh,Ny))
allocate(ukx_dum(Nh,Ny))
allocate(uky_dum(Nh,Ny))

allocate(dt_omegak_new(Nh,Ny))
allocate(dt_omegak_old(Nh,Ny))

open(unit = 20, file = 'output/Energy.dat', status = 'unknown')

!====================== USER INPUTS ====================================

call cfg%get("length","Lx",Lx)
call cfg%get("length","Ly",Ly)

call cfg%get("resistivity","nu",nu)

call cfg%get("time","time_min",time_min)
call cfg%get("time","time_max",time_max)
call cfg%get("time","dt",dt)

call cfg%get("initialProfile","W0",W0)
call cfg%get("initialProfile","m0",m0)
call cfg%get("initialProfile","d",d)

Lx = Lx*pi
Ly = Ly*pi

dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)
h = 1.0d0/(dx*dy)

d = d*pi

do i = 1,Np
  xp(i) = rand()*Lx
  yp(i) = rand()*Ly
enddo

do i = 1,Nx
  x(i) = (i-1)*dx
  do j = 1,Ny
    y(j) = (j-1)*dy
    omega(i,j) = W0/dcosh((y(j)+0.5d0*pi)/d)**2.0 - W0/dcosh((y(j)-0.5d0*pi)/d)**2.0
    omega(i,j) = omega(i,j) - W0/dcosh((y(j)+1.50d0*pi)/d)**2.0+W0/dcosh((y(j)-1.50d0*pi)/d)**2.0
    omega(i,j) = omega(i,j) + 0.10d0*dcos(m0*x(i))
  enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega, omegak, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

!==================== MAIN PROGRAM ========================================

do time = time_min, time_max, dt

t = nint(time/dt) - int(time_min/dt)

  call derive (Nx, Ny, Nh, nu, time, omegak, omegak_new, dt_omegak_old, dt_omegak_new)
  call ab (Nx, Ny, Nh, nu, time, dt, omegak, omegak_new, dt_omegak_old, dt_omegak_new)

do i = 1, Nh
  do j = 1, Ny
    dt_omegak_old(i,j) = dt_omegak_new(i,j)
    omegak(i,j) = omegak_new(i,j)
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

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

do i = 1,Nx
  do j = 1,Ny
  omega(i,j) = omega(i,j)/dfloat(Nx*Ny)
  ux(i,j) = ux(i,j)/dfloat(Nx*Ny)
  uy(i,j) = uy(i,j)/dfloat(Nx*Ny)
  if (mod(dfloat(t),1000.0d0) == 0.0d0) then
  write(filename, '("output/fort.",I8.8)') t+100
  open(unit=t+100,file=filename,status='unknown')
  write(t+100,*) x(i), y(j), omega(i,j), ux(i,j), uy(i,j)
  endif
  enddo
enddo

if (mod(dfloat(t),1000.0d0) == 0.0d0) then
  close(t+100)
endif

  do i = 1,Np

  m = aint(xp(i)/dx) + 1

  n = aint(yp(i)/dy) + 1

  if (m /= Nx .and. n /= Ny) then

  vxp(i) = ux(m,n) * (x(m+1) - xp(i)) * (y(n+1) - yp(i)) * h &
         + ux(m+1,n) * (xp(i) - x(m)) * (y(n+1) - yp(i)) * h &
         + ux(m,n+1) * (x(m+1) - xp(i)) * (yp(i) - y(n)) * h &
         + ux(m+1,n+1) * (xp(i) - x(m)) * (yp(i) - y(n)) * h

  vyp(i) = uy(m,n) * (x(m+1) - xp(i)) * (y(n+1) - yp(i)) * h &
         + uy(m+1,n) * (xp(i) - x(m)) * (y(n+1) - yp(i)) * h &
         + uy(m,n+1) * (x(m+1) - xp(i)) * (yp(i) - y(n)) * h &
         + uy(m+1,n+1) * (xp(i) - x(m)) * (yp(i) - y(n)) * h

  endif

  if (m == Nx .and. n == Ny) then

  vxp(i) = ux(m,n) * (Lx - xp(i)) * (Ly - yp(i)) * h &
         + ux(m+1,n) * (xp(i) - x(m)) * (Ly - yp(i)) * h &
         + ux(m,n+1) * (Lx - xp(i)) * (yp(i) - y(n)) * h &
         + ux(m+1,n+1) * (xp(i) - x(m)) * (yp(i) - y(n)) * h

  vyp(i) = uy(m,n) * (Lx - xp(i)) * (Ly - yp(i)) * h &
         + uy(m+1,n) * (xp(i) - x(m)) * (Ly - yp(i)) * h &
         + uy(m,n+1) * (Lx - xp(i)) * (yp(i) - y(n)) * h &
         + uy(m+1,n+1) * (xp(i) - x(m)) * (yp(i) - y(n)) * h

  endif

  if (m == Nx .and. n /= Ny) then

  vxp(i) = ux(m,n) * (Lx - xp(i)) * (y(n+1) - yp(i)) * h &
         + ux(m+1,n) * (xp(i) - x(m)) * (y(n+1) - yp(i)) * h &
         + ux(m,n+1) * (Lx - xp(i)) * (yp(i) - y(n)) * h &
         + ux(m+1,n+1) * (xp(i) - x(m)) * (yp(i) - y(n)) * h

  vyp(i) = uy(m,n) * (Lx - xp(i)) * (y(n+1) - yp(i)) * h &
         + uy(m+1,n) * (xp(i) - x(m)) * (y(n+1) - yp(i)) * h &
         + uy(m,n+1) * (Lx - xp(i)) * (yp(i) - y(n)) * h &
         + uy(m+1,n+1) * (xp(i) - x(m)) * (yp(i) - y(n)) * h

  endif

  if (m /= Nx .and. n == Ny) then

  vxp(i) = ux(m,n) * (x(m+1) - xp(i)) * (Ly - yp(i)) * h &
         + ux(m+1,n) * (xp(i) - x(m)) * (Ly - yp(i)) * h &
         + ux(m,n+1) * (x(m+1) - xp(i)) * (yp(i) - y(n)) * h &
         + ux(m+1,n+1) * (xp(i) - x(m)) * (yp(i) - y(n)) * h

  vyp(i) = uy(m,n) * (x(m+1) - xp(i)) * (Ly - yp(i)) * h &
         + uy(m+1,n) * (xp(i) - x(m)) * (Ly - yp(i)) * h &
         + uy(m,n+1) * (x(m+1) - xp(i)) * (yp(i) - y(n)) * h &
         + uy(m+1,n+1) * (xp(i) - x(m)) * (yp(i) - y(n)) * h

  endif

  enddo

  do i = 1,Np
  xp(i) = xp(i) + vxp(i) * dt
  yp(i) = yp(i) + vyp(i) * dt

  xp(i) = xp(i) - (int(xp(i)/Lx)) * Lx
    if (xp(i) .le. 0.0d0) then
      xp(i) = xp(i) + Lx
    endif

  yp(i) = yp(i) - (int(yp(i)/Ly)) * Ly
    if (yp(i) .le. 0.0d0) then
      yp(i) = yp(i) + Ly
    endif

  enddo

  if (mod(dfloat(t),1000.0d0) == 0.0d0) then
  do i = 1,Np
  write(filename, '("output/fort.",I8.8)') t+105
  open(unit=t+105,file=filename,status='unknown')
  write(t+105,*) xp(i),yp(i)
  enddo
  endif

enddo !time

end subroutine tara2dHydSer
