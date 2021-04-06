subroutine tara1dHydSer(arg)

use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

real (kind=8), parameter :: pi = 3.14159265358979323846d0

include "fftw3.f"

integer (kind=4) N,Nh

real (kind = 8), dimension (:), allocatable :: x, ux, dux, ddux

complex (kind=8), dimension(:), allocatable :: uk,duk,dduk

integer*8 :: plan
real*8 :: L

integer :: i,t
real*8 dx,k,time,time_min,time_max,dt,nu
real (kind = 8), dimension (:), allocatable :: force_u_n, force_u_o

character (len=90) :: filename
character (len=32) :: arg

cfg = parse_cfg(arg)

call cfg%get("grid","Nx",N)

Nh = (N/2) + 1

allocate(x(N))
allocate(ux(N))
allocate(dux(N))
allocate(ddux(N))

allocate(force_u_n(N))
allocate(force_u_o(N))

allocate(uk(Nh))
allocate(duk(Nh))
allocate(dduk(Nh))

call cfg%get("length","Lx",L)

call cfg%get("resistivity","nu",nu)

call cfg%get("time","time_min",time_min)
call cfg%get("time","time_max",time_max)
call cfg%get("time","dt",dt)

L  = L*pi
dx = L/dfloat(N)

do i = 1,N
  x(i) = dfloat(i-1)*dx
  ux(i) = dsin(x(i))
enddo

do time = time_min, time_max, dt

t = nint(time/dt) - int(time_min/dt)

call dfftw_plan_dft_r2c_1d(plan, N, ux, uk, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, ux, uk)
call dfftw_destroy_plan(plan)

do i = 1,Nh
  k = 2.0d0*pi*dfloat(i-1)/L
  duk(i) = (0.0d0,1.0d0) * k * uk(i)
  dduk(i) = - k * k * uk(i)
enddo

call dfftw_plan_dft_c2r_1d(plan, N, duk, dux, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, duk, dux)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_c2r_1d(plan, N, dduk, ddux, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, dduk, ddux)
call dfftw_destroy_plan(plan)

do i = 1,N
dux(i) = dux(i)/dfloat(N)
ddux(i) = ddux(i)/dfloat(N)
  if (t .ge. 1000 .and. mod(t,500) == 0) then
    write(filename, '("output/fort.",I8.8)') t
    open(unit=t,file=filename,status='unknown')
    write(t,*) t, x(i), ux(i)
  endif
enddo

do i = 1,N
  force_u_n(i) = -ux(i)*dux(i) + nu*ddux(i)
    if (t==0) then
    ux(i) = ux(i) + dt * force_u_n(i)
    else
    ux(i) = ux(i) + dt* ( (3.0/2.0) * force_u_n(i) - (1.0/2.0) * force_u_o(i))
    endif
 force_u_o(i) = force_u_n(i)
enddo

enddo

end subroutine tara1dHydSer
