subroutine tara2dScrSer(arg)

use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3.f"

integer (kind=4) Np,Nx,Ny,Nh

integer ( kind = 4 ) i,j,t
real ( kind = 8 ) Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,alpha,nu,lambda,C0,C1,C2,C3,W0,m,d,G,Gb
REAL ( kind = 8 ) Energy,Energy0,yEnergy,yEnergy0
real (kind = 8), dimension (:), allocatable :: x,y
real (kind=8), dimension(:, :), allocatable :: psi,omega,omega_dum,ux,uy,ux_dum,uy_dum
complex (kind=8), dimension(:, :), allocatable :: omegak,omegak_dum,omegak_new,dt_omegak_old,dt_omegak_new
complex (kind=8), dimension(:, :), allocatable :: psik,psik_dum,ukx,uky,ukx_dum,uky_dum
integer ( kind = 8 ) plan_forward,plan_backward

character (len=90) :: filename
character (len=32) :: arg

cfg = parse_cfg(arg)

call cfg%get("particle","Np",Np)
call cfg%get("grid","Nx",Nx)
call cfg%get("grid","Ny",Ny)

Nh = (Nx/2) + 1


allocate(x(Nx))
allocate(y(Ny))

allocate(psi(Nx,Ny))
allocate(omega(Nx,Ny))
allocate(omega_dum(Nx,Ny))
allocate(ux(Nx,Ny))
allocate(uy(Nx,Ny))
allocate(ux_dum(Nx,Ny))
allocate(uy_dum(Nx,Ny))

allocate(omegak(Nh,Ny))
allocate(omegak_dum(Nh,Ny))
allocate(omegak_new(Nh,Ny))
allocate(dt_omegak_old(Nh,Ny))
allocate(dt_omegak_new(Nh,Ny))
allocate(psik(Nh,Ny))
allocate(psik_dum(Nh,Ny))
allocate(ukx(Nh,Ny))
allocate(uky(Nh,Ny))
allocate(ukx_dum(Nh,Ny))
allocate(uky_dum(Nh,Ny))

open(unit=10,file='output/Initial_Condition_Omega_REP.dat',status='unknown')
open(unit=15,file='output/Initial_Condition_Velocity.dat',status='unknown')
open(unit=40,file='output/Energy.dat',status='unknown')

!===================== USER INPUTS ============================================

call cfg%get("length","Lx",Lx)
call cfg%get("length","Ly",Ly)

call cfg%get("resistivity","nu",nu)

call cfg%get("time","time_min",time_min)
call cfg%get("time","time_max",time_max)
call cfg%get("time","dt",dt)

call cfg%get("initialProfile","W0",W0)
call cfg%get("initialProfile","m0",m)
call cfg%get("initialProfile","d",d)

Lx = Lx*pi
Ly = Ly*pi

dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

d = d*pi

alpha = 1.0d0
lambda = 10.0d0

do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0+real(j-1)*dy
    ux(i,j) = 0.0d0
    uy(i,j) = 0.0d0
    !read(25,*) x(i),y(j),omega(i,j)
    !omega(i,j) = 2.0d0*sin(x(i))*cos(y(j))
    !read(50,*) x(i),y(j),omega(i,j)
    omega(i,j) = W0/dcosh((y(j)+0.50d0*pi)/d)**2.0-W0/dcosh((y(j)-0.50d0*pi)/d)**2.0
    omega(i,j) = omega(i,j) - W0/dcosh((y(j)+1.50d0*pi)/d)**2.0+W0/dcosh((y(j)-1.50d0*pi)/d)**2.0
    omega(i,j) = omega(i,j)+0.01*dcos(m*x(i))
    omega_dum(i,j) = omega(i,j)
    write(15,*) x(i),y(j),omega(i,j),ux(i,j),uy(i,j)
  end do
end do

!===================== INITIAL TIME DATA ===============================================

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega_dum, omegak, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak(i,j)/( kx*kx + ky*ky )
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
      endif
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/( kx*kx + ky*ky )
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
    psik_dum(i,j) = psik(i,j)
    omegak_dum(i,j) = omegak(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, psik_dum, psi, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

Energy0 = 0.0d0
 C0 = 0.0d0
 C1 = 0.0d0
 C2 = 0.0d0
 C3 = 0.0d0

do i = 1, Nx
  do j = 1, Ny
  psi(i,j) = psi(i,j)/(dfloat(Nx)*dfloat(Ny))
  omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
  ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
  uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
  Energy0 = Energy0 + (ux(i,j)**2 + uy(i,j)**2)
  yEnergy0 = yEnergy0 + (uy(i,j)**2)
  C0 = C0 + omega(i,j)**0.0d0
  C1 = C1 + omega(i,j)**1.0d0
  C2 = C2 + omega(i,j)**2.0d0
  C3 = C3 + omega(i,j)**3.0d0
  write(10,*) x(i),y(j),omega(i,j),psi(i,j),ux(i,j),uy(i,j)
  end do
end do

!write(40,*) 0, Energy0, yEnergy0, C0,C1,C2,C3

!======================= MAIN PROGRAM =====================================================

do time = time_min,time_max,dt

t = nint(time/dt) - int(time_min/dt)

!if (mod(t,100) == 0) then
!print*, "time = ", time!, t/100
!endif

!====================================================================================
  call derive2css (Nx,Ny,Nh,alpha,nu,time,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
!  call rk42css (Nx,Ny,Nh,alpha,nu,time,dt,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
  call ab2css (Nx,Ny,Nh,alpha,nu,time,dt,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
!====================================================================================

do i = 1,Nx/2+1
  do j = 1,Ny
    dt_omegak_old(i,j) = dt_omegak_new(i,j)
    omegak(i,j) = omegak_new(i,j)
  enddo
enddo

do i = 1,Nx/2+1
  do j = 1,Ny/2
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat(j-1)/Ly
      if (i == 1 .and. j == 1) then
      psik(i,j) = (0.0d0,0.0d0)
      ukx(i,j) = (0.0d0,0.0d0)
      uky(i,j) = (0.0d0,0.0d0)
      else
      psik(i,j) = omegak(i,j)/( kx*kx + ky*ky + lambda)
      ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
      uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
      endif
    omegak_dum(i,j) = (kx*kx + ky*ky) * psik(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
  do j = Ny/2+1,Ny
    kx = 2.0d0*pi*dfloat(i-1)/Lx
    ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
    psik(i,j) = omegak(i,j)/( kx*kx + ky*ky )
    ukx(i,j) = + (0.0d0,1.0d0) * ky * psik(i,j)
    uky(i,j) = - (0.0d0,1.0d0) * kx * psik(i,j)
    omegak_dum(i,j) = (kx*kx + ky*ky) * psik(i,j)
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

Energy = 0.0d0
 C0 = 0.0d0
 C1 = 0.0d0
 C2 = 0.0d0
 C3 = 0.0d0

do i = 1,Nx
  do j = 1,Ny
  omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
  ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
  uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
  Energy = Energy + (ux(i,j)**2 + uy(i,j)**2)
  yEnergy = yEnergy + (uy(i,j)**2)
  C0 = C0 + omega(i,j)**0.0d0
  C1 = C1 + omega(i,j)**1.0d0
  C2 = C2 + omega(i,j)**2.0d0
  C3 = C3 + omega(i,j)**3.0d0
    if (mod(dfloat(t),100.0) == 0.0) then
      write(filename, '("output/fort.",I8.8)') t+100
      open(unit=t+100,file=filename,status='unknown')
      write(t+100,*) x(i),y(j),ux(i,j),uy(i,j),omega(i,j)
    endif
  enddo
enddo
    if (mod(dfloat(t),100.0) == 0.0) then
    close(t+100)
    endif

WRITE(40,*) time,Energy, yEnergy!,C0,C1,C2,C3

enddo ! time

!====================================================================================

end subroutine tara2dScrSer
