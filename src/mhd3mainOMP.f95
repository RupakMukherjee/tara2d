subroutine tara3dMHDOMP(arg)

use omp_lib
use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3.f"

integer (kind=4) Np,Nx,Ny,Nz,Nh

integer ( kind = 4 ) i,j,k,l,n,o,t,P_s,P_max,q,q_max
real ( kind = 8 ) Lx,Ly,Lz,dx,dy,dz,dr,delta,kx,ky,kz,time,time_min,time_max,dt,G1,G2,G3,A,B,C
real ( kind = 8 ) spheat,rho0,ms,U0,mu,kf,mode,sigma,uy0,M,CS,Re,Rm,PM,P0,eta,mu_0,MA,VA,B0,ba
real ( kind = 8 ) Pressure,Energy,Energy0,T_Energy,y_Energy0,B_Field,C0,C1,C2,C3,div_B_tot,Hf,Hm,HT

real (kind = 8), dimension (:), allocatable :: x,y,z
real (kind=8), dimension(:, :, :), allocatable :: r,rho,ux,uy,uz
real (kind=8), dimension(:, :, :), allocatable :: rho_ux,rho_uy,rho_uz
real (kind=8), dimension(:, :, :), allocatable :: jx,jy,jz,Ax,Ay,Az
real (kind=8), dimension(:, :, :), allocatable :: omega_x,omega_y,omega_z,omega2
real (kind=8), dimension(:, :, :), allocatable :: P,E,Bx,By,Bz,B2,div_B
real (kind=8), dimension(:, :, :), allocatable :: rho_dum,ux_dum,uy_dum,uz_dum
real (kind=8), dimension(:, :, :), allocatable :: rho_ux_dum,rho_uy_dum,rho_uz_dum
real (kind=8), dimension(:, :, :), allocatable :: P_dum,E_dum,Bx_dum,By_dum,Bz_dum
real ( kind = 8 ) d_rho(1024),d_ux(1024),d_uy(1024),d_uz(1024),d_Bx(1024),d_By(1024),d_Bz(1024)
real ( kind = 8 ) S_rho(10,1024),S_u(10,1024),S_B(10,1024)

complex (kind=8), dimension(:, :, :), allocatable :: rho_k,ux_k,uy_k,uz_k
complex (kind=8), dimension(:, :, :), allocatable :: rho_ux_k,rho_uy_k,rho_uz_k
complex (kind=8), dimension(:, :, :), allocatable :: omega_x_k,omega_y_k,omega_z_k
complex (kind=8), dimension(:, :, :), allocatable :: jx_k,jy_k,jz_k,Ax_k,Ay_k,Az_k
complex (kind=8), dimension(:, :, :), allocatable :: P_k,E_k,Ek,Bk
complex (kind=8), dimension(:, :, :), allocatable :: Bx_k,By_k,Bz_k,div_B_k
complex (kind=8), dimension(:, :, :), allocatable :: rho_k_dum,ux_k_dum,uy_k_dum,uz_k_dum
complex (kind=8), dimension(:, :, :), allocatable :: omega_x_k_dum,omega_y_k_dum,omega_z_k_dum
complex (kind=8), dimension(:, :, :), allocatable :: rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum
complex (kind=8), dimension(:, :, :), allocatable :: E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum
complex (kind=8), dimension(:, :, :), allocatable :: rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new
complex (kind=8), dimension(:, :, :), allocatable :: E_k_new,Bx_k_new,By_k_new,Bz_k_new

complex (kind=8), dimension(:, :, :), allocatable :: d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old
complex (kind=8), dimension(:, :, :), allocatable :: d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old
complex (kind=8), dimension(:, :, :), allocatable :: d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new
complex (kind=8), dimension(:, :, :), allocatable :: d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new

integer ( kind = 8 ) plan_forward,plan_backward ! FFTW

integer ( kind = 4 ) iret,thread_num,num_thread,proc_num ! OMP
real ( kind = 8 ) t1,t2

common/comm/iret,thread_num,time_max,Lx,Ly,Lz,dx,dy,dz,spheat,mu,ms,CS,mu_0,eta,dt,kf

integer,parameter :: seed = 99999999

character (len=90) :: filename
character (len=32) :: arg

cfg = parse_cfg(arg)

call cfg%get("grid","Nx",Nx)
call cfg%get("grid","Ny",Ny)
call cfg%get("grid","Nz",Nz)

Nh = (Nx/2) + 1

call srand(seed)

allocate(x(Nx))
allocate(y(Ny))
allocate(z(Nz))

allocate(r(Nx,Ny,Nz))
allocate(rho(Nx,Ny,Nz))
allocate(ux(Nx,Ny,Nz))
allocate(uy(Nx,Ny,Nz))
allocate(uz(Nx,Ny,Nz))
allocate(rho_ux(Nx,Ny,Nz))
allocate(rho_uy(Nx,Ny,Nz))
allocate(rho_uz(Nx,Ny,Nz))
allocate(jx(Nx,Ny,Nz))
allocate(jy(Nx,Ny,Nz))
allocate(jz(Nx,Ny,Nz))
allocate(Ax(Nx,Ny,Nz))
allocate(Ay(Nx,Ny,Nz))
allocate(Az(Nx,Ny,Nz))
allocate(omega_x(Nx,Ny,Nz))
allocate(omega_y(Nx,Ny,Nz))
allocate(omega_z(Nx,Ny,Nz))
allocate(omega2(Nx,Ny,Nz))
allocate(P(Nx,Ny,Nz))
allocate(E(Nx,Ny,Nz))
allocate(Bx(Nx,Ny,Nz))
allocate(By(Nx,Ny,Nz))
allocate(Bz(Nx,Ny,Nz))
allocate(B2(Nx,Ny,Nz))
allocate(div_B(Nx,Ny,Nz))
allocate(rho_dum(Nx,Ny,Nz))
allocate(ux_dum(Nx,Ny,Nz))
allocate(uy_dum(Nx,Ny,Nz))
allocate(uz_dum(Nx,Ny,Nz))
allocate(rho_ux_dum(Nx,Ny,Nz))
allocate(rho_uy_dum(Nx,Ny,Nz))
allocate(rho_uz_dum(Nx,Ny,Nz))
allocate(P_dum(Nx,Ny,Nz))
allocate(E_dum(Nx,Ny,Nz))
allocate(Bx_dum(Nx,Ny,Nz))
allocate(By_dum(Nx,Ny,Nz))
allocate(Bz_dum(Nx,Ny,Nz))
allocate(rho_k(Nh,Ny,Nz))
allocate(ux_k(Nh,Ny,Nz))
allocate(uy_k(Nh,Ny,Nz))
allocate(uz_k(Nh,Ny,Nz))
allocate(rho_ux_k(Nh,Ny,Nz))
allocate(rho_uy_k(Nh,Ny,Nz))
allocate(rho_uz_k(Nh,Ny,Nz))
allocate(omega_x_k(Nh,Ny,Nz))
allocate(omega_y_k(Nh,Ny,Nz))
allocate(omega_z_k(Nh,Ny,Nz))
allocate(jx_k(Nh,Ny,Nz))
allocate(jy_k(Nh,Ny,Nz))
allocate(jz_k(Nh,Ny,Nz))
allocate(Ax_k(Nh,Ny,Nz))
allocate(Ay_k(Nh,Ny,Nz))
allocate(Az_k(Nh,Ny,Nz))
allocate(P_k(Nh,Ny,Nz))
allocate(E_k(Nh,Ny,Nz))
allocate(Ek(Nh,Ny,Nz))
allocate(Bk(Nh,Ny,Nz))
allocate(Bx_k(Nh,Ny,Nz))
allocate(By_k(Nh,Ny,Nz))
allocate(Bz_k(Nh,Ny,Nz))
allocate(div_B_k(Nh,Ny,Nz))
allocate(rho_k_dum(Nh,Ny,Nz))
allocate(ux_k_dum(Nh,Ny,Nz))
allocate(uy_k_dum(Nh,Ny,Nz))
allocate(uz_k_dum(Nh,Ny,Nz))
allocate(omega_x_k_dum(Nh,Ny,Nz))
allocate(omega_y_k_dum(Nh,Ny,Nz))
allocate(omega_z_k_dum(Nh,Ny,Nz))
allocate(rho_ux_k_dum(Nh,Ny,Nz))
allocate(rho_uy_k_dum(Nh,Ny,Nz))
allocate(rho_uz_k_dum(Nh,Ny,Nz))
allocate(E_k_dum(Nh,Ny,Nz))
allocate(Bx_k_dum(Nh,Ny,Nz))
allocate(By_k_dum(Nh,Ny,Nz))
allocate(Bz_k_dum(Nh,Ny,Nz))
allocate(rho_k_new(Nh,Ny,Nz))
allocate(rho_ux_k_new(Nh,Ny,Nz))
allocate(rho_uy_k_new(Nh,Ny,Nz))
allocate(rho_uz_k_new(Nh,Ny,Nz))
allocate(E_k_new(Nh,Ny,Nz))
allocate(Bx_k_new(Nh,Ny,Nz))
allocate(By_k_new(Nh,Ny,Nz))
allocate(Bz_k_new(Nh,Ny,Nz))
allocate(d_rho_k_dt_old(Nh,Ny,Nz))
allocate(d_rho_ux_k_dt_old(Nh,Ny,Nz))
allocate(d_rho_uy_k_dt_old(Nh,Ny,Nz))
allocate(d_rho_uz_k_dt_old(Nh,Ny,Nz))
allocate(d_E_k_dt_old(Nh,Ny,Nz))
allocate(d_Bx_k_dt_old(Nh,Ny,Nz))
allocate(d_By_k_dt_old(Nh,Ny,Nz))
allocate(d_Bz_k_dt_old(Nh,Ny,Nz))
allocate(d_rho_k_dt_new(Nh,Ny,Nz))
allocate(d_rho_ux_k_dt_new(Nh,Ny,Nz))
allocate(d_rho_uy_k_dt_new(Nh,Ny,Nz))
allocate(d_rho_uz_k_dt_new(Nh,Ny,Nz))
allocate(d_E_k_dt_new(Nh,Ny,Nz))
allocate(d_Bx_k_dt_new(Nh,Ny,Nz))
allocate(d_By_k_dt_new(Nh,Ny,Nz))
allocate(d_Bz_k_dt_new(Nh,Ny,Nz))


!===================== FILENAMES ==============================================

open(unit=5,file='output/System_information.dat',status='unknown')
open(unit=15,file='output/Initial_Grid_Data.dat',status='unknown')
!open(unit=20,file='INPUT.dat',status='old')  ! This input file is the file generated from Vorticity code.
!open(unit=25,file='Initial_Grid_Data_Reproduced.dat',status='unknown')
!open(unit=30,file='Initial_Energy_Spectra.dat',status='unknown')
open(unit=35,file='output/Energy_Spectra.dat',status='unknown')
open(unit=40,file='output/Energy.dat',status='unknown')
open(unit=50,file='output/Structure_Function.dat',status='unknown')
open(unit=60,file='output/State.dat',status='unknown')

!===================== USER INPUTS ============================================

! Define Number of Threads.
proc_num = omp_get_num_procs()
thread_num = 40
call omp_set_num_threads (thread_num)

write(5,*) "Number of Processors=", proc_num, "Number of Threads=", thread_num

! System Size.
call cfg%get("length","Lx",Lx)
call cfg%get("length","Ly",Ly)
call cfg%get("length","Lz",Lz)
Lx = Lx*pi; Ly = Ly*pi; Lz = Lz*pi
!Lx = 1.0d0; Ly = 1.0d0; Lz = 1.0d0

! Grid Resolution.
dx = Lx/dfloat(Nx); dy = Ly/dfloat(Ny); dz = Lz/dfloat(Nz)

! Runtime Details and Time Step Width.
call cfg%get("time","time_min",time_min)
call cfg%get("time","time_max",time_max)
call cfg%get("time","dt",dt)

! Ratio of Specific Heats/ Adiabetic Index/ Heat Capacity Ratio/ Poisson Constant.
spheat = 1.0d0

! Co-efficient of Viscosity.
call cfg%get("resistivity","nu",mu)

! Co-efficient of Resistivity.
eta = 0.005d0

! Magnetic Permeability.
mu_0 = 1.0d0

! Mass of the Fluid Element.
ms = 1.0d0

! Background Density.
rho0 = 1.0d0

! Initial Pressure.
P0 = 1.0d0

! Mach Number.
!M = 0.50d0
M = 0.10d0

! Alfven Mach Number.
MA = 1.0d0

! Maximum Velocity.
!U0 = 0.14732247502913884d0
U0 = 0.10d0

! Sound Speed.
! CS = dsqrt(spheat*P0/rho0)
! U0 = M * CS

! Sound Speed.
 CS = 1.0d0/M

! Alfven Speed.
VA = U0/MA

! Initial Magnetic Field.
B0 = VA*dsqrt(mu_0*rho0)

write(5,*) "Sound Speed=", CS, "Initial Velocity=", U0
write(5,*) "Alfven Speed=", VA, "Initial Magnetic Field=", B0

! Forcing Length Scale.
kf = 1.0d0

 A = 1.0d0
 B = 1.0d0
 C = 1.0d0

! Raynold's Number.
Re = 450.0d0 !ms*rho0*U0*Lx/mu; write(5,*) "Raynold's Number =", Re
Rm = 450.0d0 !U0*Lx/eta; write(5,*) "Magnetic Raynold's Number =", Rm
PM = mu/eta; write(5,*) "Prandtl Number =", PM

! Structure Function Order.
P_max = 7

! Maximum Length Scale.
q_max = int(dsqrt(Lx*Lx+Ly*Ly+Lz*Lz)) + 1

! Width in Binning.
delta = 1.0d0

do P_s = 1,P_max,1
  do q = 1,q_max,1
    S_rho(P_s,q) = 0.0d0
    S_u(P_s,q) = 0.0d0
    S_B(P_s,q) = 0.0d0
  enddo !q
enddo !P_s

! Grid Generation.
do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0d0+real(j-1)*dy
    do k = 1, Nz
      z(k)=0.0d0+real(k-1)*dz
      ! Initial Density Distribution.
      rho(i,j,k) = rho0
      ! Initial Velocity Distribution.
     ux(i,j,k) = U0*(A*dsin(kf*z(k)) + C*dcos(kf*y(j)))
     uy(i,j,k) = U0*(B*dsin(kf*x(i)) + A*dcos(kf*z(k)))
     uz(i,j,k) = U0*(C*dsin(kf*y(j)) + B*dcos(kf*x(i)))
      !ux(i,j,k) = 2.0d0*U0/dsqrt(3.0d0) * rand() ! * (dtanh((y(j)-Ly/3.0d0)/a) - dtanh((y(j)-2.0d0*Ly/3.0d0)/a) - 1.0d0)
      !ux(i,j,k) = U0*tanh((y(j)-Ly/2.0d0)/a)
      !uy(i,j,k) = 2.0d0*U0/dsqrt(3.0d0) * rand() !0.0d0
      ! Initial Velocity Perturbation.
      !uy(i,j,k) = uy(i,j,k) + uy0*dsin(mode*x(i))*dexp(-((y(j)-Ly/3.0d0)**2)/(sigma**2))
      !uy(i,j,k) = uy(i,j,k) + uy0*dsin(mode*x(i))*dexp(-((y(j)-2.0d0*Ly/3.0d0)**2)/(sigma**2))
      !uy(i,j,k) = uy(i,j,k) + uy0*dsin(mode*x(i))*exp(-((y(j)-Ly/2.0d0)**2)/(sigma**2))
      !uz(i,j,k) = 2.0d0*U0/dsqrt(3.0d0) * rand() !0.0d0
      ! Initial Velocity Distribution.
      !read(20,*) G1,G2,G3,ux(i,j,k),uy(i,j,k),uz(i,j,k)  ! You may define your own function here.
      ! Initial Pressure Distribution.
      P(i,j,k) = P0
      ! Initial Magnetic Field Profile.
      Bx(i,j,k) = B0!*(A*dsin(kf*z(k)) + C*dcos(kf*y(j)))
      By(i,j,k) = B0!*(B*dsin(kf*x(i)) + A*dcos(kf*z(k)))
      Bz(i,j,k) = B0!*(C*dsin(kf*y(j)) + B*dcos(kf*x(i)))
      !Bx(i,j,k) = 2.0d0*B0/dsqrt(3.0d0) * rand() ! * dtanh((y(j)-Ly/3.0d0)/ba)
      !By(i,j,k) = 2.0d0*B0/dsqrt(3.0d0) * rand() !0.0d0
      !Bz(i,j,k) = 2.0d0*B0/dsqrt(3.0d0) * rand() !0.0d0
      ! Initial Energy Distribution.
      E(i,j,k) = P(i,j,k) + 0.50d0*rho(i,j,k) * (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)&
                 + 0.50d0*(Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)
      ! Initial Combined Variables.
      rho_ux(i,j,k) = rho(i,j,k) * ux(i,j,k)
      rho_uy(i,j,k) = rho(i,j,k) * uy(i,j,k)
      rho_uz(i,j,k) = rho(i,j,k) * uz(i,j,k)
      ! Keep Backup of the Arrays for FFTW.
      rho_dum(i,j,k) = rho(i,j,k)
      ux_dum(i,j,k) = ux(i,j,k)
      uy_dum(i,j,k) = uy(i,j,k)
      uz_dum(i,j,k) = uz(i,j,k)
      rho_ux_dum(i,j,k) = rho_ux(i,j,k)
      rho_uy_dum(i,j,k) = rho_uy(i,j,k)
      rho_uz_dum(i,j,k) = rho_uz(i,j,k)
      P_dum(i,j,k) = P(i,j,k)
      E_dum(i,j,k) = E(i,j,k)
      Bx_dum(i,j,k) = Bx(i,j,k)
      By_dum(i,j,k) = By(i,j,k)
      Bz_dum(i,j,k) = Bz(i,j,k)
      ! Write Initial Density and Velocity Distribution in File.
      !write(15,*) x(i),y(j),z(k),rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k),P(i,j,k),E(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    enddo
  end do
end do

 close(15)

!===================== INITIAL TIME DATA ===============================================

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_dum, rho_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

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
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_ux_dum, rho_ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_uy_dum, rho_uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_uz_dum, rho_uz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, P_dum, P_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, E_dum, E_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Bx_dum, Bx_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, By_dum, By_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, Bz_dum, Bz_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

! Evaluate Initial Vorticity Spectra.
do i = 1, Nx/2+1
  do j = 1, Ny/2
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)


! Evaluate Initial Energy Spectra.
!do i = 1, Nx/2+1
  !do j = 1, Ny
    !do k = 1, Nz
      !write(30,*) i-1,j-1,k-1,abs(nk(i,j,k)),abs(ux_k(i,j,k)),abs(uy_k(i,j,k)),abs(uz_k(i,j,k))
    !enddo
  !end do
!end do

!Energy0 = 0.0d0; y_Energy0 = 0.0d0

do i = 1, Nx
  do j = 1, Ny
    do k = 1, Nz
      ! FFTW Normalisation.
       omega_x(i,j,k) = omega_x(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
       omega_y(i,j,k) = omega_y(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
       omega_z(i,j,k) = omega_z(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Initial Kinetic Energy.
      !Energy0 = Energy0 + (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)
      !y_Energy0 = y_Energy0 + (uy(i,j,k)**2)
      ! Store Grid Data Reproduced from FFTW.
      !write(25,*) x(i),y(j),z(k),omega_x(i,j,k),omega_y(i,j,k),omega_z(i,j,k), &
      !             rho(i,j,k),ux(i,j,k),uy(i,j,k),uz(i,j,k), &
      !             P(i,j,k),E(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    end do
  end do
enddo

! Write the Initial Energy in File.
!write(40,*) 0.0d0, Energy0,y_Energy0

!======================= MAIN PROGRAM =====================================================
t1 = omp_get_wtime()

! Time Loop Starts.
do time = time_min,time_max,dt

t = nint(time/dt) - dint(time_min/dt)

! Online Check of the Progress.
!if (mod(t,int(1.0/dt)) == 0) then
  !print*, t/int(1.0/dt)
!endif

!======================================================================================

! Non-Linear Term Evaluator.
  call derive (Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
               d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
               d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
               d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)

! Time Solvers.
! Adams-Bashforth
  call ab (Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
           rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new, &
           d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
           d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
           d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)

!======================================================================================

!$OMP PARALLEL SHARED(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old,d_rho_uz_k_dt_old),&
!$OMP & SHARED(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old,d_Bz_k_dt_old),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new,d_rho_uz_k_dt_new),&
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new,d_Bz_k_dt_new),&
!$OMP & SHARED(rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k),&
!$OMP & SHARED(rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new),&
!$OMP & SHARED(rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,rho_uz_k_dum,E_k_dum,Bx_k_dum,By_k_dum,Bz_k_dum) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
    do k = 1,Nz
      ! Set the Variables in Proper Format for Next Time Iteration.
      d_rho_k_dt_old(i,j,k) = d_rho_k_dt_new(i,j,k)

      d_rho_ux_k_dt_old(i,j,k) = d_rho_ux_k_dt_new(i,j,k)
      d_rho_uy_k_dt_old(i,j,k) = d_rho_uy_k_dt_new(i,j,k)
      d_rho_uz_k_dt_old(i,j,k) = d_rho_uz_k_dt_new(i,j,k)

      d_E_k_dt_old(i,j,k) = d_E_k_dt_new(i,j,k)

      d_Bx_k_dt_old(i,j,k) = d_Bx_k_dt_new(i,j,k)
      d_By_k_dt_old(i,j,k) = d_By_k_dt_new(i,j,k)
      d_Bz_k_dt_old(i,j,k) = d_Bz_k_dt_new(i,j,k)

      rho_k(i,j,k) = rho_k_new(i,j,k)

      rho_ux_k(i,j,k) = rho_ux_k_new(i,j,k)
      rho_uy_k(i,j,k) = rho_uy_k_new(i,j,k)
      rho_uz_k(i,j,k) = rho_uz_k_new(i,j,k)

      E_k(i,j,k) = E_k_new(i,j,k)

      Bx_k(i,j,k) = Bx_k_new(i,j,k)
      By_k(i,j,k) = By_k_new(i,j,k)
      Bz_k(i,j,k) = Bz_k_new(i,j,k)

      ! Keep Backup of the Arrays for FFTW.
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

! Evaluate Divergence of B in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,Lz,Bx_k,By_k,Bz_k,div_B_k,jx_k,jy_k,jz_k,Ax_k,Ay_k,Az_k) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
      jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
      jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
        if (i == 1 .and. j == 1 .and. k == 1) then
          Ax_k(i,j,k) = jx_k(i,j,k)
          Ay_k(i,j,k) = jy_k(i,j,k)
          Az_k(i,j,k) = jz_k(i,j,k)
        else
          Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
          Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
          Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
        endif
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
      jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
      jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
      Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
      Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
    enddo
  enddo
  do j = Ny/2+1,Ny
    do k = 1,Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
      jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
      jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
      Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
      Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
    enddo
    do k = Nz/2+1,Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      div_B_k(i,j,k) = (0.0d0,1.0d0)*kx*Bx_k(i,j,k) + (0.0d0,1.0d0)*ky*By_k(i,j,k) + (0.0d0,1.0d0)*kz*Bz_k(i,j,k)
      jx_k(i,j,k) = (0.0d0,1.0d0)*ky*Bz_k(i,j,k) - (0.0d0,1.0d0)*kz*By_k(i,j,k)
      jy_k(i,j,k) = (0.0d0,1.0d0)*kz*Bx_k(i,j,k) - (0.0d0,1.0d0)*kx*Bz_k(i,j,k)
      jz_k(i,j,k) = (0.0d0,1.0d0)*kx*By_k(i,j,k) - (0.0d0,1.0d0)*ky*Bx_k(i,j,k)
      Ax_k(i,j,k) = jx_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
      Ay_k(i,j,k) = jy_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
      Az_k(i,j,k) = jz_k(i,j,k)/(kx*kx+ky*ky+kz*kz)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jx_k, jx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jy_k, jy, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, jz_k, jz, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Ax_k, Ax, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Ay_k, Ay, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, Az_k, Az, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

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

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, div_B_k, div_B, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

!$OMP PARALLEL SHARED(spheat,rho,rho_ux,rho_uy,rho_uz,E,Bx,By,Bz,div_B), &
!$OMP & SHARED(ux,uy,uz,rho_dum,ux_dum,uy_dum,uz_dum,P) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      ! FFTW Normalisation.
      jx(i,j,k) = jx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      jy(i,j,k) = jy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      jz(i,j,k) = jz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

      Ax(i,j,k) = Ax(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Ay(i,j,k) = Ay(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Az(i,j,k) = Az(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

      rho(i,j,k) = rho(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

        !if (rho(i,j,k) .lt. 0.0d0) then
        !  rho(i,j,k) = abs(rho(i,j,k))
        !endif

      rho_ux(i,j,k) = rho_ux(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      rho_uy(i,j,k) = rho_uy(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      rho_uz(i,j,k) = rho_uz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

      E(i,j,k) = E(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

      Bx(i,j,k) = Bx(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      By(i,j,k) = By(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Bz(i,j,k) = Bz(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

      ! Evaluate Divergence of B.
      div_B(i,j,k) = div_B(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))

      ! Evaluate Velocity in Real Space.
      ux(i,j,k) = rho_ux(i,j,k)/rho(i,j,k)
      uy(i,j,k) = rho_uy(i,j,k)/rho(i,j,k)
      uz(i,j,k) = rho_uz(i,j,k)/rho(i,j,k)

      ! Evaluate Square of Magnetic Field.
      B2(i,j,k) = Bx(i,j,k)*Bx(i,j,k) + By(i,j,k)*By(i,j,k) + Bz(i,j,k)*Bz(i,j,k)

      ! Keep Backup of the Arrays for FFTW.
      rho_dum(i,j,k) = rho(i,j,k)
      ux_dum(i,j,k) = ux(i,j,k)
      uy_dum(i,j,k) = uy(i,j,k)
      uz_dum(i,j,k) = uz(i,j,k)

      ! Evaluate Pressure
      P(i,j,k) = CS*CS*rho(i,j,k)!( spheat - 1.0d0 ) * ( E(i,j,k) - 0.50d0 * &
                 !( rho_ux(i,j,k)*ux(i,j,k)+rho_uy(i,j,k)*uy(i,j,k)+rho_uz(i,j,k)*uz(i,j,k) - B2(i,j,k) ) )
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

! Evaluate of Structure Functions.
!if (t >= int((time_max-time_max/10.0d0)/dt)) then
!if (t == int(time_max/dt)) then

!!$OMP PARALLEL SHARED(t,x,y,z,delta,P_max,dr), &
!!$OMP & SHARED(rho,ux,uy,uz,Bx,By,Bz) PRIVATE(i,j,k,l,n,o,P_s,q,d_rho,d_ux,d_uy,d_uz,d_Bx,d_By,d_Bz)
!!$OMP DO REDUCTION (+:S_rho,S_u,S_B)

!  do i = 1,Nx,1
!    do j = 1,Ny,1
!      do k = 1,Nz,1
!        do l = 1,Nx,1
!          do n = 1,Ny,1
!            do o = 1,Nz,1
!              dr = dsqrt((x(l)-x(i))*(x(l)-x(i))+(y(n)-y(j))*(y(n)-y(j))+(z(o)-z(k))*(z(o)-z(k)))
!              q = int(dr/delta) + 1
!              d_rho(q) = rho(l,n,o) - rho(i,j,k)
!              d_ux(q) = ux(l,n,o) - ux(i,j,k)
!              d_uy(q) = uy(l,n,o) - uy(i,j,k)
!              d_uz(q) = uz(l,n,o) - uz(i,j,k)
!              d_Bx(q) = Bx(l,n,o) - Bx(i,j,k)
!              d_By(q) = By(l,n,o) - By(i,j,k)
!              d_Bz(q) = Bz(l,n,o) - Bz(i,j,k)
!                do P_s = 1,P_max,1
!                  S_rho(P_s,q) = S_rho(P_s,q) + abs(d_rho(q))**P_s
!                  s_u(P_s,q) = S_u(P_s,q) + abs(dsqrt(d_ux(q)*d_ux(q)+d_uy(q)*d_uy(q)+d_uz(q)*d_uz(q)))**P_s
!                  s_B(P_s,q) = S_B(P_s,q) + abs(dsqrt(d_Bx(q)*d_Bx(q)+d_By(q)*d_By(q)+d_Bz(q)*d_Bz(q)))**P_s
!                enddo !P_s
!            enddo !o
!          enddo !n
!        enddo !l
!      enddo !k
!    enddo !j
!  enddo !i

!!$OMP END PARALLEL

!endif

  !call dfftw_init_threads(iret)
  !call dfftw_plan_with_nthreads(thread_num)
  !call dfftw_plan_dft_r2c_3d_ (plan_forward, Nx, Ny, Nz, rho_dum, rho_k, FFTW_ESTIMATE)
  !call dfftw_execute_ (plan_forward)
  !call dfftw_destroy_plan_ (plan_forward)

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

! Evaluate Vorticity in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,Lz,ux_k,uy_k,uz_k,omega_x_k,omega_y_k,omega_z_k),&
!$OMP & SHARED(omega_x_k_dum,omega_y_k_dum,omega_z_k_dum) PRIVATE(i,j,k,kx,ky,kz)
!$OMP DO

do i = 1, Nx/2+1
  do j = 1, Ny/2
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat(j-1)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
  do j = Ny/2+1, Ny
    do k = 1, Nz/2
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat(k-1)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
    do k = Nz/2+1, Nz
      kx = 2.0d0*pi*dfloat(i-1)/Lx
      ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
      kz = 2.0d0*pi*dfloat((k-1)-Nz)/Lz
      omega_x_k(i,j,k) = (0.0d0,1.0d0)*ky*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*uy_k(i,j,k)
      omega_y_k(i,j,k) = (0.0d0,1.0d0)*kx*uz_k(i,j,k) - (0.0d0,1.0d0)*kz*ux_k(i,j,k)
      omega_z_k(i,j,k) = (0.0d0,1.0d0)*kx*uy_k(i,j,k) - (0.0d0,1.0d0)*ky*ux_k(i,j,k)

      omega_x_k_dum(i,j,k) = omega_x_k(i,j,k)
      omega_y_k_dum(i,j,k) = omega_y_k(i,j,k)
      omega_z_k_dum(i,j,k) = omega_z_k(i,j,k)
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_x_k_dum, omega_x, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_y_k_dum, omega_y, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_3d_ (plan_backward, Nx, Ny, Nz, omega_z_k_dum, omega_z, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)


! FFTW Normalisation and omega^2 Evaluation.

!$OMP PARALLEL SHARED(omega_x,omega_y,omega_z,omega2) PRIVATE(i,j,k)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      omega_x(i,j,k) = omega_x(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega_y(i,j,k) = omega_y(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega_z(i,j,k) = omega_z(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      omega2(i,j,k) = omega_x(i,j,k)**2 + omega_y(i,j,k)**2 + omega_z(i,j,k)**2
    enddo
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

! Evaluate Energy Spectra at a Specified Time.
!if (t >= int((time_max-time_max/10.0d0)/dt)) then

!!! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file. !!!
! !!$OMP PARALLEL SHARED(ux_k,uy_k,uz_k) PRIVATE(i,j,k)
! !!$OMP DO REDUCTION (+:Ek)

!  do i = 1,Nx/2+1
!    do j = 1,Ny
!      do k = 1,Nz
!        Ek(i,j,k) = Ek(i,j,k) + sqrt(abs(ux_k(i,j,k))**2 + abs(uy_k(i,j,k))**2 + abs(uz_k(i,j,k))**2)&
!                    /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!        Bk(i,j,k) = Bk(i,j,k) + sqrt(abs(Bx_k(i,j,k))**2 + abs(By_k(i,j,k))**2 + abs(Bz_k(i,j,k))**2)&
!                   /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
!      enddo
!    enddo
!  enddo

! !!$OMP END PARALLEL

!endif

Pressure = 0.0d0; Energy = 0.0d0; T_Energy = 0.0d0; B_Field = 0.0d0
 C0 = 0.0d0; C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0; div_B_Tot = 0.0d0
 Hf = 0.0d0; Hm = 0.0d0; HT = 0.0d0

! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file.
! !$OMP PARALLEL SHARED(mu_0,t,dt,x,y,z,omega_x,omega_y,omega_z,omega2,rho,ux,uy,uz) PRIVATE(i,j,k)
! !$OMP DO REDUCTION (+:Pressure,Energy,y_Energy,B_Field,C0,C1,C2,C3,div_B_Tot,Hf,Hm,HT)

do i = 1,Nx
  do j = 1,Ny
    do k = 1,Nz
      ! Evaluate Pressure
      Pressure = Pressure + P(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Energy
      !Energy = Energy + E(i,j,k)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Energy = Energy + (ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Growth Rate.
      T_Energy = T_Energy + E(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Magnetic Field.
      B_Field = B_Field + (Bx(i,j,k)**2+By(i,j,k)**2+Bz(i,j,k)**2)/(2.0d0*mu_0*dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Evaluate Casimirs.
      C0 = C0 + dsqrt(omega2(i,j,k))**0.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C1 = C1 + dsqrt(omega2(i,j,k))**1.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C2 = C2 + dsqrt(omega2(i,j,k))**2.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      C3 = C3 + dsqrt(omega2(i,j,k))**3.0d0/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Hf = Hf + ( ux(i,j,k)*omega_x(i,j,k) + uy(i,j,k)*omega_y(i,j,k) + uz(i,j,k)*omega_z(i,j,k) ) &
                /(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      Hm = Hm + ( Ax(i,j,k)*Bx(i,j,k) + Ay(i,j,k)*By(i,j,k) + Az(i,j,k)*Bz(i,j,k) )/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      HT = HT + ( (ux(i,j,k)+Ax(i,j,k))*(omega_x(i,j,k)+Bx(i,j,k)) + (uy(i,j,k)+Ay(i,j,k))*(omega_y(i,j,k)+By(i,j,k)) &
                + (uz(i,j,k)+Az(i,j,k))*(omega_z(i,j,k)+Bz(i,j,k)) )/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Check for Div B = 0
      div_B_Tot = div_B_Tot + div_B(i,j,k)/(dfloat(Nx)*dfloat(Ny)*dfloat(Nz))
      ! Write Grid Data in State Files.
        !if (mod(float(t),100.0) == 0.0) then
       ! if (mod(float(t),10000.0) == 0.0) then
        !if (t == int(time_max/dt) .and. k == Nz/2) then
        !write(filename, '("output/fort.",I8.8)') t+100
        !open(unit=t+100,file=filename,status='unknown')
        !write(t+100,*) x(i),y(j),z(k),ux(i,j,k),uy(i,j,k),uz(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
        !endif
    enddo
  enddo
enddo

! !$OMP END PARALLEL

write(40,*) time,Pressure,Energy,T_Energy,B_Field,C1,C2,Hf,Hm,HT,div_B_Tot
  call flush(40)

if (mod(float(t),10000.0) == 0.0) then
  close(t+100)
endif

enddo ! time

t2 = omp_get_wtime()

write(5,*) "Time taken for the run =",(t2 - t1)/(60.0d0*60.0d0),"Hours"

 close(5)
 close(35)
 close(40)

contains

!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive(Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
                  d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
                  d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
                  d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nz,Nh,iret,thread_num
real ( kind = 8 ) pi,time,dt,Lx,Ly,Lz,dx,dy,dz,kx,ky,kz,mu,ms,CS,mu_0,eta,spheat,kf,A,B,C,time_max

real ( kind = 8 ) rho(Nx,Ny,Nz),ux(Nx,Ny,Nz),uy(Nx,Ny,Nz),uz(Nx,Ny,Nz),E(Nx,Ny,Nz),P(Nx,Ny,Nz)
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

end subroutine derive

!===================================================================================
!=================== SUBROUTINE ADAMS BASHFORTH ====================================
!===================================================================================

subroutine ab(Nx,Ny,Nz,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,rho_uz_k,E_k,Bx_k,By_k,Bz_k, &
              rho_k_new,rho_ux_k_new,rho_uy_k_new,rho_uz_k_new,E_k_new,Bx_k_new,By_k_new,Bz_k_new, &
              d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
              d_rho_uz_k_dt_old,d_rho_uz_k_dt_new,d_E_k_dt_old,d_E_k_dt_new, &
              d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new,d_Bz_k_dt_old,d_Bz_k_dt_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nz,Nh,iret,thread_num
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

end subroutine ab

!====================================================================================

end subroutine tara3dMHDOMP
