subroutine tara2dMHDOMP(arg)

use omp_lib
use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3.f"

integer (kind=4) Np,Nx,Ny,Nh

integer ( kind = 4 ) i,j,t
real ( kind = 8 ) Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,G1,G2,G3
real ( kind = 8 ) spheat,rho0,ms,U0,mu,k,kf,A,sigma,uy0,M,CS,Re,P0,eta,mu_0,MA,VA,B0,ba
real ( kind = 8 ) Pressure,Energy,Energy0,y_Energy,y_Energy0,B_Field,C0,C1,C2,C3,div_B_tot

real (kind = 8), dimension (:), allocatable :: x,y
real (kind=8), dimension(:, :), allocatable :: rho,ux,uy,omega,j_B,rho_ux,rho_uy
real (kind=8), dimension(:, :), allocatable :: P,E,Bx,By,B2,div_B
real (kind=8), dimension(:, :), allocatable :: rho_dum,ux_dum,uy_dum,rho_ux_dum,rho_uy_dum
real (kind=8), dimension(:, :), allocatable :: P_dum,E_dum,Bx_dum,By_dum

complex (kind=8), dimension(:, :), allocatable :: rho_k,ux_k,uy_k,omega_k,j_k,rho_ux_k,rho_uy_k
complex (kind=8), dimension(:, :), allocatable :: P_k,E_k,Bx_k,By_k,div_B_k,Ek
complex (kind=8), dimension(:, :), allocatable :: rho_k_dum,ux_k_dum,uy_k_dum,omega_k_dum,j_k_dum
complex (kind=8), dimension(:, :), allocatable :: rho_ux_k_dum,rho_uy_k_dum,E_k_dum,Bx_k_dum,By_k_dum
complex (kind=8), dimension(:, :), allocatable :: rho_k_new,rho_ux_k_new,rho_uy_k_new,E_k_new,Bx_k_new,By_k_new

complex (kind=8), dimension(:, :), allocatable :: d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old
complex (kind=8), dimension(:, :), allocatable :: d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old
complex (kind=8), dimension(:, :), allocatable :: d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new
complex (kind=8), dimension(:, :), allocatable :: d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new

integer ( kind = 8 ) plan_forward,plan_backward ! FFTW

integer ( kind = 4 ) thread_num,num_thread,proc_num ! OMP
real ( kind = 8 ) t1,t2

common/comm/Lx,Ly,spheat,mu,ms,CS,mu_0,eta,dt,kf,A

integer,parameter :: seed = 99999999

character (len=90) :: filename
character (len=32) :: arg

cfg = parse_cfg(arg)

call cfg%get("grid","Nx",Nx)
call cfg%get("grid","Ny",Ny)

Nh = (Nx/2) + 1

call srand(seed)

allocate(x(Nx))
allocate(y(Ny))

allocate(rho(Nx,Ny))
allocate(ux(Nx,Ny))
allocate(uy(Nx,Ny))
allocate(omega(Nx,Ny))
allocate(j_B(Nx,Ny))
allocate(rho_ux(Nx,Ny))
allocate(rho_uy(Nx,Ny))
allocate(P(Nx,Ny))
allocate(E(Nx,Ny))
allocate(Bx(Nx,Ny))
allocate(By(Nx,Ny))
allocate(B2(Nx,Ny))
allocate(div_B(Nx,Ny))
allocate(rho_dum(Nx,Ny))
allocate(ux_dum(Nx,Ny))
allocate(uy_dum(Nx,Ny))
allocate(rho_ux_dum(Nx,Ny))
allocate(rho_uy_dum(Nx,Ny))
allocate(P_dum(Nx,Ny))
allocate(E_dum(Nx,Ny))
allocate(Bx_dum(Nx,Ny))
allocate(By_dum(Nx,Ny))

allocate(rho_k(Nh,Ny))
allocate(ux_k(Nh,Ny))
allocate(uy_k(Nh,Ny))
allocate(omega_k(Nh,Ny))
allocate(j_k(Nh,Ny))
allocate(rho_ux_k(Nh,Ny))
allocate(rho_uy_k(Nh,Ny))
allocate(P_k(Nh,Ny))
allocate(E_k(Nh,Ny))
allocate(Bx_k(Nh,Ny))
allocate(By_k(Nh,Ny))
allocate(div_B_k(Nh,Ny))
allocate(Ek(Nh,Ny))
allocate(rho_k_dum(Nh,Ny))
allocate(ux_k_dum(Nh,Ny))
allocate(uy_k_dum(Nh,Ny))
allocate(omega_k_dum(Nh,Ny))
allocate(j_k_dum(Nh,Ny))
allocate(rho_ux_k_dum(Nh,Ny))
allocate(rho_uy_k_dum(Nh,Ny))
allocate(E_k_dum(Nh,Ny))
allocate(Bx_k_dum(Nh,Ny))
allocate(By_k_dum(Nh,Ny))
allocate(rho_k_new(Nh,Ny))
allocate(rho_ux_k_new(Nh,Ny))
allocate(rho_uy_k_new(Nh,Ny))
allocate(E_k_new(Nh,Ny))
allocate(Bx_k_new(Nh,Ny))
allocate(By_k_new(Nh,Ny))
allocate(d_rho_k_dt_old(Nh,Ny))
allocate(d_rho_ux_k_dt_old(Nh,Ny))
allocate(d_rho_uy_k_dt_old(Nh,Ny))
allocate(d_E_k_dt_old(Nh,Ny))
allocate(d_Bx_k_dt_old(Nh,Ny))
allocate(d_By_k_dt_old(Nh,Ny))
allocate(d_rho_k_dt_new(Nh,Ny))
allocate(d_rho_ux_k_dt_new(Nh,Ny))
allocate(d_rho_uy_k_dt_new(Nh,Ny))
allocate(d_E_k_dt_new(Nh,Ny))
allocate(d_Bx_k_dt_new(Nh,Ny))
allocate(d_By_k_dt_new(Nh,Ny))

!===================== FILENAMES ==============================================

open(unit=5,file='output/System_information.dat',status='unknown')
open(unit=15,file='output/Initial_Grid_Data.dat',status='unknown')
!open(unit=20,file='INPUT.dat',status='old')  ! This input file is the file generated from Vorticity code.
!open(unit=25,file='Initial_Grid_Data_Reproduced.dat',status='unknown')
!open(unit=30,file='Initial_Energy_Spectra.dat',status='unknown')
!open(unit=35,file='Energy_Spectra.dat',status='unknown')
open(unit=40,file='output/Energy.dat',status='unknown')

!===================== USER INPUTS ============================================

! Define Number of Threads.
proc_num = omp_get_num_procs()
thread_num = 40
call omp_set_num_threads (thread_num)

write(5,*) "Number of Processors=", proc_num, "Number of Threads=", thread_num

! System Size.
call cfg%get("length","Lx",Lx)
call cfg%get("length","Ly",Ly)

Lx = Lx*pi; Ly = Ly*pi

! Grid Resolution.
dx = Lx/dfloat(Nx); dy = Ly/dfloat(Ny)

! Runtime Details and Time Step Width.
call cfg%get("time","time_min",time_min)
call cfg%get("time","time_max",time_max)
call cfg%get("time","dt",dt)

! Ratio of Specific Heats/ Adiabetic Index/ Heat Capacity Ratio/ Poisson Constant.
spheat = 5.0d0/3.0d0

! Co-efficient of Viscosity.
call cfg%get("resistivity","nu",mu)

! Co-efficient of Resistivity.
eta = 0.000100d0 !0.00010d0; 0.000030d0; 0.000010d0

! Magnetic Permiability.
mu_0 = 1.0d0

! Mass of the Fluid Element.
ms = 1.0d0

! Background Density.
rho0 = 1.0d0

! Initial Pressure.
P0 = 1.0d0

! Mach Number.
!M = 0.050d0
M = 0.010d0

! Alfven Mach Number.
MA = 1.0d0

! Maximum Velocity.
!U0 = 0.14732247502913884d0
!U0 = 0.6450d0

! Sound Speed.
 CS = 1.0d0/M!dsqrt(spheat*P0/rho0)
 U0 = M * CS

! Sound Speed.
! CS = U0/M

! Alfven Speed.
VA = U0/MA

! Initial Magnetic Field.
B0 = VA*dsqrt(mu_0*rho0)

write(5,*) "Sound Speed=", CS, "Initial Velocity=", U0
write(5,*) "Alfven Speed=", VA, "Initial Magnetic Field=", B0

! Forcing Scale.
kf = 1.0d0

! Mode Number of Perturbation.
!k = 2.0d0*pi

! Amplitude of Forcing
A = 1.0d0

! Magnetic Shear Width.
!ba = 0.010d0

! Perturbation Width.
!sigma = 4.0d0*a

! Perturbation Amplitude.
!uy0 = 0.00010d0

! Raynold's Number.
Re = ms*rho0*U0*Lx/mu; write(5,*) "Re =", Re

! Grid Generation.
do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0d0+real(j-1)*dy
    ! Initial Density Distribution.
    rho(i,j) = rho0
    ! Initial Velocity Distribution.
    ux(i,j) = -dsin(kf*y(j))!U0 * (dtanh((y(j)-Ly/3.0d0)/a) - dtanh((y(j)-2.0d0*Ly/3.0d0)/a) - 1.0d0)
    !ux(i,j) = U0*tanh((y(j)-Ly/2.0d0)/a)
    uy(i,j) = dsin(kf*x(i))!0.0d0
    ! Initial Velocity Perturbation.
    !uy(i,j) = uy(i,j) + uy0*dsin(k*x(i))*dexp(-((y(j)-Ly/3.0d0)**2)/(sigma**2))
    !uy(i,j) = uy(i,j) + uy0*dsin(k*x(i))*dexp(-((y(j)-2.0d0*Ly/3.0d0)**2)/(sigma**2))
    !uy(i,j) = uy(i,j) + uy0*dsin(k*x(i))*exp(-((y(j)-Ly/2.0d0)**2)/(sigma**2))
    ! Initial Velocity Distribution.
    !read(20,*) G1,G2,G3,ux(i,j),uy(i,j)  ! You may define your own function here.
    ! Initial Pressure Distribution.
    P(i,j) = CS*CS*rho(i,j)!P0 !* (1.0d0 + 0.10d0 * (dtanh((y(j)-Ly/3.0d0)/ba) - dtanh((y(j)-2.0d0*Ly/3.0d0)/ba) - 1.0d0) )
    ! Initial Magnetic Field Profile.
    Bx(i,j) = B0!-dsin(y(j))!B0 * (dtanh((y(j)-Ly/3.0d0)/ba) - dtanh((y(j)-2.0d0*Ly/3.0d0)/ba) - 1.0d0)
    By(i,j) = B0!dsin(2.0d0*x(i))!0.0d0
    ! Initial Energy Distribution.
    E(i,j) = P(i,j)/(spheat -1.0d0) + 0.50d0*rho(i,j) * (ux(i,j)**2 + uy(i,j)**2) + 0.50d0*(Bx(i,j)**2 + By(i,j)**2)
    ! Initial Combined Variables.
    rho_ux(i,j) = rho(i,j) * ux(i,j)
    rho_uy(i,j) = rho(i,j) * uy(i,j)
    ! Keep Backup of the Arrays for FFTW.
    rho_dum(i,j) = rho(i,j)
    ux_dum(i,j) = ux(i,j)
    uy_dum(i,j) = uy(i,j)
    rho_ux_dum(i,j) = rho_ux(i,j)
    rho_uy_dum(i,j) = rho_uy(i,j)
    P_dum(i,j) = P(i,j)
    E_dum(i,j) = E(i,j)
    Bx_dum(i,j) = Bx(i,j)
    By_dum(i,j) = By(i,j)
    ! Write Initial Density and Velocity Distribution in File.
    write(15,*) x(i),y(j),rho(i,j),ux(i,j),uy(i,j),P(i,j),E(i,j),Bx(i,j),By(i,j)
  end do
end do

 close(15)

!===================== INITIAL TIME DATA ===============================================

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, rho_dum, rho_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, rho_ux_dum, rho_ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, rho_uy_dum, rho_uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, P_dum, P_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, E_dum, E_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, Bx_dum, Bx_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, By_dum, By_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

! Evaluate Initial Vorticity Spectra.
do i = 1, Nx/2+1
  do j = 1, Ny/2
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  ky = 2.0d0*pi*dfloat(j-1)/Ly
  omega_k(i,j) = (0.0d0,1.0d0)*kx*uy_k(i,j) - (0.0d0,1.0d0)*ky*ux_k(i,j)
  omega_k_dum(i,j) = omega_k(i,j)
  enddo
  do j = Ny/2+1,Ny
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
  omega_k(i,j) = (0.0d0,1.0d0)*kx*uy_k(i,j) - (0.0d0,1.0d0)*ky*ux_k(i,j)
  omega_k_dum(i,j) = omega_k(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omega_k_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

! Evaluate Initial Energy Spectra.
!do i = 1, Nx/2+1
  !do j = 1, Ny
  !write(30,*) i-1,j-1,abs(nk(i,j)),abs(ukx(i,j)),abs(uky(i,j))
  !end do
!end do

!Energy0 = 0.0d0; y_Energy0 = 0.0d0

do i = 1, Nx
  do j = 1, Ny
  ! FFTW Normalisation.
   omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
  ! Evaluate Initial Kinetic Energy.
  !Energy0 = Energy0 + (ux(i,j)**2 + uy(i,j)**2)
  !y_Energy0 = y_Energy0 + (uy(i,j)**2)
  ! Store Grid Data Reproduced from FFTW.
  !write(25,*) x(i),y(j),omega(i,j),n(i,j),ux(i,j),uy(i,j)
  end do
end do

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
  call derive (Nx,Ny,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k, &
               d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
               d_E_k_dt_old,d_E_k_dt_new,d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new)

! Time Solvers.
! Adams-Bashforth
  call ab (Nx,Ny,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k, &
           rho_k_new,rho_ux_k_new,rho_uy_k_new,E_k_new,Bx_k_new,By_k_new, &
           d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
           d_E_k_dt_old,d_E_k_dt_new,d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new)

!======================================================================================

!$OMP PARALLEL SHARED(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old),&
!$OMP & SHARED(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new),&
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new),&
!$OMP & SHARED(rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k),&
!$OMP & SHARED(rho_k_new,rho_ux_k_new,rho_uy_k_new,E_k_new,Bx_k_new,By_k_new),&
!$OMP & SHARED(rho_k_dum,rho_ux_k_dum,rho_uy_k_dum,E_k_dum,Bx_k_dum,By_k_dum) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
    ! Set the Variables in Proper Format for Next Time Iteration.
    d_rho_k_dt_old(i,j) = d_rho_k_dt_new(i,j)

    d_rho_ux_k_dt_old(i,j) = d_rho_ux_k_dt_new(i,j)
    d_rho_uy_k_dt_old(i,j) = d_rho_uy_k_dt_new(i,j)

    d_E_k_dt_old(i,j) = d_E_k_dt_new(i,j)

    d_Bx_k_dt_old(i,j) = d_Bx_k_dt_new(i,j)
    d_By_k_dt_old(i,j) = d_By_k_dt_new(i,j)

    rho_k(i,j) = rho_k_new(i,j)

    rho_ux_k(i,j) = rho_ux_k_new(i,j)
    rho_uy_k(i,j) = rho_uy_k_new(i,j)

    E_k(i,j) = E_k_new(i,j)

    Bx_k(i,j) = Bx_k_new(i,j)
    By_k(i,j) = By_k_new(i,j)

    ! Keep Backup of the Arrays for FFTW.
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

! Evaluate Divergence of B in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,Bx_k,By_k,div_B_k) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny/2
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  ky = 2.0d0*pi*dfloat(j-1)/Ly
  div_B_k(i,j) = (0.0d0,1.0d0)*kx*Bx_k(i,j) + (0.0d0,1.0d0)*ky*By_k(i,j)
  enddo
  do j = Ny/2+1,Ny
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
  div_B_k(i,j) = (0.0d0,1.0d0)*kx*Bx_k(i,j) + (0.0d0,1.0d0)*ky*By_k(i,j)
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

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, div_B_k, div_B, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

!$OMP PARALLEL SHARED(spheat,rho,rho_ux,rho_uy,E,Bx,By,div_B,ux,uy,ux_dum,uy_dum,P) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
  ! FFTW Normalisation.
  rho(i,j) = rho(i,j)/(dfloat(Nx)*dfloat(Ny))

  rho_ux(i,j) = rho_ux(i,j)/(dfloat(Nx)*dfloat(Ny))
  rho_uy(i,j) = rho_uy(i,j)/(dfloat(Nx)*dfloat(Ny))

  E(i,j) = E(i,j)/(dfloat(Nx)*dfloat(Ny))

  Bx(i,j) = Bx(i,j)/(dfloat(Nx)*dfloat(Ny))
  By(i,j) = By(i,j)/(dfloat(Nx)*dfloat(Ny))

  ! Evaluate Divergence of B.
  div_B(i,j) = div_B(i,j)/(dfloat(Nx)*dfloat(Ny))

  ! Evaluate Velocity in Real Space.
  ux(i,j) = rho_ux(i,j)/rho(i,j)
  uy(i,j) = rho_uy(i,j)/rho(i,j)

  ! Evaluate Square of Magnetic Field.
  B2(i,j) = Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j)

  ! Keep Backup of the Arrays for FFTW.
  ux_dum(i,j) = ux(i,j)
  uy_dum(i,j) = uy(i,j)

  ! Evaluate Pressure
  P(i,j) = CS*CS*rho(i,j)!( spheat - 1.0d0 ) * ( E(i,j) - 0.50d0 * ( rho_ux(i,j)*ux(i,j)+rho_uy(i,j)*uy(i,j) - B2(i,j) ) )
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_dum, ux_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_dum, uy_k, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

! Evaluate Vorticity in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,ux_k,uy_k,omega_k,omega_k_dum,j_k,j_k_dum) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1, Nx/2+1
  do j = 1, Ny/2
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  ky = 2.0d0*pi*dfloat(j-1)/Ly
  omega_k(i,j) = (0.0d0,1.0d0)*kx*uy_k(i,j) - (0.0d0,1.0d0)*ky*ux_k(i,j)
  j_k(i,j) = (0.0d0,1.0d0)*kx*By_k(i,j) - (0.0d0,1.0d0)*ky*Bx_k(i,j)
  omega_k_dum(i,j) = omega_k(i,j)
  j_k_dum(i,j) = j_k(i,j)
  enddo
  do j = Ny/2+1,Ny
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
  omega_k(i,j) = (0.0d0,1.0d0)*kx*uy_k(i,j) - (0.0d0,1.0d0)*ky*ux_k(i,j)
  j_k(i,j) = (0.0d0,1.0d0)*kx*By_k(i,j) - (0.0d0,1.0d0)*ky*Bx_k(i,j)
  omega_k_dum(i,j) = omega_k(i,j)
  j_k_dum(i,j) = j_k(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omega_k_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, j_k_dum, j_B, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_destroy_plan_ (plan_forward)
  call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisation

!$OMP PARALLEL SHARED(omega,j_B) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
  omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
  j_B(i,j) = j_B(i,j)/(dfloat(Nx)*dfloat(Ny))
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

! Evaluate Energy Spectra at a Specified Time.
!if (t == 10000) then

!!! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file. !!!
! !!$OMP PARALLEL SHARED(ukx,uky) PRIVATE(i,j)
! !!$OMP DO REDUCTION (+:Ek)

  !do i = 1,Nx/2+1
  !  do j = 1,Ny
  !  Ek(i,j) = Ek(i,j) + sqrt(abs(ukx(i,j))**2 + abs(uky(i,j))**2)
  !  write(35,*) i,j,sqrt(float((i-1)**2)+float((j-1)**2)),abs(Ek(i,j))
  !  enddo
  !enddo

! !!$OMP END PARALLEL

!endif

Pressure = 0.0d0; Energy = 0.0d0; y_Energy = 0.0d0; B_Field = 0.0d0
 C0 = 0.0d0; C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0; div_B_Tot = 0.0d0

! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file.
! !$OMP PARALLEL SHARED(t,dt,x,y,omega,j,n,ux,uy) PRIVATE(i,j)
! !$OMP DO REDUCTION (+:Energy,y_Energy)

do i = 1,Nx
  do j = 1,Ny
  ! Evaluate Pressure
  !Pressure = Pressure + P(i,j)/(dfloat(Nx)*dfloat(Ny))
  ! Evaluate Energy
  Energy = Energy + (ux(i,j)**2+uy(i,j)**2)/(dfloat(Nx)*dfloat(Ny))!rho(i,j)*(ux(i,j)**2 + uy(i,j)**2)
  ! Evaluate Growth Rate.
  !y_Energy = y_Energy + rho(i,j)*(uy(i,j)**2)/(2.0d0*dfloat(Nx)*dfloat(Ny))
  ! Evaluate Magnetic Field.
  B_Field = B_Field + (Bx(i,j)**2+By(i,j)**2)/(dfloat(Nx)*dfloat(Ny))
  ! Evaluate Casimirs.
  !C0 = C0 + omega(i,j)**0.0d0/(dfloat(Nx)*dfloat(Ny))
  !C1 = C1 + omega(i,j)**1.0d0/(dfloat(Nx)*dfloat(Ny))
  C2 = C2 + (omega(i,j)**2.0d0)/(dfloat(Nx)*dfloat(Ny))
  C3 = C3 + (j_B(i,j)**2.0d0)/(dfloat(Nx)*dfloat(Ny))
  ! Check for Div B = 0
  div_B_Tot = div_B_Tot + div_B(i,j)/(dfloat(Nx)*dfloat(Ny))
    ! Write Grid Data in State Files.
    if (mod(float(t),100000.0) == 0.0) then
      write(filename, '("output/fort.",I8.8)') t+100000
      open(unit=t+100000,file=filename,status='unknown')
      write(t+100000,*) x(i),y(j),rho(i,j),ux(i,j),uy(i,j),P(i,j),E(i,j),Bx(i,j),By(i,j)
    endif
  enddo
enddo

! !$OMP END PARALLEL

if (mod(float(t),100.0) == 0.0) then
write(40,*) time,Energy,B_Field,C2,C3,div_B_Tot
endif
  call flush (40)

if (mod(float(t),100000.0) == 0.0) then
  close(t+100000)
endif
enddo ! time

t2 = omp_get_wtime()

write(5,*) "Time taken for the run = ",t2 - t1

 close(5)
 close(40)

contains

!===================================================================================
!================== SUBROUTINE NONLINEAR DERIVATIVE ================================
!===================================================================================

subroutine derive(Nx,Ny,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k, &
                  d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
                  d_E_k_dt_old,d_E_k_dt_new,d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh
real ( kind = 8 ) pi,time,dt,Lx,Ly,kx,ky,mu,ms,CS,mu_0,eta,spheat,kf,A

real ( kind = 8 ) rho(Nx,Ny),ux(Nx,Ny),uy(Nx,Ny),E(Nx,Ny),P(Nx,Ny),Bx(Nx,Ny),By(Nx,Ny),B2(Nx,Ny),Source(Nx,Ny)
real ( kind = 8 ) rho_ux(Nx,Ny),rho_uy(Nx,Ny),Mom_x_1(Nx,Ny),Mom_x_2(Nx,Ny),Mom_y_1(Nx,Ny),Mom_y_2(Nx,Ny)

real ( kind = 8 ) d_ux_dx(Nx,Ny),d_uy_dy(Nx,Ny),Fx(Nx,Ny),Fy(Nx,Ny)
real ( kind = 8 ) Energy_x(Nx,Ny),Energy_y(Nx,Ny),E_Visc(Nx,Ny),Mag_x(Nx,Ny),Mag_y(Nx,Ny)
real ( kind = 8 ) d_Bx_dy(Nx,Ny),d_By_dx(Nx,Ny),curl_B(Nx,Ny)

complex ( kind = 8 ) rho_k(Nh,Ny),ux_k(Nh,Ny),uy_k(Nh,Ny),rho_ux_k(Nh,Ny),rho_uy_k(Nh,Ny)
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

end subroutine derive

!===================================================================================
!=================== SUBROUTINE ADAMS BASHFORTH ====================================
!===================================================================================

subroutine ab(Nx,Ny,Nh,pi,time,rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k, &
              rho_k_new,rho_ux_k_new,rho_uy_k_new,E_k_new,Bx_k_new,By_k_new, &
              d_rho_k_dt_old,d_rho_k_dt_new,d_rho_ux_k_dt_old,d_rho_ux_k_dt_new,d_rho_uy_k_dt_old,d_rho_uy_k_dt_new, &
              d_E_k_dt_old,d_E_k_dt_new,d_Bx_k_dt_old,d_Bx_k_dt_new,d_By_k_dt_old,d_By_k_dt_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh
real ( kind = 8 ) pi,time,dt,Lx,Ly,spheat,mu,ms,CS,mu_0,eta,kf,A

complex ( kind = 8 ) rho_k(Nh,Ny),rho_ux_k(Nh,Ny),rho_uy_k(Nh,Ny),E_k(Nh,Ny),Bx_k(Nh,Ny),By_k(Nh,Ny)
complex ( kind = 8 ) rho_k_new(Nh,Ny),rho_ux_k_new(Nh,Ny),rho_uy_k_new(Nh,Ny),E_k_new(Nh,Ny),Bx_k_new(Nh,Ny),By_k_new(Nh,Ny)

complex ( kind = 8 ) d_rho_k_dt_old(Nh,Ny),d_rho_ux_k_dt_old(Nh,Ny),d_rho_uy_k_dt_old(Nh,Ny),d_E_k_dt_old(Nh,Ny)
complex ( kind = 8 ) d_Bx_k_dt_old(Nh,Ny),d_By_k_dt_old(Nh,Ny)
complex ( kind = 8 ) d_rho_k_dt_new(Nh,Ny),d_rho_ux_k_dt_new(Nh,Ny),d_rho_uy_k_dt_new(Nh,Ny),d_E_k_dt_new(Nh,Ny)
complex ( kind = 8 ) d_Bx_k_dt_new(Nh,Ny),d_By_k_dt_new(Nh,Ny)

common/comm/Lx,Ly,spheat,mu,ms,CS,mu_0,eta,dt,kf,A

!$OMP PARALLEL SHARED(dt,rho_k,rho_ux_k,rho_uy_k,E_k,Bx_k,By_k),&
!$OMP & SHARED(rho_k_new,rho_ux_k_new,rho_uy_k_new,E_k_new,Bx_k_new,By_k_new),&
!$OMP & SHARED(d_rho_k_dt_new,d_rho_ux_k_dt_new,d_rho_uy_k_dt_new),&
!$OMP & SHARED(d_E_k_dt_new,d_Bx_k_dt_new,d_By_k_dt_new),&
!$OMP & SHARED(d_rho_k_dt_old,d_rho_ux_k_dt_old,d_rho_uy_k_dt_old),&
!$OMP & SHARED(d_E_k_dt_old,d_Bx_k_dt_old,d_By_k_dt_old) PRIVATE(i,j)
!$OMP DO

do i = 1,Nh
  do j = 1,Ny
    ! Density Equation Evolution.
    rho_k_new(i,j) = rho_k(i,j) + ( (3.0d0/2.0d0)*d_rho_k_dt_new(i,j) - (1.0d0/2.0d0)*d_rho_k_dt_old(i,j) )*dt
    ! Momentum Equation Evolution.
    rho_ux_k_new(i,j) = rho_ux_k(i,j) + ( (3.0d0/2.0d0)*d_rho_ux_k_dt_new(i,j) - (1.0d0/2.0d0)*d_rho_ux_k_dt_old(i,j) )*dt
    rho_uy_k_new(i,j) = rho_uy_k(i,j) + ( (3.0d0/2.0d0)*d_rho_uy_k_dt_new(i,j) - (1.0d0/2.0d0)*d_rho_uy_k_dt_old(i,j) )*dt
    ! Energy Equation Evolution.
    E_k_new(i,j) = E_k(i,j) + ( (3.0d0/2.0d0)*d_E_k_dt_new(i,j) - (1.0d0/2.0d0)*d_E_k_dt_old(i,j) )*dt
    ! Energy Equation Evolution.
    Bx_k_new(i,j) = Bx_k(i,j) + ( (3.0d0/2.0d0)*d_Bx_k_dt_new(i,j) - (1.0d0/2.0d0)*d_Bx_k_dt_old(i,j) )*dt
    By_k_new(i,j) = By_k(i,j) + ( (3.0d0/2.0d0)*d_By_k_dt_new(i,j) - (1.0d0/2.0d0)*d_By_k_dt_old(i,j) )*dt
  end do
end do

!$OMP END DO
!$OMP END PARALLEL

return

end subroutine ab

!====================================================================================


end subroutine tara2dMHDOMP
