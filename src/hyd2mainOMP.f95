subroutine tara2dHydOMP(arg)

use omp_lib
use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

real ( kind = 8 ), parameter :: pi=3.14159265358979323846d0

include "fftw3.f"

integer (kind=4) Np,Nx,Ny,Nh

integer ( kind = 4 ) i,j,t,PVN,NPV,PV
real ( kind = 8 ) Lx,Ly,dx,dy,kx,ky,time,time_min,time_max,dt,mux,muy,ms,Re,CS,M,U0,n0,d,theta,C0,C1,C2,C3,G1,G2,G3
real ( kind = 8 ) Energy,Energy0,y_Energy,y_Energy0
real (kind = 8), dimension (:), allocatable :: x,y
real (kind=8), dimension(:, :), allocatable :: omega,omega_dum,n,n_dum,ux,uy,ux_dum,uy_dum
complex (kind=8), dimension(:, :), allocatable :: dn_dt_old,dn_dt_new
complex (kind=8), dimension(:, :), allocatable :: dux_dt_old,duy_dt_old,dux_dt_new,duy_dt_new
complex (kind=8), dimension(:, :), allocatable :: psik,psik_dum,omegak,omegak_dum,nk,nk_dum
complex (kind=8), dimension(:, :), allocatable :: ukx,ukx_dum,uky,uky_dum,Ek
complex (kind=8), dimension(:, :), allocatable :: nk_new,ukx_new,uky_new
integer ( kind = 8 ) plan_forward,plan_backward
integer ( kind = 4 ) thread_num,num_thread,proc_num
real ( kind = 8 ) t1,t2
common/comm/Lx,Ly,mux,muy,ms,CS,dt

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

allocate(omega(Nx,Ny))
allocate(omega_dum(Nx,Ny))
allocate(n(Nx,Ny))
allocate(n_dum(Nx,Ny))
allocate(ux(Nx,Ny))
allocate(uy(Nx,Ny))
allocate(ux_dum(Nx,Ny))
allocate(uy_dum(Nx,Ny))

allocate(dn_dt_old(Nh,Ny))
allocate(dn_dt_new(Nh,Ny))

allocate(dux_dt_old(Nh,Ny))
allocate(duy_dt_old(Nh,Ny))
allocate(dux_dt_new(Nh,Ny))
allocate(duy_dt_new(Nh,Ny))
allocate(psik(Nh,Ny))
allocate(psik_dum(Nh,Ny))
allocate(omegak(Nh,Ny))
allocate(omegak_dum(Nh,Ny))
allocate(nk(Nh,Ny))
allocate(nk_dum(Nh,Ny))
allocate(ukx(Nh,Ny))
allocate(ukx_dum(Nh,Ny))
allocate(uky(Nh,Ny))
allocate(uky_dum(Nh,Ny))
allocate(Ek(Nh,Ny))
allocate(nk_new(Nh,Ny))
allocate(ukx_new(Nh,Ny))
allocate(uky_new(Nh,Ny))

!===================== FILENAMES ==============================================

open(unit=5,file='output/System_information.dat',status='unknown')
open(unit=15,file='output/Initial_Grid_Data.dat',status='unknown')
!open(unit=20,file='INPUT.dat',status='old')  ! This input file is the file generated from Vorticity code.
open(unit=25,file='output/Initial_Grid_Data_Reproduced.dat',status='unknown')
!open(unit=30,file='Initial_Energy_Spectra.dat',status='unknown')
!open(unit=35,file='Energy_Spectra.dat',status='unknown')
open(unit=40,file='output/Energy.dat',status='unknown')

!===================== USER INPUTS ============================================

! Define Number of Threads.
proc_num = omp_get_num_procs()
thread_num = 28
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

! Co-efficient of Viscosity.
call cfg%get("resistivity","nu",mux)
mux = 0.00020d0
muy = mux

! Mass of the Fluid Element.
ms = 1.0d0

! Maximum Velocity.
U0 = 1.8850d0

! Mach Number.
M = 0.50d0

! Sound Speed.
 CS = U0/M

! Background Density.
n0 = 1.0d0

! Raynold's Number.
Re = ms*n0*U0*Lx/mux; write(5,*) "Re =", Re

! Position of Point Vortices.
d = 0.40d0*2.0d0*pi

! Number of Point Vortices.
PVN = 5

! Grid Generation.
do i = 1, Nx
  x(i)=0.0d0+real(i-1)*dx
  do j = 1, Ny
    y(j)=0.0d0+real(j-1)*dy
    ! Initial Density Distribution.
    n(i,j) = n0
    ! Initial Velocity Distribution.
    !read(20,*) G1,G2,G3,ux(i,j),uy(i,j)  ! You may define your own function here.
    if ( (x(i) - Lx/2.0d0)**2 + (y(j) -Ly/2.0d0)**2 <= (0.30d0*2.0d0*pi)**2 ) then
    omega(i,j) = 1.0d0
    else
    omega(i,j) = 0.0d0
    endif

      do PV = 1,PVN
      NPV = NPV-1
      theta = 2.0d0*pi/PVN
        if ( (x(i) - (Lx/2.0d0+d*dcos(NPV*theta)))**2 + (y(j) - (Ly/2.0d0+d*dsin(NPV*theta)))**2 <= (0.20d0)**2 ) then
        omega(i,j) = 10.0d0
        endif
      enddo

    ! Keep Backup of the Arrays for FFTW.
    n_dum(i,j) = n(i,j)
    !ux_dum(i,j) = ux(i,j)
    !uy_dum(i,j) = uy(i,j)
    omega_dum(i,j) = omega(i,j)
    ! Write Initial Density and Velocity Distribution in File.
    write(15,*) x(i),y(j),omega(i,j),n(i,j),ux(i,j),uy(i,j)
  end do
end do
 close (15)
!===================== INITIAL TIME DATA ===============================================

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, n_dum, nk, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

!  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux_dum, ukx, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_forward)

!  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, uy_dum, uky, FFTW_ESTIMATE)
!  call dfftw_execute_ (plan_forward)

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, omega_dum, omegak, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

! Evaluate Initial Vorticity Spectra.
do i = 1, Nx/2+1
  do j = 1, Ny/2
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
  !omegak(i,j) = (0.0d0,1.0d0)*kx*uky(i,j) - (0.0d0,1.0d0)*ky*ukx(i,j)
  !omegak_dum(i,j) = omegak(i,j)
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
  !omegak(i,j) = (0.0d0,1.0d0)*kx*uky(i,j) - (0.0d0,1.0d0)*ky*ukx(i,j)
  !omegak_dum(i,j) = omegak(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, uky_dum, uy, FFTW_ESTIMATE)
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
   ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
   uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
  ! Evaluate Initial Kinetic Energy.
  !Energy0 = Energy0 + (ux(i,j)**2 + uy(i,j)**2)
  !y_Energy0 = y_Energy0 + (uy(i,j)**2)
  ! Store Grid Data Reproduced from FFTW.
  write(25,*) x(i),y(j),omega(i,j),n(i,j),ux(i,j),uy(i,j)
  end do
end do
 close (25)
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
  call derive2coh (Nx,Ny,Nh,pi,time,nk,ukx,uky, &
  dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new)

! Time Solvers.
! Adams-Bashforth
  call ab2coh (Nx,Ny,Nh,pi,time,nk,ukx,uky,nk_new,ukx_new,uky_new, &
  dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new)
! Euler Method
  !call euler2coh (Nx,Ny,Nh,pi,time,nk,ukx,uky,nk_new,ukx_new,uky_new,dn_dt_new,dux_dt_new,duy_dt_new)
! Predictor-Corrector
  !call pc_nk2coh (Nx,Ny,Nh,pi,time,nk,nk_new,ukx,uky,dn_dt_new,dux_dt_new,duy_dt_new)
  !call pc_ukx2coh (Nx,Ny,Nh,pi,time,nk,ukx,ukx_new,uky,dn_dt_new,dux_dt_new,duy_dt_new)
  !call pc_uky2coh (Nx,Ny,Nh,pi,time,nk,ukx,uky,uky_new,dn_dt_new,dux_dt_new,duy_dt_new)
! Runge-Kutta (4)
  !call rk4_nk2coh (Nx,Ny,Nh,pi,time,nk,nk_new,ukx,uky,dn_dt_new,dux_dt_new,duy_dt_new) ! dux_dt_new,duy_dt_new NOT NEEDED
  !call rk4_ukx2coh (Nx,Ny,Nh,pi,time,nk,ukx,ukx_new,uky,dn_dt_new,dux_dt_new,duy_dt_new)
  !call rk4_uky2coh (Nx,Ny,Nh,pi,time,nk,ukx,uky,uky_new,dn_dt_new,dux_dt_new,duy_dt_new)

!======================================================================================

!$OMP PARALLEL SHARED(dn_dt_old,dux_dt_old,duy_dt_old,dn_dt_new,dux_dt_new,duy_dt_new),&
!$OMP & SHARED(nk,ukx,uky,nk_new,ukx_new,uky_new,nk_dum,ukx_dum,uky_dum) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx/2+1
  do j = 1,Ny
    ! Set the Variables in Proper Format for Next Time Iteration.
    dn_dt_old(i,j) = dn_dt_new(i,j)
    dux_dt_old(i,j) = dux_dt_new(i,j)
    duy_dt_old(i,j) = duy_dt_new(i,j)
    nk(i,j) = nk_new(i,j)
    ukx(i,j) = ukx_new(i,j)
    uky(i,j) = uky_new(i,j)
    ! Keep Backup of the Arrays for FFTW.
    nk_dum(i,j) = nk(i,j)
    ukx_dum(i,j) = ukx(i,j)
    uky_dum(i,j) = uky(i,j)
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

! Evaluate Vorticity in Spectral Space.

!$OMP PARALLEL SHARED(Lx,Ly,ukx,uky,omegak,omegak_dum) PRIVATE(i,j,kx,ky)
!$OMP DO

do i = 1, Nx/2+1
  do j = 1, Ny/2
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  ky = 2.0d0*pi*dfloat(j-1)/Ly
  omegak(i,j) = (0.0d0,1.0d0)*kx*uky(i,j) - (0.0d0,1.0d0)*ky*ukx(i,j)
  omegak_dum(i,j) = omegak(i,j)
  enddo
  do j = Ny/2+1,Ny
  kx = 2.0d0*pi*dfloat(i-1)/Lx
  ky = 2.0d0*pi*dfloat((j-1)-Ny)/Ly
  omegak(i,j) = (0.0d0,1.0d0)*kx*uky(i,j) - (0.0d0,1.0d0)*ky*ukx(i,j)
  omegak_dum(i,j) = omegak(i,j)
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

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, omegak_dum, omega, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

  call dfftw_destroy_plan_ (plan_backward)

! FFTW Normalisation.

!$OMP PARALLEL SHARED(n,ux,uy,omega) PRIVATE(i,j)
!$OMP DO

do i = 1,Nx
  do j = 1,Ny
  n(i,j) = n(i,j)/(dfloat(Nx)*dfloat(Ny))
  ux(i,j) = ux(i,j)/(dfloat(Nx)*dfloat(Ny))
  uy(i,j) = uy(i,j)/(dfloat(Nx)*dfloat(Ny))
  omega(i,j) = omega(i,j)/(dfloat(Nx)*dfloat(Ny))
  enddo
enddo

!$OMP END DO
!$OMP END PARALLEL

! Evaluate Energy Spectra at a Specified Time.
!if (t == 10000) then

! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file.
! !$OMP PARALLEL SHARED(ukx,uky) PRIVATE(i,j)
! !$OMP DO REDUCTION (+:Ek)

  !do i = 1,Nx/2+1
  !  do j = 1,Ny
  !  Ek(i,j) = Ek(i,j) + sqrt(abs(ukx(i,j))**2 + abs(uky(i,j))**2)
  !  write(35,*) i,j,sqrt(float((i-1)**2)+float((j-1)**2)),abs(Ek(i,j))
  !  enddo
  !enddo

! !$OMP END PARALLEL

!endif

Energy = 0.0d0; y_Energy = 0.0d0
 C0 = 0.0d0
 C1 = 0.0d0
 C2 = 0.0d0
 C3 = 0.0d0

! Try to avoid this OMP loop since it alters the sequence of i and j in Outout file.
! !$OMP PARALLEL SHARED(t,dt,x,y,omega,n,ux,uy) PRIVATE(i,j)
! !$OMP DO REDUCTION (+:Energy,y_Energy,C0,C1,C2,C3)

do i = 1,Nx
  do j = 1,Ny
  ! Evaluate Energy
  Energy = Energy + (ux(i,j)**2 + uy(i,j)**2)
  ! Evaluate Growth Rate.
  y_Energy = y_Energy + (uy(i,j)**2)
  C0 = C0 + omega(i,j)**0.0d0
  C1 = C1 + omega(i,j)**1.0d0
  C2 = C2 + omega(i,j)**2.0d0
  C3 = C3 + omega(i,j)**3.0d0
    ! Write Grid Data in State Files.
    if (mod(float(t),1.0d0/(10.0d0*dt)) == 0.0) then
      write(filename, '("output/fort.",I8.8)') t+100
      open(unit=t+100,file=filename,status='unknown')
      write(t+100,*) x(i),y(j),omega(i,j),n(i,j),ux(i,j),uy(i,j)
    endif
  enddo
enddo

! !$OMP END PARALLEL

if (mod(float(t),1.0d0/(10.0d0*dt)) == 0.0) then
  close(t+100)
endif

write(40,*) time,Energy/(2.0d0*dfloat(Nx)*dfloat(Ny)),(y_Energy)/(2.0d0*dfloat(Nx)*dfloat(Ny)), &
 C0/(dfloat(Nx)*dfloat(Ny)),C1/(dfloat(Nx)*dfloat(Ny)),C2/(dfloat(Nx)*dfloat(Ny)),C3/(dfloat(Nx)*dfloat(Ny))

enddo ! time

t2 = omp_get_wtime()

write(5,*) "Time taken for the run = ",t2 - t1

!====================================================================================

end subroutine tara2dHydOMP
