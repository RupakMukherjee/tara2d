subroutine tara2dHydOMP(arg)

use omp_lib
use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

real (kind=8), parameter :: pi = 3.14159265358979323846d0

include "fftw3.f"

integer (kind=4) Np,Nx,Ny,Nh

integer (kind=4) i,j
integer (kind=4) thread_num,num_thread,proc_num,iret
real dx,dy,Lx,Ly

real (kind = 8), dimension (:), allocatable :: x,y

real (kind=8), dimension(:, :), allocatable :: ux,ux_dum

complex (kind=8), dimension(:, :), allocatable :: ukx, ukx_dum

integer (kind=8) plan_forward, plan_backward
integer ( kind = 4 ) seed

character (len=90) :: filename
character (len=32) :: arg

cfg = parse_cfg(arg)

call cfg%get("grid","Nx",Nx)
call cfg%get("grid","Ny",Ny)

Nh = (Nx/2) + 1

seed = 123456789

allocate(x(Nx))
allocate(y(Ny))
allocate(ux(Nx,Ny))
allocate(ux_dum(Nx,Ny))
allocate(ukx(Nh,Ny))
allocate(ukx_dum(Nh,Ny))

  proc_num = omp_get_num_procs()
  thread_num = proc_num
  call omp_set_num_threads (thread_num)
  print*, thread_num

  call cfg%get("length","Lx",Lx)
  call cfg%get("length","Ly",Ly)

dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

do i = 1,Nx
x(i) = (i-1)*dx
  do j = 1,Ny
  y(j) = (j-1)*dy
     ux(i,j) = ran(seed)
     seed=ux(i,j)*seed
  enddo
enddo

  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ux, ukx, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)

do i = 1,Nh
  do j = 1,Ny
  ukx_dum(i,j) = ukx(i,j)
!  write(*,*) i,j,ukx_dum(i,j)
  enddo
enddo
  call dfftw_init_threads(iret)
  call dfftw_plan_with_nthreads(thread_num)
  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ukx_dum, ux_dum, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

do i = 1, Nx
  do j = 1,Ny
    write(filename, '("output/fort.",I8.8)') 100
    open(unit=100,file=filename,status='unknown')
    write(100,*) i,j, ux(i,j), ux_dum(i,j)/(dfloat(Nx)*dfloat(Ny))
  enddo
enddo


end subroutine tara2dHydOMP
