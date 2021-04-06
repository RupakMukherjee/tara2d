program TARA2D

use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

integer dim, iargc, numb

character (len=90) :: filename
character (len=90) :: arch, mode, domain
character (len=32) :: arg

numb = iargc()

call getarg(1, arg)

if (numb == 0) then
  arg = "input.ini"
  write(*,*) "No input file provided. Default input.ini is being used."
endif

cfg = parse_cfg(arg)

call cfg%get("dimension","dim",dim)
call cfg%get("architecture","arch",arch)
call cfg%get("parallel","mode",mode)
call cfg%get("type","domain",domain)

!DIMENSION = 1
if (dim == 1) then
  if (arch == "cpu") then
    if (mode == "serial") then
      write(*,*) "Congrats! Your ONE dimensional SERIAL HYDRODYNAMIC code is now running!"
      call tara1dHydSer(arg)
    elseif (mode == "openmp") then
      write(*,*) "Under preparation. Please wait!"
    elseif (mode == "mpi") then
      write(*,*) "Under preparation. Please wait!"
    endif
  elseif (arch == "single-gpu") then
    write(*,*) "Under preparation. Please wait!"
  elseif (arch == "multi-gpu") then
    write(*,*) "Under preparation. Please wait!"
  endif
!DIMENSION = 2
elseif (dim == 2) then
  if (arch == "cpu") then
    if (mode == "serial") then
      write(*,*) "Congrats! Your TWO dimensional SERIAL HYDRODYNAMIC code is now running!"
      call tara2dHydSer(arg)
    elseif (mode == "openmp" .and. domain == "hydro") then
      write(*,*) "Congrats! Your TWO dimensional OPENMP parallel HYDRODYNAMIC code is now running!"
      call tara2dHydOMP(arg)
    elseif (mode == "openmp" .and. domain == "mhd") then
      write(*,*) "Congrats! Your TWO dimensional OPENMP parallel MHD code is now running!"
      call tara2dMHDOMP(arg)
    elseif (mode == "mpi") then
      write(*,*) "Under preparation. Please wait!"
    endif
  elseif (arch == "single-gpu") then
    write(*,*) "Under preparation. Please wait!"
  elseif (arch == "multi-gpu") then
    write(*,*) "Under preparation. Please wait!"
  endif
!DIMENSION = 3
elseif (dim == 3) then
  if (arch == "cpu") then
    if (mode == "serial") then
      write(*,*) "Under preparation. Please wait!"
    elseif (mode == "openmp" .and. domain == "mhd") then
      write(*,*) "Congrats! Your THREE dimensional OPENMP parallel MHD code is now running!"
      call tara3dMHDOMP(arg)
    elseif (mode == "mpi") then
      write(*,*) "Under preparation. Please wait!"
    endif
  elseif (arch == "single-gpu") then
    write(*,*) "Under preparation. Please wait!"
  elseif (arch == "multi-gpu") then
    write(*,*) "Under preparation. Please wait!"
  endif
else
  write(*,*) "Incompatible Input file. Please check TARA user-manual."
endif

end program TARA2D
