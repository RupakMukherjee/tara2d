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

if (dim == 1 .and. arch == "cpu" .and. mode == "serial") then
  write(*,*) "Congrats! Your code is running!"
  call tara1dHydSer(arg)
elseif (dim == 1 .and. arch == "cpu" .and. mode == "openmp") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 1 .and. arch == "cpu" .and. mode == "mpi") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 1 .and. arch == "single-gpu") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 1 .and. arch == "multi-gpu") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 2 .and. arch == "cpu" .and. mode == "serial") then
  write(*,*) "Congrats! Your code is running!"
  call tara2dHydSer(arg)
elseif (dim == 2 .and. arch == "cpu" .and. mode == "openmp" .and. domain == "hydro") then
  write(*,*) "Congrats! Your code is running!"
  call tara2dHydOMP(arg)
elseif (dim == 2 .and. arch == "cpu" .and. mode == "openmp" .and. domain == "mhd") then
  write(*,*) "Congrats! Your code is running!"
  call tara2dMHDOMP(arg)
elseif (dim == 2 .and. arch == "cpu" .and. mode == "mpi") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 2 .and. arch == "single-gpu") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 2 .and. arch == "multi-gpu") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == "cpu" .and. mode == "serial") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == "cpu" .and. mode == "openmp" .and. domain == "mhd") then
  write(*,*) "Congrats! Your code is running!"
  call tara3dMHDOMP(arg)
elseif (dim == 3 .and. arch == "cpu" .and. mode == "mpi") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == "single-gpu") then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == "multi-gpu") then
  write(*,*) "Under preparation. Please wait!"
else
  write(*,*) "Incompatible Input file. Please check TARA user-manual."
endif

end program TARA2D
