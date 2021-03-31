program TARA2D

use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

integer dim, arch, mode

character (len=90) :: filename

cfg = parse_cfg("input.ini")

call cfg%get("dimension","dim",dim)
call cfg%get("architecture","arch",arch)
call cfg%get("parallel","mode",mode)

if (dim == 1 .and. arch == 1 .and. mode == 1) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 1 .and. arch == 1 .and. mode == 2) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 1 .and. arch == 1 .and. mode == 3) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 1 .and. arch == 2) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 1 .and. arch == 3) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 2 .and. arch == 1 .and. mode == 1) then
  write(*,*) "Congrats! Your code is running!"
  call tara2dHydSer()
elseif (dim == 2 .and. arch == 1 .and. mode == 2) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 2 .and. arch == 1 .and. mode == 3) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 2 .and. arch == 2) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 2 .and. arch == 3) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == 1 .and. mode == 1) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == 1 .and. mode == 2) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == 1 .and. mode == 3) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == 2) then
  write(*,*) "Under preparation. Please wait!"
elseif (dim == 3 .and. arch == 3) then
  write(*,*) "Under preparation. Please wait!"
else
  write(*,*) "Incompatible Input file. Please check TARA user-manual."
endif

end program TARA2D
