program TARA2D

use cfgio_mod, only: cfg_t, parse_cfg

implicit none

type(cfg_t):: cfg

integer mode

character (len=90) :: filename

cfg = parse_cfg("input.ini")

call cfg%get("architecture","mode",mode)

if (mode == 1) then
  call tara2dHydSer()
endif

end program TARA2D
