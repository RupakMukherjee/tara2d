!=================== SUBROUTINE ADAMS-BASHFORTH ==========================

subroutine ab2csh (Nx, Ny, Nh, nu, time, dt, omegak, omegak_new, dt_omegak_old, dt_omegak_new)
implicit none
integer (kind = 4) Nx, Ny, Nh, i, j
real (kind = 8) time, dt, nu
complex (kind=8) omegak(Nh,Ny), omegak_new(Nh, Ny)
complex (kind=8) dt_omegak_new(Nh, Ny),  dt_omegak_old(Nh, Ny)

do i = 1,Nh
  do j = 1, Ny
    omegak_new(i,j) = omegak(i,j) + ( (3.0d0/2.0d0) * dt_omegak_new(i,j) - (1.0d0/2.0d0) * dt_omegak_old(i,j)) * dt
  enddo
enddo

end subroutine ab2csh
