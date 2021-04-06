
!===================================================================================
!=================== SUBROUTINE RUNGE - KUTTA 4 ====================================
!===================================================================================

subroutine rk42css(Nx,Ny,Nh,alpha,nu,time,dt,omegak,omegak_new,dt_omegak_old,dt_omegak_new)
implicit none
integer ( kind = 4 ) Nx,Ny,Nh, i, j
real ( kind = 8 ) time,dt,alpha,nu,lambda
complex ( kind = 8 ) komega1(Nh,Ny),komega2(Nh,Ny),komega3(Nh,Ny),komega4(Nh,Ny)
complex ( kind = 8 ) omegak(Nh,Ny),omegak_new(Nh,Ny),dum_omegak(Nh,Ny),dt_omegak(Nh,Ny)
complex ( kind = 8 ) dt_omegak_old(Nh,Ny)
complex ( kind = 8 ) dt_omegak_new(Nh,Ny)

do i = 1,Nh
  do j = 1,Ny
    komega1(i,j) = dt_omegak(i,j)
    dum_omegak(i,j) = omegak(i,j) + komega1(i,j)*dt/2.0
  end do
end do

  call derive2css(Nx,Ny,Nh,alpha,nu,time+dt/2.0,dum_omegak,omegak_new,dt_omegak_old,dt_omegak_new)

do i = 1,Nh
  do j = 1,Ny
    komega2(i,j) = dt_omegak(i,j)
    dum_omegak(i,j) = omegak(i,j) + komega2(i,j)*dt/2.0
  end do
end do

  call derive2css(Nx,Ny,Nh,alpha,nu,time+dt/2.0,dum_omegak,omegak_new,dt_omegak_old,dt_omegak_new)

do i = 1,Nh
  do j = 1,Ny
    komega3(i,j) = dt_omegak(i,j)
    dum_omegak(i,j) = omegak(i,j) + komega3(i,j)*dt/2.0
  end do
end do

  call derive2css(Nx,Ny,Nh,alpha,nu,time+dt,dum_omegak,omegak_new,dt_omegak_old,dt_omegak_new)

do i = 1,Nh
  do j = 1,Ny
    komega4(i,j) = dt_omegak(i,j)
    omegak_new(i,j) = omegak(i,j) + dt/6.0d0*(komega1(i,j) + 2.0d0*komega2(i,j) + 2.0d0*komega3(i,j) + komega4(i,j))
!    omegak_new(i,j) = omegak(i,j) + komega1(i,j)*dt     ! EULER SOLVER
  end do
end do

return

end subroutine rk42css
