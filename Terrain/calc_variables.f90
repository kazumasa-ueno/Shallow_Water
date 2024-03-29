!********************************************
! 色々な変数を計算するためのモジュール
!  calc_u:       zの値からuを計算
!  calc_v:       zの値からvを計算
!  calc_gamma:       zの値からgammaを計算
!  calc_Au:      Auの値を計算
!  calc_Az:      Azの値を計算
!  calc_b:        bの値を計算
!  calc_Fu:      移流元のuの値を計算
!  calc_Fv:      移流元のvの値を計算
!  inner_u:      uの内装を計算
!  inner_v:      vの内装を計算
!  z_frac:        分数インデックスでの値を計算
!  calc_res:      残差を計算
!********************************************

module calc_variables_mod
  use constant
  use boundary_mod
  implicit none
  
contains

  subroutine calc_u(u,v,z,gamma,times,Nx)
    implicit none
    
    integer, intent(in) :: Nx, times
    real(8), intent(in) :: v(0:Nx+1), z(0:Nx+1), gamma(0:Nx+1)
    real(8), intent(inout) :: u(0:Nx)

    integer :: i
    real(8) :: Fu, Fv, u_b(0:Nx) !u before update

    u_b(:) = u(:)

    do i = 0, Nx
      Fu = calc_Fuu(i,u_b,dt,dx,dtau,Nx)
      Fv = calc_Fvu(i,u_b,v,dt,dx,dtau,Nx)
      u(i) = (Fu - g*(dt/dx)*(z(cir_i(i+1,Nx))-z(cir_i(i,Nx))) &
      & + f*dt*Fv + dt*XForce) / (1+(gamma(cir_i(i,Nx))+gamma(cir_i(i+1,Nx)))*0.5d0*dt)
    end do

    ! call boundary_u(u,Nx,times)      

  end subroutine calc_u

  subroutine calc_v(u,v,gamma,Nx)
    implicit none

    integer, intent(in) :: Nx
    real(8), intent(in) :: u(0:Nx), gamma(0:Nx+1)
    real(8), intent(inout) :: v(0:Nx+1)

    integer :: i
    real(8) :: Fu, Fv, v_b(0:Nx+1) !v before update

    v_b(:) = v(:)

    do i = 0, Nx+1
      Fu = calc_Fuc(i,u,dt,dx,dtau,Nx)
      Fv = calc_Fvc(i,u,v,dt,dx,dtau,Nx)
      v(i) = (Fv - f*dt*Fu) / (1+gamma(cir_i(i,Nx))*dt)
    end do

    ! call boundary_v(v,Nx)
    
  end subroutine calc_v

  subroutine calc_gamma(u,v,z,h,gamma,Nx)
    implicit none

    integer, intent(in) :: Nx
    real(8), intent(in) :: u(0:Nx), v(0:Nx+1), z(0:Nx+1), h(0:Nx+1)
    real(8), intent(out) :: gamma(0:Nx+1)

    integer :: i

    do i = 0, Nx+1
      gamma(i) = g*sqrt(((u(cir_i(i-1,Nx))+u(cir_i(i,Nx)))*0.5d0)**2+v(i)**2) &
      & /(Cz**2*(h(i)+z(i)))
    enddo
    ! gamma = 0.d0
    
  end subroutine calc_gamma

  subroutine calc_Au(z,gamma,h,g,dt,dx,Nx,Au)
    implicit none
    
    integer, intent(in) :: Nx
    real(8), intent(in) :: z(0:Nx+1), gamma(0:Nx+1), h(0:Nx+1), g, dt, dx
    real(8), intent(out) :: Au(0:Nx)

    integer :: i

    do i = 0, Nx
      Au(i) = g * (dt/dx)**2 * (z_frac(z(i:i+1)) + z_frac(h(i:i+1))) / (1 + z_frac(gamma(i:i+1))*dt)
    end do

  end subroutine calc_Au

  subroutine calc_Az(Au,Nx,Az)
    implicit none
    
    integer, intent(in) :: Nx
    real(8), intent(in) :: Au(0:Nx)
    real(8), intent(out) :: Az(1:Nx)

    integer :: i, j

    do i = 1, Nx
      Az(i) = 1 + Au(i) + Au(i-1)
    end do

  end subroutine calc_Az

  subroutine calc_b(u,v,z,gamma,h,f,dt,dx,dtau,Nx,b)
    implicit none
    
    integer, intent(in) :: Nx
    real(8), intent(in) :: u(0:Nx), v(0:Nx+1), z(0:Nx+1), h(0:Nx+1), f, gamma(0:Nx+1)
    real(8), intent(in) :: dt, dx, dtau
    real(8), intent(out) :: b(1:Nx)

    integer :: i
    real(8) :: Fu(3), Fv(3)

    do i = 1, Nx
      Fu(1) = calc_Fuu(i,u,dt,dx,dtau,Nx)
      Fu(2) = calc_Fuu(i-1,u,dt,dx,dtau,Nx)
      Fv(1) = calc_Fvu(i,u,v,dt,dx,dtau,Nx)
      Fv(2) = calc_Fvu(i-1,u,v,dt,dx,dtau,Nx)
      b(i) = z(i) - (dt/dx) * ( &
      & (z_frac(z(i:i+1))+z_frac(h(i:i+1))) / (1+z_frac(gamma(i:i+1))*dt) &
      & * (Fu(1)+f*dt*Fv(1)+dt*XForce) &
      & - (z_frac(z(i-1:i))+z_frac(h(i-1:i))) / (1+z_frac(gamma(i-1:i))*dt) &
      & * (Fu(2)+f*dt*Fv(2)+dt*XForce) )
    end do

  end subroutine calc_b

  !----------------------------!
  !!!!!!逆流の時は修正必要!!!!!!!!!
  !----------------------------!

  ! calculate Fu(i+1/2) actuarlly(u(i))
  ! 0 <= i <= Nx
  real(8) function calc_Fuu(i,u,dt,dx,dtau,Nx)
    implicit none
    
    integer, intent(in) :: Nx, i
    real(8), intent(in) :: u(0:Nx), dt, dx, dtau

    integer :: s, smax
    real(8) :: x, u_s

    smax = int(dt/dtau)
    x = (i+0.5d0)*dx
    u_s = u(cir_i(i,Nx))
    do s = 1, smax
      x = x - dtau*u_s
      ! if(x<0) then
      !   u_s = u_upstream
      !   exit
      ! end if
      call inner_u(x,u_s,u,dx,Nx)
    end do
    calc_Fuu = u_s

    ! calc_Fuu = u(i)
  end function calc_Fuu

  ! calculate Fu(i) actually(u(i-1/2))
  ! 1 <= i <= Nx
  real(8) function calc_Fuc(i,u,dt,dx,dtau,Nx)
    implicit none
    
    integer, intent(in) :: Nx, i
    real(8), intent(in) :: u(0:Nx), dt, dx, dtau

    integer :: s, smax
    real(8) :: x, u_s

    smax = int(dt/dtau)
    x = i*dx
    ! u_s = z_frac(u(i-1:i))
    u_s = (u(cir_i(i-1,Nx))+cir_i(i,Nx))*0.5d0
    do s = 1, smax
      x = x - dtau*u_s
      ! if(x<0) then
      !   u_s = u_upstream
      !   exit
      ! end if
      call inner_u(x,u_s,u,dx,Nx)
    end do
    calc_Fuc = u_s

    ! calc_Fuc = u(i)
  end function calc_Fuc

  ! calculate Fv(i)
  ! 1 <= i <= Nx
  real(8) function calc_Fvc(i,u,v,dt,dx,dtau,Nx)
    implicit none
    
    integer, intent(in) :: Nx, i
    real(8), intent(in) :: u(0:Nx), v(0:Nx+1), dt, dx, dtau

    integer :: s, smax
    real(8) :: x, u_s, v_s

    smax = int(dt/dtau)
    x = i*dx
    u_s = (u(cir_i(i,Nx))+u(cir_i(i-1,Nx)))*0.5d0
    v_s = v(i)
    do s = 1, smax
      x = x - dtau*u_s
      ! if(x<0) then
      !   v_s = v_upstream
      !   exit
      ! endif
      call inner_u(x,u_s,u,dx,Nx)
      call inner_v(x,v_s,u,v,dx,Nx)
    end do
    calc_Fvc = v_s

    ! calc_Fvc = v(i)
  end function calc_Fvc

  ! calculate Fv(i+1/2)
  ! 0 <= i <= Nx
  real(8) function calc_Fvu(i,u,v,dt,dx,dtau,Nx)
    implicit none
    
    integer, intent(in) :: Nx, i
    real(8), intent(in) :: u(0:Nx), v(0:Nx+1), dt, dx, dtau

    integer :: s, smax
    real(8) :: x, u_s, v_s

    smax = int(dt/dtau)
    x = (i+0.5d0)*dx
    u_s = u(cir_i(i,Nx))
    v_s = (v(i)+v(i+1))*0.5d0
    do s = 1, smax
      x = x - dtau*u_s
      ! if(x<0) then
      !   v_s = v_upstream
      !   exit
      ! endif
      call inner_u(x,u_s,u,dx,Nx)
      call inner_v(x,v_s,u,v,dx,Nx)
    end do
    calc_Fvu = v_s

    ! calc_Fvu = v(i)
  end function calc_Fvu

  subroutine inner_u(x,u_s,u,dx,Nx)
    implicit none

    integer, intent(in) :: Nx
    real(8), intent(in) :: x, u(0:Nx), dx
    real(8), intent(out) :: u_s
    integer :: iu, iv
    real(8) :: pu, xv, yv, pv, qv

    iu = int(x/dx-0.5d0)
    pu = x/dx - 0.5d0 - iu
    u_s = (1-pu)*u(cir_i(iu,Nx)) + pu*u(cir_i(iu+1,Nx))

  end subroutine inner_u

  subroutine inner_v(x,v_s,u,v,dx,Nx)
    implicit none

    integer, intent(in) :: Nx
    real(8), intent(in) :: x, u(0:Nx), v(0:Nx+1), dx
    real(8), intent(out) :: v_s
    integer :: iu, iv
    real(8) :: pu, pv

    iv = int(x/dx)
    pv = x/dx - iv
    v_s = (1-pv)*v(cir_i(iv,Nx)) + pv*v(cir_i(iv+1,Nx))

  end subroutine inner_v

  !calculate the value when z's index is fractional
  real(8) function z_frac(z)
    implicit none

    real(8), intent(in) :: z(2)

    z_frac = (z(1)+z(2))*0.5d0

  end function z_frac
  
  subroutine calc_res(z,Au,Az,b,Nx,Res)
    implicit none
    
    integer, intent(in) :: Nx
    real(8), intent(in) :: z(0:Nx+1), Au(0:Nx), Az(1:Nx), b(1:Nx)
    real(8), intent(out) :: Res

    integer :: i

    Res = 0.d0
    do i = 1, Nx
      Res = Res + (b(i) + Au(i-1)*z(i-1) + Au(i)*z(i+1) - Az(i)*z(i))**2
    end do
    Res = (Res**0.5d0)/Nx

  end subroutine calc_res


end module calc_variables_mod
