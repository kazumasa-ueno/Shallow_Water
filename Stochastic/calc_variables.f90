!********************************************
! 色々な変数を計算するためのモジュール
!  calc_u:       zの値からuを計算
!  calc_Au:      Auの値を計算
!  calc_Az:      Azの値を計算
!  calc_b:        bの値を計算
!  calc_Fu:      移流元のuの値を計算
!  inner_u:      uの内装を計算
!  z_frac:        分数インデックスでの値を計算
!  calc_res:      残差を計算
!********************************************

module calc_variables_mod
  use constant
  use boundary_mod
  implicit none

  interface z_frac
    module procedure z_frac_2, z_frac_list
  end interface z_frac
  
contains

  subroutine calc_u(level,u,z)
    implicit none
    
    integer, intent(in) :: level
    real(8), intent(in) :: z(:,:)
    real(8), intent(inout) :: u(:,:)

    integer :: i
    real(8) :: Fu 
    real(8), allocatable :: u_b(:) !u before update
    real(8) :: ldx
    integer :: lNx
    
    call calc_level(level,ldx,lNx)
    allocate(u_b(lNx))
    
    do i = 1, lNx
      Fu = calc_Fuu(level,i,u)
      u_b(i) = (Fu - g*(dt/ldx)*(z(cir(i+1,lNx),level)-z(cir(i,lNx),level)) &
      & + dt*XForce)
      ! & + nu*dt/dx**2*(u(cir(i+1,lNx),level)-2*u(cir(i,lNx),level)+u(cir(i-1,lNx),level)) + dt*XForce)
    end do
    u(:,level) = u_b(:)
    
    deallocate(u_b)
    
    ! call boundary_u(u,Nx,times)      
    
  end subroutine calc_u
  
  subroutine calc_Au(level,z,h,Au)
    implicit none
    
    integer, intent(in) :: level
    real(8), intent(in) :: z(:,:), h(:,:)
    real(8), intent(inout) :: Au(:,:)
    
    real(8) :: ldx
    integer :: lNx
    integer :: i

    call calc_level(level,ldx,lNx)
    
    do i = 1, lNx
      Au(i,level) = g * (dt/ldx)**2 * (z_frac(z(i:i+1,level)) + z_frac(h(i:i+1,level)))
    end do

  end subroutine calc_Au

  subroutine calc_Az(level,Au,Az)
    implicit none
    
    integer, intent(in) :: level
    real(8), intent(in) :: Au(:,:)
    real(8), intent(inout) :: Az(:,:)

    real(8) :: ldx
    integer :: lNx
    integer :: i
    
    call calc_level(level,ldx,lNx)
    
    do i = 1, lNx
      Az(i,level) = 1 + Au(cir(i,lNx),level) + Au(cir(i-1,lNx),level)
    end do
    
  end subroutine calc_Az
  
  subroutine calc_b(level,u,z,h,b)
    implicit none
    
    integer, intent(in) :: level
    real(8), intent(in) :: u(:,:), z(:,:), h(:,:)
    real(8), intent(inout) :: b(:,:)
    
    integer :: i
    real(8) :: Fu(2)
    real(8) :: ldx
    integer :: lNx
    
    call calc_level(level,ldx,lNx)
    do i = 1, lNx
      Fu(1) = calc_Fuu(level,i,u)
      Fu(2) = calc_Fuu(level,i-1,u)
      b(i,level) = z(cir(i,lNx),level) - (dt/ldx) * ( &
      & (z_frac(z(cir(i,lNx),level),z(cir(i+1,lNx),level))+z_frac(h(cir(i,lNx),level),h(cir(i+1,lNx),level))) &
      & * (Fu(1)+dt*XForce) &
      & - (z_frac(z(cir(i-1,lNx),level),z(cir(i,lNx),level))+z_frac(h(cir(i-1,lNx),level),h(cir(i,lNx),level))) &
      & * (Fu(2)+dt*XForce) )
      ! & * (Fu(1)+ nu/dx*(u(cir(i+1,lNx),level)-2*u(cir(i,lNx),level)+u(cir(i-1,lNx),level)) + dt*XForce) &
      ! & - (z_frac(z(cir(i-1,lNx),level),z(cir(i,lNx),level))+z_frac(h(cir(i-1,lNx),level),h(cir(i,lNx),level))) &
      ! & * (Fu(2)+nu/dx*(u(cir(i,lNx),level)-2*u(cir(i-1,lNx),level)+u(cir(i-2,lNx),level))+dt*XForce) )
    end do

  end subroutine calc_b


  ! calculate Fu(i+1/2) actuarlly(u(i))
  ! 0 <= i <= Nx
  real(8) function calc_Fuu(level,i,u)
    implicit none
    
    integer, intent(in) :: level, i
    real(8), intent(in) :: u(:,:)

    integer :: s, smax
    real(8) :: x, u_s

    real(8) :: ldx
    integer :: lNx

    call calc_level(level,ldx,lNx)
    
    smax = int(dt/dtau)
    x = (i+0.5d0)*ldx
    u_s = u(cir(i,lNx),level)
    do s = 1, smax
      x = x - dtau*u_s
      call inner_u(x,u_s,u(:,level),ldx,lNx)
    end do
    calc_Fuu = u_s

  end function calc_Fuu

  subroutine inner_u(x,u_s,u,dx,Nx)
    implicit none

    integer, intent(in) :: Nx
    real(8), intent(in) :: x, u(:), dx
    real(8), intent(out) :: u_s
    integer :: iu
    real(8) :: pu

    iu = int(x/dx-0.5d0)
    pu = x/dx - 0.5d0 - iu
    u_s = (1-pu)*u(cir(iu,Nx)) + pu*u(cir(iu+1,Nx))

  end subroutine inner_u

  !calculate the value when z's index is fractional

  real(8) function z_frac_2(z1,z2)
    implicit none
    
    real(8), intent(in) :: z1, z2

    z_frac_2 = (z1+z2)*0.5d0

  end function z_frac_2

  real(8) function z_frac_list(z)
    implicit none

    real(8), intent(in) :: z(2)

    z_frac_list = (z(1)+z(2))*0.5d0

  end function z_frac_list
  
  subroutine calc_res(z,Au,Az,b,Res)
    implicit none
    
    real(8), intent(in) :: z(:,:), Au(:,:), Az(:,:), b(:,:)
    real(8), intent(out) :: Res

    integer :: i

    Res = 0.d0
    do i = 1, Nx
      Res = Res + (b(i,num_levels) + Au(cir(i-1,Nx),num_levels)*z(cir(i-1,Nx),num_levels) &
      & + Au(cir(i,Nx),num_levels)*z(cir(i+1,Nx),num_levels) - Az(cir(i,Nx),num_levels)*z(cir(i,Nx),num_levels))**2
    end do
    Res = (Res**0.5d0)/Nx

  end subroutine calc_res

  subroutine calc_level(level,ldx,lNx)
    implicit none

    integer, intent(in) :: level
    real(8), intent(out) :: ldx
    integer, intent(out) :: lNx
    
    ldx = dx*2**(num_levels-level)
    lNx = Nx/2**(num_levels-level)

  end subroutine calc_level


end module calc_variables_mod
