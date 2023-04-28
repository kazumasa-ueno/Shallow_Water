!********************************************
! 境界条件を計算
!********************************************

module boundary_mod
  use constant
  implicit none
  
contains

  subroutine boundary(z,Nx,dt,dx,gamma,f,g,Fu,Fv,u)
    implicit none

    integer, intent(in) :: Nx
    real(8), intent(in) :: dt,dx,gamma(0:Nx+1),f,Fu,Fv,g,u(0:Nx)
    real(8), intent(inout) :: z(0:Nx+1)

    ! z(0) = z(1) - dx/(g*dt)*((1+(gamma(0)+gamma(1))*0.5d0*dt)*u(0)-Fu-f*dt*Fv)
    ! z(0) = z(1)
    z(Nx+1) = z(0)
    
  end subroutine boundary
  
  subroutine boundary_u(u,Nx,times)
    implicit none
    
    integer, intent(in) :: Nx, times
    real(8), intent(inout) :: u(0:Nx)
    
    u(0) = u_upstream
    
  end subroutine boundary_u
  
  subroutine boundary_v(v,Nx)
    implicit none

    integer, intent(in) :: Nx
    real(8), intent(inout) :: v(0:Nx+1)

    v(0) = v_upstream
    v(Nx+1) = v(Nx)
    
  end subroutine boundary_v

  integer function cir_i(i,Nx)
    implicit none
    
    integer, intent(in) :: i, Nx

    if(i>Nx) then
      cir_i = i - Nx - 1
    else if(i<0) then
      cir_i = i + Nx + 1
    else
      cir_i = i
    endif

  end function cir_i

end module boundary_mod