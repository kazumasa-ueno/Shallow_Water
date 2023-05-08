module initialization
  use constant
  use boundary_mod
  use calc_variables_mod
  implicit none
  
contains
  subroutine initialize(u,v,z,gamma,h)
    implicit none

    real(8), intent(out) :: u(0:Nx), v(0:Nx+1), z(0:Nx+1), gamma(0:Nx+1), h(0:Nx+1)

    integer :: i

    u(:) = 0.d0
    v(:) = 0.d0
    z(:) = 0.d0
    h(:) = 1.d0
    z(:) = 0.d0 
    do i = Nx/5*2, Nx/5*3
      h(i) = 1.d3 - 10.d0
    end do
    call calc_gamma(u,v,z,h,gamma,Nx)
    call boundary_u(u,Nx,1)
    
  end subroutine initialize

end module initialization