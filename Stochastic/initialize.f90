module initialization
  use constant
  use boundary_mod
  use calc_variables_mod
  implicit none
  
contains
  subroutine initialize(u,z,h)
    implicit none

    real(8), intent(inout) :: u(:,:), z(:,:), h(:,:)

    u(:,:) = 0.d0
    z(:,:) = 0.d0
    h(:,:) = 1.d4
    
  end subroutine initialize

end module initialization