module initialization
  use constant
  use boundary_mod
  use calc_variables_mod
  implicit none
  
contains
  subroutine initialize(u,z,h)
    implicit none

    real(8), intent(inout) :: u(:,:), z(:,:), h(:,:)
    integer :: i

    u(:,:) = 0.d0
    z(:,:) = 0.d0
    ! h(:,:) = 1.d4
    h(:,:) = 1.d3
    do i = Nx/5*2, Nx/5*3
      h(i,num_levels) = 1.d3 - 10.d0
    end do
    
  end subroutine initialize

end module initialization