module initialization
  use constant
  use boundary_mod
  use calc_variables_mod
  implicit none
  
contains
  subroutine initialize(u,z,h)
    implicit none

    real(8), intent(inout) :: u(:,:), z(:,:), h(:,:)
    integer :: i, l

    do l = 1, num_levels
      do i = 1, Nx
        u(i,l) = 0.d0
        z(i,l) = 0.d0
        h(i,l) = 1.d3
      enddo
    enddo
    do i = Nx/5*2, Nx/5*3
      h(i,num_levels) = 1.d3 - 10.d0
    end do
    
  end subroutine initialize

  subroutine init_coef(Au,Az,b,residual)
    implicit none
    
    real(8), intent(inout) :: Au(:,:), Az(:,:), b(:,:), residual(:,:)
    integer :: i, l

    do l = 1, num_levels
      do i = 1, Nx
        Au(i,l) = 0.d0
        Az(i,l) = 0.d0
        b(i,l) = 0.d0
        residual(i,l) = 0.d0
      enddo
    enddo

  end subroutine init_coef

end module initialization