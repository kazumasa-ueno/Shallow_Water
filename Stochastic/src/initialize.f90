module initialization
  use constant
  use structure
  use boundary_mod
  use transfer_mod
  use calc_variables_mod
  implicit none
  
contains
  subroutine initialize()
    implicit none

    integer :: i, l

    real(8) :: ldx
    integer :: lNx

    do l = 1, num_levels
      do i = 1, Nx
        u(i,l) = 0.d0
        z(i,l) = 0.d0
        h(i,l) = 1.d4
        ! XForce(i,l) = 0.d0
      enddo
      call calc_level(l,ldx,lNx)
      ! do i = 1, lNx
      !   h(i,l) = 1.d4 + 10.d0*sin(2*pi*i/lNx) + 10.d0*sin(2.d0*2*pi*i/lNx) + 10.d0*sin(3.d0*2*pi*i/lNx)
      ! enddo
      call calc_XForce(l)
    enddo
    do i = Nx/256+Nx/2, Nx/256*2+Nx/2
      h(i,num_levels) = 1.d4 + 1000.d0
      ! z(i,num_levels) = 0.01
    end do
    ! do i = Nx/256+Nx/2+4, Nx/256*2+Nx/2+4
    !   ! h(i,num_levels) = 1.d4 + 1000.d0
    !   z(i,num_levels) = -0.01
    ! end do
    ! call calc_XForce(2)
    
  end subroutine initialize

  subroutine init_coef()
    implicit none
    
    integer :: i, l

    do l = 1, num_levels
      do i = 1, Nx
        Au(i,l) = 0.d0
        Az(i,l) = 0.d0
        b(i,l) = 0.d0
        residual(i,l) = 0.d0
        rhs(i,l) = 0.d0
      enddo
    enddo

  end subroutine init_coef

end module initialization