!********************************************
! ProlongationやInterpolationを計算
!********************************************

module transfer_mod
  use boundary_mod
  use calc_variables_mod
  implicit none
  
contains

  ! level -> level-1
  subroutine Prolongation(level,u,z,h)
    implicit none

    integer, intent(in) :: level 
    real(8), intent(inout) :: u(:,:), z(:,:), h(:,:)

    integer :: ic, iff

    real(8) :: ldx
    integer :: lNx

    call calc_level(level,ldx,lNx)

    do ic = 1, lNx
      iff = 2*ic
      u(ic,level-1) = u(cir(iff,lNx),level)
      z(ic,level-1) = (z(cir(iff,lNx),level) + z(cir(iff-1,lNx),level))*0.5d0
      h(ic,level-1) = (h(cir(iff,lNx),level) + h(cir(iff-1,lNx),level))*0.5d0
    end do

  end subroutine Prolongation

  ! level -> level-1
  subroutine Prolongation_defect(level,residual)
    implicit none

    integer :: level
    real(8), intent(inout) :: residual(:,:)

    integer :: ic, iff
    real(8) :: ldx !細かい方の格子幅
    integer :: lNx !細かい方の格子数

    call calc_level(level,ldx,lNx)

    do ic = 1, lNx/2
      iff = 2*ic
      residual(ic,level-1) = (residual(cir(iff,lNx),level) + residual(cir(iff-1,lNx),level))*0.5d0
    end do

  end subroutine Prolongation_defect

  ! level -> level+1
  subroutine Interpolation_defect(level,residual)
    implicit none

    integer, intent(in) :: level !fine grid number
    real(8), intent(inout)  :: residual(:,:)

    integer :: ic, iff
    real(8) :: dc1, dc2, dc3

    real(8) :: ldx
    integer :: lNx

    call calc_level(level,ldx,lNx)

    do iff = 2, lNx*2, 2
      ic = iff/2
      dc1 = residual(cir(ic  ,lNx),level)
      dc2 = residual(cir(ic-1,lNx),level)
      dc3 = residual(cir(ic+1,lNx),level)
      residual(iff,level+1) = (3*dc1+dc3)*0.25
      residual(iff-1,level+1) = (3*dc1+dc2)*0.25
    end do
    
  end subroutine Interpolation_defect
end module transfer_mod