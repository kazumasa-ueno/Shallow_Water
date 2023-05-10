module solve
  use calc_variables_mod
  use transfer_mod
  implicit none

contains
  recursive subroutine MGCYC(k,u,z,h,Au,Az,b,residual,cyc,times)
  implicit none

  integer, intent(in) :: k
  real(8), intent(inout) :: Au(:,:), Az(:,:), b(:,:), u(:,:), h(:,:), residual(:,:)
  real(8), intent(inout) :: z(:,:)
  integer, intent(in) :: cyc, times

  integer :: nt, i

  real(8) :: tmp_res(Nx)

  real(8) :: ldx
  integer :: lNx

  call calc_level(k,ldx,lNx)

  !Presmoothing
  do nt = 1, nu1
    call smooth(k,residual,Au,Az,b)
  end do

  !Coarse grid correction
  !Compute the defect
  tmp_res(:) = residual(:,k)
  do i = 1, lNx
    residual(i,k) = b(cir(i,lNx),k) + Au(cir(i-1,lNx),k)*tmp_res(cir(i-1,lNx)) &
    & + Au(cir(i,lNx),k)*tmp_res(cir(i+1,lNx)) - Az(cir(i,lNx),k)*tmp_res(cir(i,lNx))
  end do

  !Restrict the defect
  call Prolongation(k,u,z,h)
  call Prolongation_defect(k,residual)
  call calc_Au(k,z,h,Au)
  call calc_Az(k,Au,Az)
  call calc_b(k,u,z,h,b)
  ! write(*,*) 'ok'

  !Compute an approximate solution v of the defect equation on k-1
    if(k==1) then
      do nt = 1, 1
        call smooth(k,residual,Au,Az,b)
      end do
    else
      call MGCYC(k-1,u,z,h,Au,Az,b,residual,cyc,times)
    end if

    !Interpolate the correction
    call Interpolation_defect(k,residual)

    !Compute the corrected approximation on k
    z(:,k) = z(:,k) + residual(:,k)

    !Postsmoothing
    do nt = 1, nu2
      call smooth(k,z,Au,Az,b)
    end do

  end subroutine MGCYC

  subroutine smooth(level,z,Au,Az,b)
    implicit none
    
    integer, intent(in) :: level
    real(8), intent(in) :: Au(:,:), Az(:,:), b(:,:)
    real(8), intent(inout) :: z(:,:)

    integer :: i
    real(8) :: ldx
    integer :: lNx

    call calc_level(level,ldx,lNx)

    do i = 1, lNx
    ! do i = lNx, 1, -1
      z(i,level) = (b(cir(i,lNx),level)+Au(cir(i-1,lNx),level)*z(cir(i-1,lNx),level)+Au(cir(i,lNx),level)*z(cir(i+1,lNx),level)) &
      & /Az(cir(i,lNx),level)
    end do

  end subroutine smooth
end module solve