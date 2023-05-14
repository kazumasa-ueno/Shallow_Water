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

    real(8) :: tmp_z(Nx)

    real(8) :: ldx
    integer :: lNx

    call calc_level(k,ldx,lNx)

    !Presmoothing
    do nt = 1, nu1
      call smooth(k,z,Au,Az,b)
    end do

    !Coarse grid correction
    !Compute the defect
    tmp_z(:) = residual(:,k)
    do i = 1, lNx
      residual(i,k) = b(cir(i,lNx),k) + Au(cir(i-1,lNx),k)*z(cir(i-1,lNx),k) &
      & + Au(cir(i,lNx),k)*z(cir(i+1,lNx),k) - Az(cir(i,lNx),k)*z(cir(i,lNx),k)
    end do

    !Restrict the defect
    call Prolongation(k,u,z,h)
    call Prolongation_defect(k,residual)

    !Compute an approximate solution v of the defect equation on k-1
    z(:,k-1) = z(:,k-1) + residual(:,k-1)
    if(k==2) then
      do nt = 1, 1
        call smooth(k-1,z,Au,Az,b)
      end do
    else
      call MGCYC(k-1,u,z,h,Au,Az,b,residual,cyc,times)
    end if

    !Interpolate the correction
    call Interpolation_defect(k-1,z)

    !Compute the corrected approximation on k
    ! z(:,k) = z(:,k) + residual(:,k)

    !Postsmoothing
    do nt = 1, nu2
      call smooth(k,z,Au,Az,b)
    end do

  end subroutine MGCYC

  recursive subroutine FMG(max_k, k, u, z, h, Au, Az, b, residual, cyc, times)
    implicit none

    integer, intent(in) :: max_k, k
    real(8), intent(inout) :: Au(:,:), Az(:,:), b(:,:), u(:,:), h(:,:), residual(:,:)
    real(8), intent(inout) :: z(:,:)
    integer, intent(in) :: cyc, times

    integer :: nt, i

    real(8) :: tmp_z(Nx)

    ! real(8) :: ldx
    ! integer :: lNx

    ! Now do a multigrid V-cycle or W-cycle on the fine grid
    do i  = 1, 20
      call MGCYC(k, u, z, h, Au, Az, b, residual, cyc, times)
    enddo

    if (k < max_k) then
      call Interpolation_defect(k, z)  ! Interpolate from coarse grid
      call FMG(max_k, k+1, u, z, h, Au, Az, b, residual, cyc, times)  ! Recursive call
    end if

  end subroutine FMG


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