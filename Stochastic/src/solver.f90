module solve
  use structure
  use calc_variables_mod
  use transfer_mod
  implicit none

contains
  recursive subroutine MGCYC(k)
    implicit none

    integer, intent(in) :: k

    integer :: nt, i

    real(8) :: tmp_z(Nx)

    real(8) :: ldx
    integer :: lNx

    call calc_level(k,ldx,lNx)

    !Presmoothing
    do nt = 1, nu1
      call smooth(k)
    end do

    !Coarse grid correction
    !Compute the defect
    do i = 1, lNx
      residual(i,k) = b(cir(i,lNx),k) + Au(cir(i-1,lNx),k)*z(cir(i-1,lNx),k) &
      & + Au(cir(i,lNx),k)*z(cir(i+1,lNx),k) - Az(cir(i,lNx),k)*z(cir(i,lNx),k)
    end do

    !Restrict the defect
    call Prolongation(k)
    call Prolongation_defect(k)
    tmp_z(:) = z(:,k-1)
    do i = 1, lNx
      b(i,k-1) = residual(i,k-1) - Au(cir(i-1,lNx),k-1)*z(cir(i-1,lNx),k-1) &
      & - Au(cir(i,lNx),k-1)*z(cir(i+1,lNx),k-1) + Az(cir(i,lNx),k-1)*z(cir(i,lNx),k-1)
    enddo

    !Compute an approximate solution v of the defect equation on k-1
    if(k==2) then
      do nt = 1, 1
        call smooth(k-1)
      end do
    else
      call MGCYC(k-1)
    end if

    residual(:,k-1) = z(:,k-1) - tmp_z(:)

    !Interpolate the correction
    call Interpolation_defect(k-1)

    !Compute the corrected approximation on k
    z(:,k) = z(:,k) + residual(:,k)

    !Postsmoothing
    do nt = 1, nu2
      call smooth(k)
    end do

  end subroutine MGCYC

  recursive subroutine FMG(max_k, k)
    implicit none

    integer, intent(in) :: max_k, k

    integer :: nt, i

    real(8) :: tmp_z(Nx)

    ! real(8) :: ldx
    ! integer :: lNx

    ! Now do a multigrid V-cycle or W-cycle on the fine grid
    do i  = 1, 20
      call MGCYC(k)
    enddo

    fmgz(:,k) = z(:,k)

    if (k < max_k) then
      call Interpolation(k)  ! Hyper Interpolate from coarse grid (z)
      call calc_coef(k+1)
      call FMG(max_k, k+1)  ! Recursive call
    end if

  end subroutine FMG


  subroutine smooth(level)
    implicit none
    
    integer, intent(in) :: level

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