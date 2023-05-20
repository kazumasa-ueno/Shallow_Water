!********************************************
! Stochastic forcingに対する応答
! 浅水波方程式をFMGで解くプログラム
! Semi-implicit Semi-Lagtange法を使用して離散化
!********************************************

program main
  use constant
  use structure
  use boundary_mod
  use calc_variables_mod
  use transfer_mod
  use initialization
  use solve
  use output
  use,intrinsic :: iso_fortran_env
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer(int32) :: time_begin_c,time_end_c, CountPerSec, CountMax !時間測定用

  real(8) :: Res, difference  !Resは残差のl2ノルム、differenceは前の時間との残差  
  integer :: times, cyc !時間ループ用と収束までの繰り返し用
  integer :: ios !ファイル開く用

  !for debug
  real(8) :: Prev(Nx), tmp(Nx) !前の値を格納しておくための配列

  integer :: i, l

  !u,z,hの初期化
  call initialize
  call init_coef
  
  !メインのループ
  do times = 0, ntmax-1
    !時間計測スタート
    ! call system_clock(time_begin_c, CountPerSec, CountMax)
    
    !係数計算
    call calc_coef
    
    difference = 100 !大きな値にセット
    Res = 100 !同上
    cyc = 0
    
    !収束するまで繰り返し
    ! do while(Res>1.e-19)
    !   cyc = cyc + 1
    do cyc = 1, 10
      Prev(:) = z(:,num_levels)
      
      !zの計算
      ! call MGCYC(num_levels)
      call FMG(num_levels, 3)
      ! call smooth(num_levels)
      
      tmp(:) = reshape(Prev(:) - z(:,1),(/(Nx)/))
      difference = dot_product(tmp,tmp)      
      call calc_res(z,Au,Az,b,Res)      
      ! if(times==4 .and. cyc<51) then
      !   write(30,*) Res
      ! end if
      ! write(*,*) 'cyc = ', cyc, Res, difference
      
    end do
    
    write(*,*) 'times = ', times, sum(z(:,num_levels)), Res
    do l = 1, num_levels
      call calc_u(l)
    enddo
    do l = 1, num_levels
      call calc_XForce(l)
    enddo
    
    !時間計測終わり
    ! call system_clock(time_end_c)
    ! print *,time_begin_c,time_end_c, CountPerSec,CountMax
    ! write(*,*) 'nt = ', times, real(time_end_c - time_begin_c)/CountPerSec,"sec"

    do l = 3, num_levels
      z(:,l) = fmgz(:,l)
    enddo
    
    !格子中心での値を記録
    if(mod(times,100)==0) then
      call outucd(times)
    endif
  end do
  
  stop
end program main
