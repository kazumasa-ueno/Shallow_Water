!********************************************
! Stochastic forcingに対する応答
! 浅水波方程式をFMGで解くプログラム
! Semi-implicit Semi-Lagtange法を使用して離散化
!********************************************

program main
  use constant
  use boundary_mod
  use calc_variables_mod
  use transfer_mod
  use initialization
  use solve
  use,intrinsic :: iso_fortran_env
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer(int32) :: time_begin_c,time_end_c, CountPerSec, CountMax !時間測定用
  
  real(8) :: h(Nx,num_levels)  !基準面からの水深
  real(8) :: u(Nx,num_levels), z(Nx,num_levels), XForce(Nx,num_levels)
  real(8) :: Au(Nx,num_levels), Az(Nx,num_levels), b(Nx,num_levels) !係数
  real(8) :: residual(Nx,num_levels)
  real(8) :: Res, difference  !Resは残差のl2ノルム、differenceは前の時間との残差  
  integer :: times, cyc !時間ループ用と収束までの繰り返し用
  integer :: ios !ファイル開く用

  !for debug
  real(8) :: Prev(Nx), tmp(Nx) !前の値を格納しておくための配列

  integer :: i, l

  ! open(unit=10, file="./output/u.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/u.txt"
  open(unit=12, file="./output/z.txt", iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file ./output/z.txt"
  

  !u,z,hの初期化
  call initialize(u,z,h,XForce)
  call init_coef(Au,Az,b,residual)
  
  !メインのループ
  do times = 0, ntmax-1
    !時間計測スタート
    ! call system_clock(time_begin_c, CountPerSec, CountMax)
    
    !係数計算
    call calc_coef(u,z,h,Au,Az,b,XForce)
    
    difference = 100 !大きな値にセット
    Res = 100 !同上
    cyc = 0
    
    !収束するまで繰り返し
    ! do while(Res>1.e-19)
    !   cyc = cyc + 1
    do cyc = 1, 50
      Prev(:) = z(:,num_levels)
      
      !zの計算
      ! call MGCYC(num_levels,u,z,h,Au,Az,b,residual,cyc,times)
      ! call FMG(num_levels, 2, u, z, h, Au, Az, b, residual, cyc, times)
      call smooth(num_levels,z,Au,Az,b)
      
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
      call calc_u(l,u,z,XForce)
    enddo
    do l = 1, num_levels
      call calc_XForce(l,z,h,XForce)
    enddo
    
    !時間計測終わり
    ! call system_clock(time_end_c)
    ! print *,time_begin_c,time_end_c, CountPerSec,CountMax
    ! write(*,*) 'nt = ', times, real(time_end_c - time_begin_c)/CountPerSec,"sec"
    
    !格子中心での値を記録
    if(mod(times,100)==0) then
      write(12,*) z(:,num_levels)
    endif
  end do
  
  stop
end program main
