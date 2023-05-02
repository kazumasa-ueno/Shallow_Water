!********************************************
!  中央部にアノマリーを置いたときのロスビー波や重力波の伝播
! 浅水波方程式をマルチグリッド法で解くプログラム
! Semi-implicit Semi-Lagtange法を使用して離散化
!********************************************

program main
  use constant
  use boundary_mod
  use calc_variables_mod
  use transfer_mod
  use initialization
  use,intrinsic :: iso_fortran_env
  implicit none

  integer(int32) :: time_begin_c,time_end_c, CountPerSec, CountMax !時間測定用
  
  real(8) :: h(0:Nx+1)  !基準面からの水深
  real(8) :: u(0:Nx), v(0:Nx+1), z(0:Nx+1), gamma(0:Nx+1), u_b(0:Nx), v_b(0:Nx+1) !u_b, v_bは移流計算実行用の一時格納配列
  real(8) :: Au(0:Nx), Az(1:Nx), b(1:Nx) !係数
  real(8) :: Fu(Nx), Fv(Nx)
  real(8) :: Res, difference  !Resは残差のl2ノルム、differenceは前の時間との残差  
  integer :: times, cyc !時間ループ用と収束までの繰り返し用
  integer :: ios !ファイル開く用
  integer :: i !空間ループ用
  real(8) :: Fu_tmp, Fv_tmp

  !for debug
  real(8) :: Prev(0:Nx+1), tmp(Nx+2) !前の値を格納しておくための配列

  ! open(unit=10, file="./output/u.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/u.txt"
  ! open(unit=11, file="./output/v.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/v.txt"
  open(unit=12, file="./output/z.txt", iostat=ios, status="replace", action="write")
  if ( ios /= 0 ) stop "Error opening file ./output/z.txt"
  ! open(unit=20, file="./output/z1.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
  ! open(unit=21, file="./output/z2.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/z2.txt"
  ! open(unit=22, file="./output/z3.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/z3.txt"
  ! open(unit=23, file="./output/z4.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
  ! open(unit=24, file="./output/z5.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/z2.txt"
  ! open(unit=25, file="./output/z6.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/z3.txt"
  ! open(unit=26, file="./output/z7.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
  ! open(unit=27, file="./output/z8.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
  ! open(unit=30, file="./output/res5.txt", iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file ./output/res.txt"
  

  !u,v,z,gamma,hの初期化
  call initialize(u,v,z,gamma,h)

  !メインのループ
  do times = 0, ntmax-1
    !時間計測スタート
    call system_clock(time_begin_c, CountPerSec, CountMax)
    
    !係数計算
    call calc_Au(z,gamma,h,g,dt,dx,Nx,Au)
    call calc_Az(Au,Nx,Az)
    call calc_b(u,v,z,gamma,h,f,dt,dx,dtau,Nx,b)
    
    difference = 100 !大きな値にセット
    Res = 100 !同上
    cyc = 0
    
    !収束するまで繰り返し
    ! do while(Res>1.e-17)
    !   cyc = cyc + 1
    do cyc = 1, 500
      Prev(:) = z(:)
      
      !zの計算
      call MGCYC(l,z,Au,Az,b,Nx,Nx/2,Res,z(1:Nx),cyc,times)
      ! call smooth(z,Au,Az,b,Nx)
      Fu_tmp = calc_Fuu(0,u,dt,dx,dtau,Nx)
      Fv_tmp = calc_Fvu(0,u,v,dt,dx,dtau,Nx)
      call boundary(z,Nx,dt,dx,gamma,Fu_tmp,Fv_tmp,u)
      
      tmp(:) = reshape(Prev(:) - z(:),(/(Nx+2)/))
      difference = dot_product(tmp,tmp)
      
      call calc_res(z,Au,Az,b,Nx,Res)
      
      ! if(times==4 .and. cyc<51) then
      !   write(30,*) Res
      ! end if
      ! write(*,*) 'cyc = ', cyc, Res, difference
      
    end do
    ! call boundary(z,Nx,Ny)
    
    write(*,*) 'times = ', times, sum(z)
    u_b(:) = u(:)
    v_b(:) = v(:)
    call calc_u(u,v_b,z,gamma,times,Nx)
    call calc_v(u_b,v,gamma,Nx)
    call calc_gamma(u,v,z,h,gamma,Nx)
    

    !時間計測終わり
    call system_clock(time_end_c)
    ! print *,time_begin_c,time_end_c, CountPerSec,CountMax
    ! write(*,*) 'nt = ', times, real(time_end_c - time_begin_c)/CountPerSec,"sec"
    
    !格子中心での値を記録
    if(mod(times,10)==0) then
      write(12,*) z(0:Nx)
    endif
  end do

  stop
contains

  recursive subroutine MGCYC(k,z,Au,Az,b,Nx,Nxc,Res,origin_zf,cyc,times)
    implicit none

    integer, intent(in) :: k, Nx, Nxc
    real(8), intent(in) :: Au(0:Nx), Az(1:Nx), b(1:Nx)
    real(8), intent(inout) :: z(0:Nx+1)
    real(8), intent(out) :: Res
    real(8), intent(in) :: origin_zf(1:Nx)
    integer, intent(in) :: cyc, times

    integer :: nt, i, ntmax
    real(8) :: df(1:Nx), dc(1:Nxc), wf(0:Nx+1), wc(0:Nxc+1) !defects and errors etc.
    real(8) :: Auc(0:Nxc), Azc(1:Nxc)
    real(8) :: origin_zc(1:Nxc)


    !Presmoothing
    do nt = 1, nu1
      call smooth(z,Au,Az,b,Nx)
    end do

    !Coarse grid correction
    !Compute the defect
    do i = 1, Nx
      df(i) = b(i) + Au(i-1)*z(i-1) + Au(i)*z(i+1) - Az(i)*z(i)
    end do

    !Restrice the defect
    dc(:) = 0.d0
    call Prolongation(Au,Az,Auc,Azc,Nxc)
    call Prolongation_defect(df,dc,Nxc)
    call Prolongation_defect(origin_zf,origin_zc,Nxc)

    !Compute an approximate solution v of the defect equation on k-1
    wc(:) = 0.d0
    if(k==1) then
      do nt = 1, 1
        call smooth(wc,Auc,Azc,dc,Nxc)
      end do
    else
      call MGCYC(k-1,wc,Auc,Azc,dc,Nxc,Nxc/2,Res,origin_zc,cyc,times)
    end if

    !Interpolate the correction
    wf(:) = 0.d0
    call Interpolation_defect(wc,wf,Nx)

    !Compute the corrected approximation on k
    z(:) = z(:) + wf(:)

    !Postsmoothing
    do nt = 1, nu2
      call smooth(z,Au,Az,b,Nx)
    end do


    ! if(mod(times,10)==0 .and. cyc==20) then
    ! ! if(cyc==350) then
    !   select case(k)
    !   case(1)
    !     write(20,*) origin_zf(:) + z(1:Nx)
    !   case(2)
    !     write(21,*) origin_zf(:) + z(1:Nx)
    !   case(3)
    !     write(22,*) origin_zf(:) + z(1:Nx)
    !   case(4)
    !     write(23,*) origin_zf(:) + z(1:Nx)
    !   case(5)
    !     write(24,*) origin_zf(:) + z(1:Nx)
    !   case(6)
    !     write(25,*) origin_zf(:) + z(1:Nx)
    !   case(7)
    !     write(26,*) origin_zf(:) + z(1:Nx)
    !   case(8)
    !     write(27,*) origin_zf(:) + z(1:Nx)
    !   case default
    !   end select
    ! end if

  end subroutine MGCYC

  subroutine smooth(z,Au,Az,b,Nx)
    implicit none
    
    integer, intent(in) :: Nx
    real(8), intent(in) :: Au(0:Nx), Az(1:Nx), b(1:Nx)
    real(8), intent(inout) :: z(0:Nx+1)

    integer :: i

    do i = 1, Nx
    ! do i = Nx, 1, -1
      z(i) = (b(i)+Au(i-1)*z(i-1)+Au(i)*z(i+1))/Az(i)
    end do

  end subroutine smooth

end program main
