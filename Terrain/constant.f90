module constant

  integer, parameter :: l = 3                !グリッドの深さ
  integer, parameter :: Nx = 1024
  integer, parameter :: ntmax = 3000        !時間ステップ
  integer, parameter :: nu1 = 2, nu2 = 1    !マルチグリッドサイクル内のsmooth回数
  real(8), parameter :: X = 2.d6            !領域サイズ
  real(8), parameter :: u_upstream = 0.5d0
  real(8), parameter :: v_upstream = 0.d0
  real(8), parameter :: g = 9.81d0           !重力定数
  real(8), parameter :: Cz = 80.d0          !Chezy 摩擦係数
  real(8), parameter :: pi = 4*atan(1.d0)    !円周率
  real(8), parameter :: f0 = 4*pi/86400      !コリオリパラメータf0
  ! real(8), parameter :: dt = 60.d0*4        !時間間隔
  real(8), parameter :: dt = 90.d0          !時間間隔
  real(8), parameter :: dtau = dt/10.d0      !移流計算用小時間間隔

  real(8), parameter :: f=f0   !コリオリパラメータ
  real(8), parameter :: dx=X/Nx  !格子間隔
end module constant