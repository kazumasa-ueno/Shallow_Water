module constant

  integer, parameter :: num_levels = 2                !グリッドの深さ
  integer, parameter :: Nx = 128
  integer, parameter :: ntmax = 9000        !時間ステップ
  integer, parameter :: nu1 = 2, nu2 = 1    !マルチグリッドサイクル内のsmooth回数
  real(8), parameter :: X = 4.d6            !領域サイズ
  real(8), parameter :: g = 9.81d0           !重力定数
  real(8), parameter :: nu = 1.d11/86400     !粘性係数
  ! real(8), parameter :: Cz = 80.d0          !Chezy 摩擦係数
  real(8), parameter :: pi = 4*atan(1.d0)    !円周率
  real(8), parameter :: f0 = 4*pi/86400      !コリオリパラメータf0
  real(8), parameter :: dt = 86400.d-4          !時間間隔
  real(8), parameter :: dtau = dt/10.d0      !移流計算用小時間間隔

  ! real(8), parameter :: f=f0   !コリオリパラメータ
  real(8), parameter :: dx=X/Nx  !格子間隔

  real(8), parameter :: XForce = 1.d-5
end module constant