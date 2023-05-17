module constant

  integer, parameter :: num_levels = 6                !グリッドの深さ
  integer, parameter :: Nx = 512
  integer, parameter :: ntmax = 15000        !時間ステップ
  integer, parameter :: nu1 = 2, nu2 = 20    !マルチグリッドサイクル内のsmooth回数
  real(8), parameter :: X = 1.d7            !領域サイズ
  real(8), parameter :: g = 9.81d0           !重力定数
  real(8), parameter :: nu = 1.d-5/86400     !粘性係数
  ! real(8), parameter :: Cz = 80.d0          !Chezy 摩擦係数
  real(8), parameter :: pi = 4*atan(1.d0)    !円周率
  real(8), parameter :: f0 = 4*pi/86400      !コリオリパラメータf0
  real(8), parameter :: dt = 86400.d-4          !時間間隔
  real(8), parameter :: dtau = dt/10.d0      !移流計算用小時間間隔

  ! stochastic force用定数
  real(8), parameter :: mu = 1.d11/(86400**1.5)

  ! real(8), parameter :: f=f0   !コリオリパラメータ
  real(8), parameter :: dx=X/Nx  !格子間隔

  real(8), parameter :: XForce_const = 1.d-5
  ! real(8), parameter :: XForce_const = 0.d0
end module constant