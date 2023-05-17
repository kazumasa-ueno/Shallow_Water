!********************************************
! 境界条件を計算
!********************************************

module boundary_mod
  use constant
  implicit none
  
contains

  integer function cir(i,Nx)
    implicit none
    
    integer, intent(in) :: i, Nx

    cir = mod(i,Nx)
    if(cir==0) cir=Nx

  end function cir

end module boundary_mod