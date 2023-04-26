!********************************************
! 境界条件を計算
!********************************************

module boundary_mod
	implicit none
	
contains

	subroutine boundary(z,Nx)
		implicit none

		integer, intent(in) :: Nx
		real(8), intent(inout) :: z(0:Nx+1)

		z(0) = z(Nx)
		z(Nx+1) = z(1)
		
	end subroutine boundary

	subroutine boundary_defect(d,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: d(0:Nx+1,0:Ny+1)

		d(0,:)    = d(1,:)
		d(Nx+1,:) = d(Nx,:)
		d(:,0)    = d(:,1)
		d(:,Ny+1) = d(:,Ny)

	end subroutine boundary_defect

	subroutine boundary_u(u,Nx)
		implicit none

		integer, intent(in) :: Nx
		real(8), intent(inout) :: u(0:Nx)

		u(0) = u(Nx)

	end subroutine boundary_u

end module boundary_mod