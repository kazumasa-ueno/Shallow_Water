!********************************************
! 境界条件を計算
!********************************************

module boundary_mod
	implicit none
	
contains

	subroutine boundary(z,Nx,Ny,Fu,Fv,f,g,dx,dy)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: Fu(:,:), Fv(:,:), f(0:Ny+1), g, dx, dy
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

		integer :: i, j

		z(:,0) = z(:,1)
		z(:,Ny+1) = z(:,Ny)
		z(0,:) = z(1,:)
		z(Nx+1,:) = z(Nx,:)
		! z(:,0) = z(:,1) - f(1)*dy/g*Fu(:,1)
		! z(:,Ny+1) = z(:,Ny) - f(Ny)*dy/g*Fu(:,2)
		! z(0,:) = z(1,:) + f(:)*dx/g*Fv(1,:)
		! z(Nx+1,:) = z(Nx,:) + f(:)*dx/g*Fv(2,:)
		
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

	subroutine boundary_u(u,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: u(0:Nx,0:Ny+1)

		u(:,0) = -u(:,1)
		u(:,Ny+1) = -u(:,Ny)

	end subroutine boundary_u

	subroutine boundary_v(v,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: v(0:Nx+1,0:Ny)

		v(0,:) = -v(1,:)
		v(Nx+1,:) = -v(Nx,:)
		
	end subroutine boundary_v

	subroutine boundary_gamma(gamma,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: gamma(0:Nx+1,0:Ny+1)

		gamma(:,0) = gamma(:,1)
		gamma(:,Ny+1) = gamma(:,Ny)
		gamma(0,:) = gamma(1,:)
		gamma(Nx+1,:) = gamma(Nx,:)
		
	end subroutine boundary_gamma
end module boundary_mod