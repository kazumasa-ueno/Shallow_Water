module boundary_mod
	implicit none
	
contains

	subroutine boundary(z,Nx,Ny,jmin,jmax)
		implicit none

		integer, intent(in) :: Nx, Ny, jmin, jmax
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

		integer :: i, j

		do j = 0, jmin-1
			z(0,j) = z(1,j)
		end do
		do j = jmax+1, Ny+1
			z(0,j) = z(1,j)
		end do

		z(:,0) = z(:,1)
		z(:,Ny+1) = z(:,Ny)
		z(Nx+1,:) = z(Nx,:)
		
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