module boundary_mod
	implicit none
	
contains

	subroutine boundary(z,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

		z(0,:)    = z(1,:)
		z(Nx+1,:) = z(Nx,:)
		z(:,0)    = z(:,1)
		z(:,Ny+1) = z(:,Ny+1)
		
	end subroutine boundary

	subroutine boundary_defect(d,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: d(0:Nx+1,0:Ny+1)

		d(0,:)    = d(1,:)
		d(Nx+1,:) = d(Nx,:)
		d(:,0)    = d(:,1)
		d(:,Ny+1) = d(:,Ny)

	end subroutine boundary

	subroutine boundary_u(u,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: u(-1:Nx+1,0:Ny+1)

		u(-1,:) = u(1,:)
		u(Nx+1,:) = u(Nx-1,:)

	end subroutine boundary_u

	subroutine boundary_v(v,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: v(0:Nx+1,-1:Ny+1)

		v(:,-1) = v(:,1)
		v(:,Ny+1) = v(:,Ny-1)
		
	end subroutine boundary_v

end module boundary