module boundary_mod
	implicit none
	
contains

	! subroutine boundary(z,u,v,f,g,dt,dx,dy,Nx,Ny)
	! 	implicit none

	! 	integer, intent(in) :: Nx, Ny
	! 	real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), f(0:Ny+1), g, dt, dx, dy
	! 	real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

	! 	integer :: i, j
	! 	real(8) :: Fu, Fv

	! 	do i = 1, Nx
	! 		Fu = calc_Fu(i,1,u,v,dt,dx,dy,Nx,Ny)
	! 		z(i,0) = z(i,1) + f(1)*dy/g*Fu
	! 		Fu = calc_Fu(i,Ny,u,v,dt,dx,dy,Nx,Ny)
	! 		z(i,Ny+1) = z(i,Ny) + f(Ny)*dy/g*Fu
	! 	end do

	! 	do j = 1, Ny
	! 		Fv = calc_Fv(1,j,u,v,dt,dx,dy,Nx,Ny)
	! 		z(0,j) = z(1,j) - z_frac(f(j:j+1))*dx/g*Fv
	! 		Fv = calc_Fv(Nx,j,u,v,dt,dx,dy,Nx,Ny)
	! 		z(Nx+1,j) = z(Nx,j) - z_frac(f(j:j+1))*dx/g*Fv
	! 	end do
		
	! end subroutine boundary

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

		! u(0,:) = 0.d0
		! u(Nx,:) = 0.d0
		u(:,0) = -u(:,1)
		u(:,Ny+1) = -u(:,Ny)

	end subroutine boundary_u

	subroutine boundary_v(v,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(inout) :: v(0:Nx+1,0:Ny)

		! v(:,0) = 0.d0
		! v(:,Ny) = 0.d0
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