program main
	use boundary_mod
	use calc_variables_mod
	use transfer_mod
	implicit none
	
	integer, parameter :: l = 8
	integer, parameter :: Nx = 2**l, Ny = 2**l
	integer, parameter :: ntmax = 100
	integer, parameter :: nu1 = 1, nu2 = 1
	real(8), parameter :: g = 9.81d0 !gravity acceleration
	real(8), parameter :: pi = 4*atan(1.d0)
	real(8), parameter :: f0 = 4*pi/86400
	real(8), parameter :: X = 1.d3, Y = 1.d3

	real(8) :: f(0:Ny+1) !corioli parameter
	real(8) :: h(0:Nx+1,0:Ny+1)
	real(8) :: u(-1:Nx+1,0:Ny+1), v(0:Nx+1,-1:Ny+1), z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1)
	real(8) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny)
	real(8) :: dt, dx, dy
	integer :: nt, cyc

	call initialize(u,v,z,gamma,h)
	f = 0.d0
	call calc_f(f,Ny,f0)
	call when_l(l,X,Y,Nx,Ny,dx,dy)

	do nt = 1, ntmax
		call calc_Au(z,gamma,h,g,dt,dx,Nx,Ny,Au)
		call calc_Av(z,gamma,h,g,dt,dy,Nx,Ny,Av)
		call calc_Az(Au,Av,Nx,Ny,Az)
		call calc_b(u,v,z,gamma,h,f,dt,dx,dy,Nx,Ny,b)
		do cyc = 1, 10
			call MGCYC(k,z,Au,Av,Az,b,nu1,nu2,Nx,Ny,Nxc,Nyc)
		end do
		call calc_u(u,v,z,f,gamma,dt,dx,g,Nx,Ny)
		call calc_v(u,v,z,f,gamma,dt,dy,g,Nx,Ny)
		! call calc_gamma()
	end do

	stop
contains

	subroutine initialize(u,v,z,gamma,h,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(out) :: u(-1:Nx+1,0:Ny+1), v(0:Nx+1,-1:Ny+1), z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), h(0:Nx+1,0:Ny+1)

		u(:,:) = 0.d0
		v(:,:) = 0.d0
		h(:,:) = 1.d3
		z(:,:) = 0.d0
		gamma(:,:) = 0.d0 !!!てきとう
		
	end subroutine initialize

	recursive subroutine MGCYC(k,z,Au,Av,Az,b,nu1,nu2,Nx,Ny,Nxc,Nyc)
		implicit none

		integer, intent(in) :: k, nu1, nu2, Nx, Ny, Nxc, Nyc
		real(8), intent(in) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny)
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

		integer :: nt, i, j
		real(8) :: df(0:Nx+1,0:Ny+1), dc(0:Nxc+1,0:Nyc+1), vf(0:Nx+1,0:Ny+1), vc(0:Nxc+1,0:Nyc+1) !defects and errors etc.
		real(8) :: Auc(0:Nxc,1:Nyc), Avc(1:Nxc,0:Nyc), Azc(1:Nxc,1:Nyc)


		!Presmoothing
		do nt = 1, nu1
			call smooth(z,Au,Av,Az,b,Nx,Ny)
			call boundary(z,Nx,Ny)
		end do

		!Coarse grid correction
		!Compute the defect
		do j = 1, Ny
			do i = 1, Nx
				df(i,j) = b(i,j) + Au(i-1,j)*z(i-1,j) + Au(i,j)*z(i+1,j) + Av(i,j-1)*z(i,j-1) + Av(i,j)*z(i,j+1) - Az(i,j)*z(i,j)
			end do
		end do
		call boundary_defect(df)

		!Restrice the defect
		dc(:,:) = 0.d0
		call Prolongation_defect(df,dc,Nxc,Nyc)
		call Prolongation(Au,Av,Az,Auc,Avc,Azc,Nxc,Nyc)

		!Compute an approximate solution v of the defect equation on k-1
		vc(:,:) = 0.d0
		if(k==2) then
			do nt = 1, 100
				call smooth(vc,Auc,Avc,Azc,dc,Nxc,Nyc)
				call boundary_defect(vc,Nxc,Nyc)
			end do
		else
			call MGCYC(k-1,vc,Auc,Avc,Azc,dc,nu1,nu2,Nxc,Nyc,Nxc/2,Nyc/2)
		end if

		!Interpolate the correction
		vf(:,:) = 0.d0
		call Interpolation_defect(vc,vf,Nx,Ny)

		!Compute the corrected approximation on k
		z(:,:) = z(:,:) + vf(:,:)

		!Postsmoothing
		do nt = 1, nu2
			call smooth(z,Au,Av,Az,b,Nx,Ny)
			call boundary(z,Nx,Ny)
		end do

	end subroutine MGCYC

	subroutine smooth(z,Au,Av,Az,b,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny)
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

		integer :: i, j

		do j = 1, Ny
			do i = 1, Nx
				z(i,j) = (b(i,j)+Au(i-1,j)*z(i-1,j)+Au(i,j)*z(i+1,j)+Av(i,j-1)*z(i,j-1)+Av(i,j)*z(i,j+1))/Az(i,j)
			end do
		end do

	end subroutine smooth

end program main