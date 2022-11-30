module calc_variables_mod
	use boundary_mod
	implicit none
	
contains

	!calculate grid size and maximum grid number
	subroutine when_l(l,X,Y,Nx,Ny,dx,dy)
		implicit none
		
		integer, intent(in) :: l
		real(8), intent(in) :: X, Y
		integer, intent(in) :: Nx, Ny
		real(8), intent(out) :: dx, dy

		dx = X/Nx
		dy = Y/Ny
	
	end subroutine when_l

	! calculate corioli parameter
	subroutine calc_f(f,Ny,f0,Y)
		implicit none

		integer, intent(in)  :: Ny
		real(8), intent(in)  :: f0, Y
		real(8), intent(out) :: f(0:Ny+1)

		real(8), parameter :: earth_R = 6400d3 !earth radius [m]
		real(8), parameter :: pi = 4*atan(1.d0)
		integer :: j
		real(8) :: dy, dy_rad

		! dy = Y / (Ny-1)
		! dy_rad = dy / earth_R * pi

		! do j = 1, Ny
		! 	f(j) = f0 * sin((j-Ny/2)*dy_rad + pi/6.d0)
		! end do
		! f(0) = f(1)
		! f(Ny+1) = f(Ny)
		
		f(:) = f0

	end subroutine calc_f

	subroutine calc_u(u,v,z,f,gamma,dt,dx,dy,dtau,g,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: v(0:Nx+1,0:Ny), z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), f(0:Ny+1), dt, dx, dy, dtau, g
		real(8), intent(inout) :: u(0:Nx,0:Ny+1)

		integer :: i, j
		real(8) :: Fu, Fv, u_b(0:Nx,0:Ny+1) !u before update

		u_b(:,:) = u(:,:)

		do j = 1, Ny-1
			do i = 1, Nx-1
				Fu = calc_Fu(i,j,u_b,v,dt,dx,dy,dtau,Nx,Ny)
				Fv = calc_Fv(i,j,u_b,v,dt,dx,dy,dtau,Nx,Ny)
				u(i,j) = (Fu - g*(dt/dx)*(z(i+1,j)-z(i,j)) - f(j)*dt*Fv) / (1+z_frac(gamma(i:i+1,j))*dt)
				! u(i,j) = (u(i,j) - g*(dt/dx)*(z(i+1,j)-z(i,j)) - f(j)*dt*v(i,j)) / (1+z_frac(gamma(i:i+1,j))*dt)
			end do
		end do
		call boundary_u(u,Nx,Ny)

	end subroutine calc_u

	subroutine calc_v(u,v,z,f,gamma,dt,dx,dy,dtau,g,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: u(0:Nx,0:Ny+1), z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), f(0:Ny+1), dt, dx, dy, dtau, g
		real(8), intent(inout) :: v(0:Nx+1,0:Ny)

		integer :: i, j
		real(8) :: Fu, Fv, v_b(0:Nx+1,0:Ny) !v before update

		v_b(:,:) = v(:,:)

		do j = 1, Ny-1
			do i = 1, Nx-1
				Fu = calc_Fu(i,j,u,v_b,dt,dx,dy,dtau,Nx,Ny)
				Fv = calc_Fv(i,j,u,v_b,dt,dx,dy,dtau,Nx,Ny)
				v(i,j) = (Fv - g*(dt/dy)*(z(i,j+1)-z(i,j)) + f(j)*dt*Fu) / (1+z_frac(gamma(i,j:j+1))*dt)
				! v(i,j) = (v(i,j) - g*(dt/dy)*(z(i,j+1)-z(i,j)) + f(j)*dt*u(i,j)) / (1+z_frac(gamma(i,j:j+1))*dt)
			end do
		end do
		call boundary_v(v,Nx,Ny)
		
	end subroutine calc_v

	subroutine calc_gamma(u,v,z,h,gamma,g,Cz,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), z(0:Nx+1,0:Ny+1), h(0:Nx+1,0:Ny+1), g, Cz
		real(8), intent(out) :: gamma(0:Nx+1,0:Ny+1)

		integer :: i, j

		do j = 1, Ny
			do i = 1, Nx
				gamma(i,j) = g*(z_frac(u(i-1:i,j))**2+z_frac(v(i,j-1:j))**2)**0.5d0 / Cz**2 / (z(i,j)+h(i,j))
			end do
		end do
		call boundary_gamma(gamma,Nx,Ny)
		
	end subroutine calc_gamma

	subroutine calc_Au(z,gamma,h,g,dt,dx,Nx,Ny,Au)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), h(0:Nx+1,0:Ny+1), g, dt, dx
		real(8), intent(out) :: Au(0:Nx,1:Ny)

		integer :: i, j

		do j = 1, Ny
			do i = 0, Nx
				Au(i,j) = g * (dt/dx)**2 * (z_frac(z(i:i+1,j)) + z_frac(h(i:i+1,j))) / (1 + z_frac(gamma(i:i+1,j))*dt)
			end do
		end do

	end subroutine calc_Au

	subroutine calc_Av(z,gamma,h,g,dt,dy,Nx,Ny,Av)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), h(0:Nx+1,0:Ny+1), g, dt, dy
		real(8), intent(out) :: Av(1:Nx,0:Ny)

		integer :: i, j

		do j = 0, Ny
			do i = 1, Nx
				Av(i,j) = g * (dt/dy)**2 * (z_frac(z(i,j:j+1)) + z_frac(h(i,j:j+1))) / (1 + z_frac(gamma(i,j:j+1))*dt)
			end do
		end do

	end subroutine calc_Av

	subroutine calc_Az(Au,Av,Nx,Ny,Az)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny)
		real(8), intent(out) :: Az(1:Nx,1:Ny)

		integer :: i, j

		do j = 1, Ny
			do i = 1, Nx
				Az(i,j) = 1 + Au(i,j) + Au(i-1,j) + Av(i,j) + Av(i,j-1)
				! write(*,*) Az(i,j)
			end do
		end do

	end subroutine calc_Az

	subroutine calc_b(u,v,z,gamma,h,f,dt,dx,dy,dtau,Nx,Ny,b)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), z(0:Nx+1,0:Ny+1), h(0:Nx+1,0:Ny+1), f(0:Ny+1), gamma(0:Nx+1,0:Ny+1)
		real(8), intent(in) :: dt, dx, dy, dtau
		real(8), intent(out) :: b(1:Nx,1:Ny)

		integer :: i, j
		real(8) :: Fu(3), Fv(3)

		do j = 1, Ny
			do i = 1, Nx
				Fu(1) = calc_Fu(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fu(2) = calc_Fu(i-1,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fu(3) = calc_Fu(i,j-1,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fv(1) = calc_Fv(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fv(2) = calc_Fv(i-1,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fv(3) = calc_Fv(i,j-1,u,v,dt,dx,dy,dtau,Nx,Ny)
				b(i,j) = z(i,j) - (dt/dx) * ( &
				& (z_frac(z(i:i+1,j))+z_frac(h(i:i+1,j))) / (1+z_frac(gamma(i:i+1,j))*dt) &
				& * (Fu(1)-z_frac(f(j:j+1))*dt*Fv(1)) &
				& - (z_frac(z(i-1:i,j))+z_frac(h(i-1:i,j))) / (1+z_frac(gamma(i-1:i,j))*dt) &
				& * (Fu(2)-z_frac(f(j:j+1))*dt*Fv(2)) ) &
				& - (dt/dy) * ( &
				& (z_frac(z(i,j:j+1))+z_frac(h(i,j:j+1))) / (1+z_frac(gamma(i,j:j+1))*dt) &
				& * (Fv(1)+f(j)*dt*Fu(1)) &
				& - (z_frac(z(i,j-1:j))+z_frac(h(i,j-1:j))) / (1+z_frac(gamma(i,j-1:j))*dt) &
				& * (Fv(3)+f(j-1)*dt*Fu(3)) )
			end do
		end do

	end subroutine calc_b

	! calculate Fu(i+1/2,j)
	real(8) function calc_Fu(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny, i, j
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), dt, dx, dy, dtau

		integer :: s, smax
		real(8) :: x, y, u_s, v_s

		if(i==0 .or. i==Nx) then
			calc_Fu = 0.d0
		else
			smax = int(dt/dtau)
			x = i*dx
			y = j*dy
			u_s = u(i,j)
			v_s = (v(i,j)+v(i,j-1)+v(i+1,j)+v(i+1,j-1))*0.25d0
			do s = 1, smax
				x = x - dtau*u_s
				y = y - dtau*v_s
				call inner_u(x,y,u_s,v_s,u,v,dx,dy,Nx,Ny)
			end do
			calc_Fu = u_s
		end if

	end function calc_Fu

	! calculate Fv(i,j+1/2)
	real(8) function calc_Fv(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny, i, j
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), dt, dx, dy, dtau

		integer :: s, smax
		real(8) :: x, y, u_s, v_s

		if(j==0 .or. j==Ny) then
			calc_Fv = 0.d0
		else
			smax = int(dt/dtau)
			x = i*dx
			y = j*dy
			u_s = (u(i,j)+u(i,j+1)+u(i-1,j)+u(i-1,j+1))*0.25d0
			v_s = v(i,j)
			do s = 1, smax
				x = x - dtau*u_s
				y = y - dtau*v_s
				call inner_v(x,y,u_s,v_s,u,v,dx,dy,Nx,Ny)
			end do
			calc_Fv = v_s
		end if

	end function calc_Fv

	subroutine inner_u(x,y,u_s,v_s,u,v,dx,dy,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: x, y, u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), dx, dy
		real(8), intent(out) :: u_s, v_s
		integer :: iu, ju, iv, jv
		real(8) :: pu, qu, xv, yv, pv, qv

		iu = int(x/dx)
		ju = int(y/dy)
		pu = x/dx - iu
		qu = y/dy - ju
		u_s = (1-pu)*((1-qu)*u(iu,ju)+qu*u(iu,ju+1)) + pu*((1-qu)*u(iu+1,ju)+qu*u(iu+1,ju+1))
		xv = x/dx + 0.5d0
		yv = y/dy - 0.5d0
		iv = int(xv)
		jv = int(yv)
		pv = xv - iv
		qv = yv - jv
		v_s = (1-pv)*((1-qv)*v(iv,jv)+qv*v(iv,jv+1)) + pv*((1-qv)*v(iv+1,jv)+qv*v(iv+1,jv+1))

	end subroutine inner_u

	subroutine inner_v(x,y,u_s,v_s,u,v,dx,dy,Nx,Ny)
		implicit none

		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: x, y, u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), dx, dy
		real(8), intent(out) :: u_s, v_s
		integer :: iu, ju, iv, jv
		real(8) :: pu, qu, xu, yu, pv, qv

		xu = x/dx - 0.5d0
		yu = y/dy + 0.5d0
		iu = int(xu)
		ju = int(yu)
		pu = xu - iu
		qu = yu - ju
		u_s = (1-pu)*((1-qu)*u(iu,ju)+qu*u(iu,ju+1)) + pu*((1-qu)*u(iu+1,ju)+qu*u(iu+1,ju+1))
		iv = int(x/dx)
		jv = int(y/dy)
		pv = x/dx - iv
		qv = y/dy - jv
		v_s = (1-pv)*((1-qv)*v(iv,jv)+qv*v(iv,jv+1)) + pv*((1-qv)*v(iv+1,jv)+qv*v(iv+1,jv+1))

	end subroutine inner_v

	!calculate the value when z's index is fractional
	real(8) function z_frac(z)
		implicit none

		real(8), intent(in) :: z(2)

		z_frac = (z(1)+z(2))*0.5d0

	end function z_frac

	subroutine calc_channel(Ny,Y,jmin,jmax)
		implicit none
		
		integer, intent(in) :: Ny
		real(8), intent(in) :: Y
		integer, intent(out) :: jmin, jmax

		real(8) :: dy
		integer :: h

		dy = Y/Ny
		h = int(150.d0/dy)
		jmin = Ny/2 - h
		jmax = Ny/2 + h

	end subroutine calc_channel

	subroutine channel_z(z,time,pi,Nx,Ny,jmin,jmax)
		implicit none
		
		integer, intent(in) :: Nx, Ny, jmin, jmax
		real(8), intent(in) :: time, pi
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

		integer :: j

		do j = jmin, jmax
			z(0,j) = 0.4d0*sin(2.d0*pi*time/12.d0/3.6d3)
		end do

	end subroutine channel_z

	subroutine channel_z_defect(df,Nx,Ny,jmin,jmax)
		implicit none
		
		integer, intent(in) :: Nx, Ny, jmin, jmax
		real(8), intent(inout) :: df(0:Nx+1,0:Ny+1)

		integer :: j

		do j = jmin, jmax
			df(0,j) = 0.d0
		end do

	end subroutine channel_z_defect

end module calc_variables_mod