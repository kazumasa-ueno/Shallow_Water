!********************************************
! 色々な変数を計算するためのモジュール
!	when_l: 			dx,dyを計算
!	calc_f: 			コリオリパラメータを計算
!	calc_u: 			zの値からuを計算
!	calc_v: 			zの値からvを計算
!	calc_gamma: 			zの値からgammaを計算
!	calc_Au:			Auの値を計算
!	calc_Av:			Avの値を計算
!	calc_Az:			Azの値を計算
!	calc_b:				bの値を計算
!	calc_Fu:			移流元のuの値を計算
!	calc_Fv:			移流元のvの値を計算
!	inner_u:			uの内装を計算
!	inner_v:			vの内装を計算
!	z_frac:				分数インデックスでの値を計算
!	calc_channel:			湾幅のインデックスを計算
!	channel_z:			潮位の時間変化を計算
!	channel_z_defect:		潮位の時間変化の残差
!	calc_res:			残差を計算

! calc_tau:     風応力を計算
!********************************************

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

		dy = Y / (Ny-1)
		dy_rad = dy / earth_R

		do j = 1, Ny
			! f(j) = f0 * sin((j-Ny/2)*dy_rad + pi/6.d0)
			f(j) = f0 * sin((j-Ny/2)*dy_rad)
		end do
		f(0) = f(1)
		f(Ny+1) = f(Ny)
		
		! f(:) = f0
		! f(:) = 0.d0

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
				Fu = calc_Fux(i,j,u_b,v,dt,dx,dy,dtau,Nx,Ny)
				Fv = calc_Fvx(i,j,u_b,v,dt,dx,dy,dtau,Nx,Ny)
				u(i,j) = (Fu - g*(dt/dx)*(z(i+1,j)-z(i,j)) + f(j)*dt*Fv) / (1+z_frac(gamma(i:i+1,j))*dt)
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
				Fu = calc_Fuy(i,j,u,v_b,dt,dx,dy,dtau,Nx,Ny)
				Fv = calc_Fvy(i,j,u,v_b,dt,dx,dy,dtau,Nx,Ny)
				v(i,j) = (Fv - g*(dt/dy)*(z(i,j+1)-z(i,j)) - f(j)*dt*Fu) / (1+z_frac(gamma(i,j:j+1))*dt)
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

		! do j = 1, Ny
		! 	do i = 1, Nx
		! 		gamma(i,j) = g*(z_frac(u(i-1:i,j))**2+z_frac(v(i,j-1:j))**2)**0.5d0 / Cz**2 / (z(i,j)+h(i,j))
		! 	end do
		! end do
		! call boundary_gamma(gamma,Nx,Ny)
		gamma = 0.d0
		
	end subroutine calc_gamma

	subroutine calc_Au(z,gamma,h,g,dt,dx,Nx,Ny,Au)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), h(0:Nx+1,0:Ny+1), g, dt, dx
		real(8), intent(out) :: Au(0:Nx,1:Ny)

		integer :: i, j

		do j = 1, Ny
			do i = 0, Nx
				! Au(i,j) = g * (dt/dx)**2 * (z_frac(z(i:i+1,j)) + z_frac(h(i:i+1,j))) / (1 + z_frac(gamma(i:i+1,j))*dt)
				Au(i,j) = g * (dt/dx)**2 * (z_frac(h(i:i+1,j))) / (1 + z_frac(gamma(i:i+1,j))*dt)
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
				! Av(i,j) = g * (dt/dy)**2 * (z_frac(z(i,j:j+1)) + z_frac(h(i,j:j+1))) / (1 + z_frac(gamma(i,j:j+1))*dt)
				Av(i,j) = g * (dt/dy)**2 * (z_frac(h(i,j:j+1))) / (1 + z_frac(gamma(i,j:j+1))*dt)
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
			end do
		end do

	end subroutine calc_Az

	! subroutine calc_b(u,v,z,gamma,h,f,dt,dx,dy,dtau,Nx,Ny,times,b)
	subroutine calc_b(u,v,z,gamma,h,f,dt,dx,dy,dtau,Nx,Ny,times,b,q)
		implicit none
		
		integer, intent(in) :: Nx, Ny, times
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), z(0:Nx+1,0:Ny+1), h(0:Nx+1,0:Ny+1), f(0:Ny+1), gamma(0:Nx+1,0:Ny+1)
		real(8), intent(in) :: dt, dx, dy, dtau
		real(8), intent(out) :: b(1:Nx,1:Ny), q(1:Nx,1:Ny)

		integer :: i, j
		real(8) :: Fu(4), Fv(4), taux
		
		real(8), parameter :: D = 80.d0
		real(8), parameter :: rho = 1.d3

		do j = 1, Ny
			do i = 1, Nx
				Fu(1) = calc_Fux(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fu(2) = calc_Fux(i-1,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				! Fu(3) = calc_Fux(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				! Fu(4) = calc_Fux(i,j-1,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fu(3) = calc_Fuy(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fu(4) = calc_Fuy(i,j-1,u,v,dt,dx,dy,dtau,Nx,Ny)
				! Fv(1) = calc_Fvy(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				! Fv(2) = calc_Fvy(i-1,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fv(1) = calc_Fvx(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fv(2) = calc_Fvx(i-1,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fv(3) = calc_Fvy(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
				Fv(4) = calc_Fvy(i,j-1,u,v,dt,dx,dy,dtau,Nx,Ny)
				taux  = calc_tau(times,dt,dx,dy,i,j,Nx,Ny)
				b(i,j) = z(i,j) - (dt/dx) * ( &
				& (z_frac(z(i:i+1,j))+z_frac(h(i:i+1,j))) / (1+z_frac(gamma(i:i+1,j))*dt) &
				& * (Fu(1)+z_frac(f(j:j+1))*dt*Fv(1)) &
				& - (z_frac(z(i-1:i,j))+z_frac(h(i-1:i,j))) / (1+z_frac(gamma(i-1:i,j))*dt) &
				& * (Fu(2)+z_frac(f(j:j+1))*dt*Fv(2)) ) &
				& - (dt/dy) * ( &
				& (z_frac(z(i,j:j+1))+z_frac(h(i,j:j+1))) / (1+z_frac(gamma(i,j:j+1))*dt) &
				& * (Fv(3)-f(j)*dt*Fu(3)) &
				& - (z_frac(z(i,j-1:j))+z_frac(h(i,j-1:j))) / (1+z_frac(gamma(i,j-1:j))*dt) &
				& * (Fv(4)-f(j-1)*dt*Fu(4)) ) !&
				! & - taux/rho/D*dt

				! b(i,j) = z(i,j) - (dt/dx) * ( &
				! & (z_frac(h(i:i+1,j))) / (1+z_frac(gamma(i:i+1,j))*dt) &
				! & * (Fu(1)+z_frac(f(j:j+1))*dt*Fv(1)) &
				! & - (z_frac(h(i-1:i,j))) / (1+z_frac(gamma(i-1:i,j))*dt) &
				! & * (Fu(2)+z_frac(f(j:j+1))*dt*Fv(2)) ) &
				! & - (dt/dy) * ( &
				! & (z_frac(h(i,j:j+1))) / (1+z_frac(gamma(i,j:j+1))*dt) &
				! & * (Fv(3)-f(j)*dt*Fu(3)) &
				! & - (z_frac(h(i,j-1:j))) / (1+z_frac(gamma(i,j-1:j))*dt) &
				! & * (Fv(4)-f(j-1)*dt*Fu(4)) ) !&
				! ! & - taux/rho/D*dt

				q(i,j) = Fv(1)

				! write(*,*) "^^"
				! if(i==Nx .and. j==Ny) then
				! 	write(*,*) taux/rho/D*dt, b(i,j)
				! end if
			end do
		end do

	end subroutine calc_b

	! calculate Fu(i+1/2,j)
	real(8) function calc_Fux(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny, i, j
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), dt, dx, dy, dtau

		integer :: s, smax
		real(8) :: x, y, u_s, v_s

		! if(i==0 .or. i==Nx) then
		! 	calc_Fux = 0.d0
		! else
		! 	smax = int(dt/dtau)
		! 	x = i*dx
		! 	y = j*dy
		! 	u_s = u(i,j)
		! 	v_s = (v(i,j)+v(i,j-1)+v(i+1,j)+v(i+1,j-1))*0.25d0
		! 	do s = 1, smax
		! 		x = x - dtau*u_s
		! 		y = y - dtau*v_s
		! 		call inner_u(x,y,u_s,v_s,u,v,dx,dy,Nx,Ny)
		! 	end do
		! 	calc_Fux = u_s
		! end if

		calc_Fux = u(i,j) !!!!!!!

	end function calc_Fux

	! calculate Fu(i,j+1/2)
	real(8) function calc_Fuy(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny, i, j
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), dt, dx, dy, dtau

		integer :: s, smax
		real(8) :: x, y, u_s, v_s

		if(j==0 .or. j==Ny) then
			calc_Fuy = 0.d0
		else
			smax = int(dt/dtau)
			x = i*dx
			y = j*dy
			u_s = (u(i,j)+u(i,j+1)+u(i-1,j)+u(i-1,j+1))*0.25d0
			v_s = v(i,j)
			do s = 1, smax
				x = x - dtau*u_s
				y = y - dtau*v_s
				call inner_u(x,y,u_s,v_s,u,v,dx,dy,Nx,Ny)
			end do
			calc_Fuy = u_s
		end if

		! calc_Fuy = (u(i,j)+u(i,j+1)+u(i-1,j)+u(i-1,j+1))*0.25d0 !!!!!!!!

	end function calc_Fuy

	! calculate Fv(i,j+1/2)
	real(8) function calc_Fvy(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny, i, j
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), dt, dx, dy, dtau

		integer :: s, smax
		real(8) :: x, y, u_s, v_s

		if(j==0 .or. j==Ny) then
			calc_Fvy = 0.d0
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
			calc_Fvy = v_s
		end if

		! calc_Fvy = v(i,j) !!!!!!!!!

	end function calc_Fvy

	! calculate Fv(i+1/2,j)
	real(8) function calc_Fvx(i,j,u,v,dt,dx,dy,dtau,Nx,Ny)
		implicit none
		
		integer, intent(in) :: Nx, Ny, i, j
		real(8), intent(in) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), dt, dx, dy, dtau

		integer :: s, smax
		real(8) :: x, y, u_s, v_s

		! if(i==0 .or. i==Nx) then
		! 	calc_Fvx = 0.d0
		! else
		! 	smax = int(dt/dtau)
		! 	x = i*dx
		! 	y = j*dy
		! 	u_s = u(i,j)
		! 	v_s = (v(i,j)+v(i,j-1)+v(i+1,j)+v(i+1,j-1))*0.25d0
		! 	do s = 1, smax
		! 		x = x - dtau*u_s
		! 		y = y - dtau*v_s
		! 		call inner_v(x,y,u_s,v_s,u,v,dx,dy,Nx,Ny)
		! 	end do
		! 	calc_Fvx = v_s
		! end if

		calc_Fvx = (v(i,j)+v(i,j-1)+v(i+1,j)+v(i+1,j-1))*0.25d0

	end function calc_Fvx

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
	
	subroutine calc_res(z,Au,Av,Az,b,Nx,Ny,Res)
		implicit none
		
		integer, intent(in) :: Nx, Ny
		real(8), intent(in) :: z(0:Nx+1,0:Ny+1), Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny)
		real(8), intent(out) :: Res

		integer :: i, j

		Res = 0.d0
		do j = 1, Ny
			do i = 1, Nx
				Res = Res + (b(i,j) + Au(i-1,j)*z(i-1,j) + Au(i,j)*z(i+1,j) + Av(i,j-1)*z(i,j-1) + Av(i,j)*z(i,j+1) - Az(i,j)*z(i,j))**2
			end do
		end do
		Res = (Res**0.5d0)/Nx/Ny

	end subroutine calc_res

	real(8) function calc_tau(times,dt,dx,dy,i,j,Nx,Ny)
		implicit none
		
		integer, intent(in) :: times, i, j, Nx, Ny
		real(8), intent(in) :: dt, dx, dy

		real(8), parameter :: tau_wwb = 0.02d0
		real(8), parameter :: Lx = 1000.d3
		real(8), parameter :: Ly = 1000.d3
		real(8), parameter :: t0 = 86400*7.d3
		real(8) :: x0
		
		x0 = int(Nx*5.d0/14.d0)
		calc_tau = tau_wwb * exp(-(dt/t0*times)**2-((j-Ny/2)*dy/Ly)**2-((i-x0)*dx/Lx)**2)

	end function calc_tau

end module calc_variables_mod
