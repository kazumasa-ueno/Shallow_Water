program main
	use timestep
	use misc
	implicit none

	real(8), parameter :: pi = 4*atan(1.d0)
	real(8), parameter :: g = 0.02d0 ![m/s^2]
	real(8), parameter :: Height = 750.d0 ![m]
	real(8), parameter :: A = 1.d3 ![m^2/s]
	real(8), parameter :: f0 = 4*pi/86400 !!colioli parameter0 [s^-1]
	real(8), parameter :: earth_R = 6400d3 !earth radius [m]

	real(8), parameter :: dt = 12*60.d0, dx_deg = 0.125, dy_deg = 0.125, xmin = 120.d0, xmax = 180.d0, ymin = 20.d0, ymax = 50.d0
	integer, parameter :: nt = 120*1100 !!3years
	integer :: icnt, i, j
	integer, parameter :: nx = int((xmax-xmin)/dx_deg)+1
	integer, parameter :: ny = int((ymax-ymin)/dy_deg)+1
	real(8) :: term1, term2, term3

	! real(8) :: time, flop, elapsed_time

	real(8) :: f(ny) !!colioli parameter [s^-1]
	real(8) :: u(nx,ny,-1:1) !!zonal flow [m/s]
	real(8) :: v(nx,ny,-1:1) !!meridional flow [m/s]
	real(8) :: h(nx,ny,-1:1) !! sea level deviations [m]
	real(8) :: dx(ny), dy !! [m]

	open(21, file = 'output.txt', form = 'formatted')
	open(22, file = 'output_u.txt', form = 'formatted')
	open(23, file = 'output_v.txt', form = 'formatted')

	!initial condition
	u(:,:,0) = 0.d0
	u(:,:,-1) = 0.d0
	v(:,:,0) = 0.d0
	v(:,:,-1) = 0.d0
	do j = 1, ny
		do i = 1, nx
			h(i,j,0) = 10*exp(-((i-nx/2)**2+(j-ny/2)**2)/2/16.d0**2) !!Gaussian
		end do
	end do
	h(:,:,-1) = h(:,:,0)
	write(21,*) h
	write(22,*) u
	write(23,*) v
	
	!boundary condition
	u(1,:,0) = 0.d0
	v(:,1,0) = 0.d0
	u(nx,:,0) = 0.d0
	v(:,ny,0) = 0.d0

	!calculating colioli parameter, dx and dy
	do j=1, ny
		f(j) = f0 * sin((ymin + (j-1)*dy_deg)*pi/180.d0)
		dx(j) = pi*earth_R*cos((ymin + (j-1)*dy_deg)*pi/180.d0) * dx_deg/180
	end do
	dy = pi*earth_R*dy_deg/180

	!each step
	do icnt = 1, nt
		! write(*,*) 'nt = ', icnt
		do j = 2, ny-1
			do i = 2, nx-1
				!calculating u
				term1 = f(j) * (v(i-1,j,0)+v(i-1,j+1,0)+v(i,j,0)+v(i,j+1,0))/4 !!後で重み付き平均に変更
				term2 = -g * (h(i,j,0)-h(i-1,j,0))/dx(j)
				term3 = A * ((u(i+1,j,-1)-2*u(i,j,-1)+u(i-1,j,-1))/(dx(j))**2 + (u(i,j+1,-1)-2*u(i,j,-1)+u(i,j-1,-1))/(dy)**2)
				u(i,j,1) = u(i,j,-1) + 2*dt*(term1+term2+term3)
				! u(i,j,1) = u(i,j,-1) + 2*dt*(term2+term3)

				!calculating v
				term1 = -f(j) * (u(i,j-1,0)+u(i+1,j-1,0)+u(i,j,0)+u(i+1,j,0))/4 !!後で重み付き平均に変更
				term2 = -g * (h(i,j,0)-h(i,j-1,0))/dy
				term3 = A * ((v(i+1,j,-1)-2*v(i,j,-1)+v(i-1,j,-1))/(dx(j))**2 + (v(i,j+1,-1)-2*v(i,j,-1)+v(i,j-1,-1))/(dy)**2)
				v(i,j,1) = v(i,j,-1) + 2*dt*(term1+term2+term3)
				! v(i,j,1) = v(i,j,-1) + 2*dt*(term2+term3)

				!calculating h
				h(i,j,1) = h(i,j,-1) + 2*dt*(-Height*((u(i+1,j,0)-u(i,j,0))/dx(j) + (v(i,j+1,0)-v(i,j,0))/dy))
			end do
		end do
		! u(1,:,1) = 0.d0
		! u(2,:,1) = 0.d0
		! u(nx,:,1) = 0.d0
		! u(nx-1,:,1) = 0.d0
		! v(:,1,1) = 0.d0
		! v(:,2,1) = 0.d0
		! v(:,ny,1) = 0.d0
		! v(:,ny-1,1) = 0.d0
		u(:,:,-1) = u(:,:,0)
		u(:,:,0) = u(:,:,1)
		v(:,:,-1) = v(:,:,0)
		v(:,:,0) = v(:,:,1)
		h(:,:,-1) = h(:,:,0)
		h(:,:,0) = h(:,:,1)
		if(mod(icnt,120) == 0 .and. icnt<=120*100) then !each 1 day
			! write(21,*) icnt/1200*1 
			write(*,*) 'day = ', icnt/120
			write(21,*) h
			write(22,*) u
			write(23,*) v
		else if(mod(icnt,6000) == 0) then !each 50days
			write(*,*) 'day = ', icnt/120
			write(21,*) h
			write(22,*) u
			write(23,*) v
		end if
	end do

	stop
end program main