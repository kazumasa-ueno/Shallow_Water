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

	real(8) :: time, flop, elapsed_time

	real(8) :: f(ny) !!colioli parameter [s^-1]
	real(8) :: u(nx,ny,3) !!zonal flow [m/s]
	real(8) :: v(nx,ny,3) !!meridional flow [m/s]
	real(8) :: h(nx,ny,3) !! sea level deviations [m]
	real(8) :: dx(ny), dy !! [m]

	time = 0.d0
  flop = 0.d0 
  elapsed_time = 0.d0

	! open(21, file = 'output.txt', form = 'formatted')
	! open(22, file = 'output_u.txt', form = 'formatted')
	! open(23, file = 'output_v.txt', form = 'formatted')

	!initial condition
	u(:,:,2) = 0.d0
	u(:,:,1) = 0.d0
	v(:,:,2) = 0.d0
	v(:,:,1) = 0.d0
	do j = 1, ny
		do i = 1, nx
			h(i,j,2) = 10*exp(-((i-nx/2)**2+(j-ny/2)**2)/2/16.d0**2) !!Gaussian
		end do
	end do
	h(:,:,1) = h(:,:,2)
	! write(21,*) h
	! write(22,*) u
	! write(23,*) v
	
	!boundary condition
	u(1,:,2) = 0.d0
	v(:,1,2) = 0.d0
	u(nx,:,2) = 0.d0
	v(:,ny,2) = 0.d0

	!calculating colioli parameter, dx and dy
	do j=1, ny
		f(j) = f0 * sin((ymin + (j-1)*dy_deg)*pi/180.d0)
		dx(j) = pi*earth_R*cos((ymin + (j-1)*dy_deg)*pi/180.d0) * dx_deg/180
	end do
	dy = pi*earth_R*dy_deg/180

!$acc data copy(u,v,h)
	call start_timer()

	!each step
	do icnt = 1, nt
		! write(*,*) 'nt = ', icnt
		call timesteps(nx, ny, dx, dy, dt, g, A, Height, f, u, v, h)

		call update(u,v,h)

		! if(mod(icnt,120) == 0 .and. icnt<=120*100) then !each 1 day
		! 	! write(21,*) icnt/1200*1 
		! 	write(*,*) 'day = ', icnt/120
		! 	write(21,*) h
		! 	write(22,*) u
		! 	write(23,*) v
		! else if(mod(icnt,6000) == 0) then !each 50days
		! 	write(*,*) 'day = ', icnt/120
		! 	write(21,*) h
		! 	write(22,*) u
		! 	write(23,*) v
		! end if
	end do

	elapsed_time = get_elapsed_time();
!$acc end data

	write(*,*) "succeed"
    
  ! write(*, "(A7,F8.3,A6)"), "Time = ",elapsed_time," [sec]"
  ! write(*, "(A13,F7.2,A9)"), "Performance= ",flop/elapsed_time*1.0e-09," [GFlops]"

	stop
end program main