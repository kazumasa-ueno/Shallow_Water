program adjustment
	implicit none
	
	real(8), parameter :: pi = 4*atan(1.d0)
	real(8), parameter :: g = 0.02d0 ![m/s^2]
	real(8), parameter :: Height = 750.d0 ![m]
	real(8), parameter :: A = 1.d3 ![m^2/s]
	real(8), parameter :: f0 = 4*pi/86400 !!colioli parameter0 [s^-1]
	real(8), parameter :: earth_R = 6400d3 !earth radius [m]

	real(8), parameter :: dtime = 12*60.d0, dx_deg = 0.125, dy_deg = 0.125, xmin = 120.d0, xmax = 180.d0, ymin = 20.d0, ymax = 50.d0
	integer, parameter :: ntmax = 120*1100 !!3years
	integer :: ntime, nt, i, j
	integer, parameter :: imax = int((xmax-xmin)/dx_deg)+1
	integer, parameter :: jmax = int((ymax-ymin)/dy_deg)+1
	real(8) :: term1, term2, term3


	real(8) :: f(jmax) !!colioli parameter [s^-1]
	real(8) :: u(imax,jmax) !!zonal flow [m/s]
	real(8) :: un(imax,jmax) !!zonal flow after one step [m/s]
	real(8) :: up(imax,jmax) !!zonal flow before one step [m/s]
	real(8) :: v(imax,jmax) !!meridional flow [m/s]
	real(8) :: vn(imax,jmax) !!meridional flow after one step [m/s]
	real(8) :: vp(imax,jmax) !!meridional flow before one step [m/s]
	real(8) :: h(imax,jmax) !! sea level deviations [m]
	real(8) :: hn(imax,jmax) !! sea level deviations after one step [m]
	real(8) :: hp(imax,jmax) !! sea level deviations before one step [m]
	real(8) :: dx(jmax), dy !! [m]

	open(21, file = 'output.txt', form = 'formatted')
	open(22, file = 'output_u.txt', form = 'formatted')
	open(23, file = 'output_v.txt', form = 'formatted')

	!initial condition
	u(:,:) = 0.d0
	up(:,:) = 0.d0
	v(:,:) = 0.d0
	vp(:,:) = 0.d0
	do j = 1, jmax
		do i = 1, imax
			h(i,j) = 10*exp(-((i-imax/2)**2+(j-jmax/2)**2)/2/16.d0**2) !!Gaussian
		end do
	end do
	hp(:,:) = h(:,:)
	write(21,*) h
	write(22,*) u
	write(23,*) v
	
	!boundary condition
	u(1,:) = 0.d0
	v(:,1) = 0.d0
	u(imax,:) = 0.d0
	v(:,jmax) = 0.d0

	!calculating colioli parameter, dx and dy
	do j=1, jmax
		f(j) = f0 * sin((ymin + (j-1)*dy_deg)*pi/180.d0)
		dx(j) = pi*earth_R*cos((ymin + (j-1)*dy_deg)*pi/180.d0) * dx_deg/180
	end do
	dy = pi*earth_R*dy_deg/180

	!each step
	do nt = 1, ntmax
		! write(*,*) 'nt = ', nt
		do j = 2, jmax-1
			do i = 2, imax-1
				!calculating u
				term1 = f(j) * (v(i-1,j)+v(i-1,j+1)+v(i,j)+v(i,j+1))/4 !!後で重み付き平均に変更
				term2 = -g * (h(i,j)-h(i-1,j))/dx(j)
				term3 = A * ((up(i+1,j)-2*up(i,j)+up(i-1,j))/(dx(j))**2 + (up(i,j+1)-2*up(i,j)+up(i,j-1))/(dy)**2)
				un(i,j) = up(i,j) + 2*dtime*(term1+term2+term3)
				! un(i,j) = up(i,j) + 2*dtime*(term2+term3)

				!calculating v
				term1 = -f(j) * (u(i,j-1)+u(i+1,j-1)+u(i,j)+u(i+1,j))/4 !!後で重み付き平均に変更
				term2 = -g * (h(i,j)-h(i,j-1))/dy
				term3 = A * ((vp(i+1,j)-2*vp(i,j)+vp(i-1,j))/(dx(j))**2 + (vp(i,j+1)-2*vp(i,j)+vp(i,j-1))/(dy)**2)
				vn(i,j) = vp(i,j) + 2*dtime*(term1+term2+term3)
				! vn(i,j) = vp(i,j) + 2*dtime*(term2+term3)

				!calculating h
				hn(i,j) = hp(i,j) + 2*dtime*(-Height*((u(i+1,j)-u(i,j))/dx(j) + (v(i,j+1)-v(i,j))/dy))
			end do
		end do
		! un(1,:) = 0.d0
		! un(2,:) = 0.d0
		! un(imax,:) = 0.d0
		! un(imax-1,:) = 0.d0
		! vn(:,1) = 0.d0
		! vn(:,2) = 0.d0
		! vn(:,jmax) = 0.d0
		! vn(:,jmax-1) = 0.d0
		up(:,:) = u(:,:)
		u(:,:) = un(:,:)
		vp(:,:) = v(:,:)
		v(:,:) = vn(:,:)
		hp(:,:) = h(:,:)
		h(:,:) = hn(:,:)
		if(mod(nt,120) == 0 .and. nt<=120*100) then !each 1 day
			! write(21,*) nt/1200*1 
			write(*,*) 'day = ', nt/120
			write(21,*) h
			write(22,*) u
			write(23,*) v
		else if(mod(nt,6000) == 0) then !each 50days
			write(*,*) 'day = ', nt/120
			write(21,*) h
			write(22,*) u
			write(23,*) v
		end if
	end do

	stop
end program adjustment