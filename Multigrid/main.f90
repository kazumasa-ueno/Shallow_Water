program main
	use boundary_mod
	use calc_variables_mod
	use transfer_mod
	implicit none
	
	integer, parameter :: l = 4
	integer, parameter :: Nx = 160, Ny = 80
	! integer, parameter :: Nx = 48, Ny = 24
	integer, parameter :: ntmax = 10
	integer, parameter :: nu1 = 2, nu2 = 1
	real(8), parameter :: g = 9.81d0 !gravity acceleration
	real(8), parameter :: Cz = 80.d0
	real(8), parameter :: pi = 4*atan(1.d0)
	real(8), parameter :: f0 = 4*pi/86400
	real(8), parameter :: X = 6.d3, Y = 3.d3
	real(8), parameter :: dt = 60.d0*12
	real(8), parameter :: dtau = dt/100.d0

	real(8) :: f(0:Ny+1) !corioli parameter
	real(8) :: h(0:Nx+1,0:Ny+1)
	real(8) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), u_b(0:Nx,0:Ny+1), v_b(0:Nx+1,0:Ny)
	real(8) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny)
	real(8) :: dx, dy
	real(8) :: Res, difference
	integer :: times, cyc, ios, jmin, jmax
	integer :: i, j

	!for debug
	real(8) :: Prev(0:Nx+1,0:Ny+1), tmp((Nx+2)*(Ny+2))

	open(unit=10, file="./output/u.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/u.txt"
	open(unit=11, file="./output/v.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/v.txt"
	open(unit=12, file="./output/z.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z.txt"
	open(unit=20, file="./output/z1.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	open(unit=21, file="./output/z2.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z2.txt"
	open(unit=22, file="./output/z3.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z3.txt"
	open(unit=23, file="./output/z4.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	open(unit=24, file="./output/z5.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z2.txt"
	open(unit=25, file="./output/z6.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z3.txt"
	open(unit=26, file="./output/z7.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/z1.txt"
	open(unit=30, file="./output/res3.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/res.txt"
	

	call calc_channel(Ny,Y,jmin,jmax)
	call initialize(u,v,z,gamma,h,Nx,Ny,jmin,jmax)
	! f(:) = 0.d0
	call calc_f(f,Ny,f0,Y)
	call when_l(l,X,Y,Nx,Ny,dx,dy)

	do times = 1, ntmax
		call channel_z(z,(dt-1)*times,pi,Nx,Ny,jmin,jmax) !calculate z(1,jmin:jmax)
		call calc_Au(z,gamma,h,g,dt,dx,Nx,Ny,Au)
		call calc_Av(z,gamma,h,g,dt,dy,Nx,Ny,Av)
		call calc_Az(Au,Av,Nx,Ny,Az)
		call calc_b(u,v,z,gamma,h,f,dt,dx,dy,dtau,Nx,Ny,b)
		difference = 100
		Res = 100
		cyc = 0
		do while(Res>1.e-5)
			cyc = cyc + 1
		! do cyc = 1, 100
			Prev(:,:) = z(:,:)
			! call MGCYC(l,l,z,Au,Av,Az,b,nu1,nu2,Y,Nx,Ny,Nx/2,Ny/2,Res)
			if(times==5) then
				write(30,*) Res
			end if

			call smooth(z,Au,Av,Az,b,Nx,Ny,jmin,jmax)
			call boundary(z,Nx,Ny,jmin,jmax)
			tmp(:) = reshape(Prev(:,:) - z(:,:),(/(Nx+2)*(Ny+2)/))
			difference = dot_product(tmp,tmp)
			! write(*,*) 'cyc = ', cyc, difference
			! if(times==100) then
			! 	write(12,*) z(1:Nx,1:Ny)
			! end if
			Res = 0.d0
			do j = 1, Ny
				do i = 1, Nx
					Res = Res + (b(i,j) + Au(i-1,j)*z(i-1,j) + Au(i,j)*z(i+1,j) + Av(i,j-1)*z(i,j-1) + Av(i,j)*z(i,j+1) - Az(i,j)*z(i,j))**2
				end do
			end do
			Res = Res**0.5d0
			write(*,*) 'cyc = ', cyc, Res
		end do

		write(*,*) 'times = ', times
		u_b(:,:) = u(:,:)
		v_b(:,:) = v(:,:)
		call calc_u(u,v_b,z,f,gamma,dt,dx,dy,dtau,g,Nx,Ny)
		call calc_v(u_b,v,z,f,gamma,dt,dx,dy,dtau,g,Nx,Ny)
		call calc_gamma(u,v,z,h,gamma,g,Cz,Nx,Ny)
		do j = 1, Ny
			do i = 1, Nx
				write(10,*) z_frac(u(i-1:i,j))
				write(11,*) z_frac(v(i,j-1:j))
			end do
		end do
		! write(10,*) u(1:Nx,1:Ny)
		! write(11,*) v(1:Nx,1:Ny)
		write(12,*) z(1:Nx,1:Ny)
	end do

	stop
contains

	subroutine initialize(u,v,z,gamma,h,Nx,Ny,jmin,jmax)
		implicit none

		integer, intent(in) :: Nx, Ny, jmin, jmax
		real(8), intent(out) :: u(0:Nx,0:Ny+1), v(0:Nx+1,0:Ny), z(0:Nx+1,0:Ny+1), gamma(0:Nx+1,0:Ny+1), h(0:Nx+1,0:Ny+1)

		integer :: i, j

		u(:,:) = 0.d0
		v(:,:) = 0.d0
		z(:,:) = 0.d0
		h(:,:) = 0.5d0
		h(:,jmin:jmax) = 5.d0
		! do j = 1, Ny
		! 	do i = 1, Nx
		! 		z(i,j) = 10*exp(-((i-Nx)**2+(j-Ny/2)**2)/2.d0/16.d0**2) !!Gaussian
		! 	end do
		! end do
		gamma(:,:) = 0.d0
		
	end subroutine initialize

	recursive subroutine MGCYC(l,k,z,Au,Av,Az,b,nu1,nu2,Y,Nx,Ny,Nxc,Nyc,Res)
		implicit none

		integer, intent(in) :: l, k, nu1, nu2, Nx, Ny, Nxc, Nyc
		real(8), intent(in) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny), Y
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)
		real(8), intent(out) :: Res

		integer :: nt, i, j, jmin, jmax, jminc, jmaxc, ntmax
		real(8) :: df(1:Nx,1:Ny), dc(1:Nxc,1:Nyc), wf(0:Nx+1,0:Ny+1), wc(0:Nxc+1,0:Nyc+1) !defects and errors etc.
		real(8) :: Auc(0:Nxc,1:Nyc), Avc(1:Nxc,0:Nyc), Azc(1:Nxc,1:Nyc)

		call calc_channel(Ny,Y,jmin,jmax)

		!Presmoothing
		do nt = 1, nu1
			call smooth(z,Au,Av,Az,b,Nx,Ny,jmin,jmax)
			! if(k==l) then
			! 	call boundary(z,Nx,Ny,jmin,jmax)
			! else
			! 	call boundary_defect(z,Nx,Ny)
			! end if

		end do
		! if (k==1) then
		! 	write(20,*) z(1:Nx,1:Ny)
		! end if
		! select case(k)
		! case(1)
		! 	write(23,*) z(1:Nx,1:Ny)
		! case(2)
		! 	write(24,*) z(1:Nx,1:Ny)
		! case(3)
		! 	write(25,*) z(1:Nx,1:Ny)
		! case default
		! end select

		!Coarse grid correction
		!Compute the defect
		do j = 1, Ny
			! i = 1
			! if (j<jmin .or. j>jmax) then
			! 	df(i,j) = b(i,j) + Au(i-1,j)*z(i-1,j) + Au(i,j)*z(i+1,j) + Av(i,j-1)*z(i,j-1) + Av(i,j)*z(i,j+1) - Az(i,j)*z(i,j)
			! end if
			do i = 1, Nx
				df(i,j) = b(i,j) + Au(i-1,j)*z(i-1,j) + Au(i,j)*z(i+1,j) + Av(i,j-1)*z(i,j-1) + Av(i,j)*z(i,j+1) - Az(i,j)*z(i,j)
			end do
		end do
		! if(k==2) then
		! 	write(20,*) Au(1:Nx,1:Ny)
		! end if

		!Restrice the defect
		dc(:,:) = 0.d0
		call Prolongation(Au,Av,Az,Auc,Avc,Azc,Nxc,Nyc)
		call Prolongation_defect(df,dc,Nxc,Nyc)
		call calc_channel(Nyc,Y,jminc,jmaxc)

		! if(k==1) then
			! write(20,*) Azc(1:Nxc,1:Nyc)
			! write(*,*) Auc(0:1,10), Avc(1,9:10), Azc(1,10), dc(1,10)
			! write(*,*) Azc(:,:)
		! 	do j = 1, Nyc
		! 		do i = 1, Nxc
		! 			write(*,*) i, j, (Auc(i,j)+Avc(i,j))/(Azc(i,j)-Auc(i-1,j)-Avc(i,j-1)) 
		! 			write(*,*) Auc(i-1:i,j), Avc(i,j-1:j), Azc(i,j)
		! 			write(*,*) Azc(i,j)
		! 		end do
		! 	end do
		! end if

		!Compute an approximate solution v of the defect equation on k-1
		wc(:,:) = 0.d0
		if(k==1) then
			do nt = 1, 1
				call smooth(wc,Auc,Avc,Azc,dc,Nxc,Nyc,jminc,jmaxc)
				! call boundary_defect(wc,Nxc,Nyc)
				! write(*,*) jminc,jmaxc
				! write(*,*) wc(1,9:11)
				! write(*,*) Auc(0:1,10), Avc(1,9:10), Azc(1,10), dc(1,10)
				! write(26,*) wc(1:Nxc,1:Nyc)
			end do
		else
			call MGCYC(l,k-1,wc,Auc,Avc,Azc,dc,nu1,nu2,Y,Nxc,Nyc,Nxc/2,Nyc/2,Res)
		end if

		!Interpolate the correction
		wf(:,:) = 0.d0
		call Interpolation_defect(wc,wf,Nx,Ny)

		!Compute the corrected approximation on k
		z(:,:) = z(:,:) + wf(:,:)

		!Postsmoothing
		do nt = 1, nu2
			call smooth(z,Au,Av,Az,b,Nx,Ny,jmin,jmax)
			! call boundary_defect(z,Nx,Ny)
			! if(k==l) then
			! 	call boundary(z,Nx,Ny,jmin,jmax)
			! else
			! 	call boundary_defect(z,Nx,Ny)
			! end if
		end do

		Res = 0.d0
		do j = 1, Ny
			do i = 1, Nx
				Res = Res + (b(i,j) + Au(i-1,j)*z(i-1,j) + Au(i,j)*z(i+1,j) + Av(i,j-1)*z(i,j-1) + Av(i,j)*z(i,j+1) - Az(i,j)*z(i,j))**2
			end do
		end do
		Res = Res**0.5d0

		! select case(k)
		! case(1)
		! 	write(20,*) z(1:Nx,1:Ny)
		! case(2)
		! 	write(21,*) z(1:Nx,1:Ny)
		! case(3)
		! 	write(22,*) z(1:Nx,1:Ny)
		! case default
		! end select

	end subroutine MGCYC

	subroutine smooth(z,Au,Av,Az,b,Nx,Ny,jmin,jmax)
		implicit none
		
		integer, intent(in) :: Nx, Ny, jmin, jmax !jmin and jmax is channel indice
		real(8), intent(in) :: Au(0:Nx,1:Ny), Av(1:Nx,0:Ny), Az(1:Nx,1:Ny), b(1:Nx,1:Ny)
		real(8), intent(inout) :: z(0:Nx+1,0:Ny+1)

		integer :: i, j

		do j = 1, Ny
		! do j = Ny, 1, -1
			! i = 1
			! if (j<jmin .or. j>jmax) then
			! 	z(i,j) = (b(i,j)+Au(i-1,j)*z(i-1,j)+Au(i,j)*z(i+1,j)+Av(i,j-1)*z(i,j-1)+Av(i,j)*z(i,j+1))/Az(i,j)
			! end if 
			do i = 1, Nx
			! do i = Nx, 1, -1
				z(i,j) = (b(i,j)+Au(i-1,j)*z(i-1,j)+Au(i,j)*z(i+1,j)+Av(i,j-1)*z(i,j-1)+Av(i,j)*z(i,j+1))/Az(i,j)
			end do
		end do
		! call boundary(z,Nx,Ny,jmin,jmax)

	end subroutine smooth

end program main